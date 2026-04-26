#!/usr/bin/env python3
"""Run SaProt mutational-effect scoring followed by GPSite binding prediction."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


GPSITE_OVERVIEW_COLUMNS = [
    "DNA_Binding",
    "RNA_Binding",
    "Peptide_Binding",
    "Protein_Binding",
    "ATP_Binding",
    "HEM_Binding",
    "ZN_Binding",
    "CA_Binding",
    "MG_Binding",
    "MN_Binding",
]


def load_config(path: str | None) -> dict:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def cfg_value(args: argparse.Namespace, config: dict, name: str, default=None):
    value = getattr(args, name)
    if value is not None:
        return value
    return config.get(name, default)


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(chunks).upper()
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        records[current_id] = "".join(chunks).upper()

    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    return records


def write_fasta(records: dict[str, str], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def ensure_pdbs(records: dict[str, str], pdb_dir: Path) -> None:
    missing = []
    for protein_id in records:
        if not (pdb_dir / f"{protein_id}.pdb").exists():
            missing.append(f"{protein_id}.pdb")
    if missing:
        raise FileNotFoundError(
            "Missing PDB files matching FASTA headers: " + ", ".join(missing[:10])
        )


def run_command(command: list[str], cwd: Path | None = None) -> None:
    print("\n[run]", " ".join(command), flush=True)
    subprocess.run(command, cwd=str(cwd) if cwd else None, check=True)


def gpsite_output_dir(gpsite_out_root: Path, fasta_path: Path) -> Path:
    run_id = fasta_path.name.rsplit(".", 1)[0].replace(" ", "_")
    return gpsite_out_root / run_id


def read_gpsite_overview(path: Path, prefix: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"GPSite overview not found: {path}")
    rows: dict[str, dict[str, str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["ID"]] = {
                f"{prefix}_{column}": row.get(column, "")
                for column in GPSITE_OVERVIEW_COLUMNS
            }
            rows[row["ID"]][f"{prefix}_Length"] = row.get("Length", "")
            rows[row["ID"]][f"{prefix}_pLDDT"] = row.get("pLDDT", "")
            rows[row["ID"]][f"{prefix}_pTM"] = row.get("pTM", "")
    return rows


def select_mutants(
    saprot_csv: Path,
    fasta_records: dict[str, str],
    top_n: int,
    mutant_fasta: Path,
    selected_csv: Path,
) -> list[dict[str, str]]:
    with open(saprot_csv, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows or top_n <= 0:
        write_fasta({}, mutant_fasta)
        return []

    rows.sort(key=lambda item: float(item["saprot_score"]))
    selected = []
    for selection_group, group_rows in (
        ("low_score", rows[:top_n]),
        ("high_score", rows[-top_n:][::-1]),
    ):
        for rank, row in enumerate(group_rows, start=1):
            protein_id = row["protein_id"]
            sequence = list(fasta_records[protein_id])
            position = int(row["position"])
            mut_aa = row["mut_aa"]
            sequence[position - 1] = mut_aa
            mutant_id = f"{protein_id}_{row['mutation']}_{selection_group}"
            selected_row = {
                **row,
                "selection_group": selection_group,
                "selection_rank": str(rank),
                "mutant_id": mutant_id,
                "gpsite_id": mutant_id,
                "mutant_sequence": "".join(sequence),
            }
            selected.append(selected_row)

    write_fasta({row["mutant_id"]: row["mutant_sequence"] for row in selected}, mutant_fasta)
    with open(selected_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0].keys()))
        writer.writeheader()
        writer.writerows(selected)
    return selected


def merge_results(
    saprot_sorted_csv: Path,
    original_overview: Path,
    mutant_overview: Path | None,
    selected_mutants: list[dict[str, str]],
    final_csv: Path,
) -> None:
    original = read_gpsite_overview(original_overview, "original")
    mutant = read_gpsite_overview(mutant_overview, "mutant") if mutant_overview else {}
    selected_by_mutation = {
        (row["protein_id"], row["mutation"], row["selection_group"]): row
        for row in selected_mutants
    }

    with open(saprot_sorted_csv, "r", encoding="utf-8") as handle:
        saprot_rows = list(csv.DictReader(handle))

    merged_rows = []
    for rank, row in enumerate(saprot_rows, start=1):
        base = {
            **row,
            "saprot_rank": str(rank),
            **original.get(row["protein_id"], {}),
        }
        matching_selected = [
            selected
            for key, selected in selected_by_mutation.items()
            if key[0] == row["protein_id"] and key[1] == row["mutation"]
        ]
        if matching_selected:
            for selected in matching_selected:
                mutant_values = mutant.get(selected["gpsite_id"], {})
                delta_values = {}
                for column in GPSITE_OVERVIEW_COLUMNS:
                    original_key = f"original_{column}"
                    mutant_key = f"mutant_{column}"
                    delta_key = f"delta_{column}"
                    if original_key in base and mutant_key in mutant_values:
                        delta_values[delta_key] = (
                            float(mutant_values[mutant_key]) - float(base[original_key])
                        )
                merged_rows.append(
                    {
                        **base,
                        "selection_group": selected["selection_group"],
                        "selection_rank": selected["selection_rank"],
                        "gpsite_mutant_id": selected["gpsite_id"],
                        **mutant_values,
                        **delta_values,
                    }
                )
        else:
            merged_rows.append(base)

    fieldnames: list[str] = []
    for row in merged_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with open(final_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", help="Optional JSON config file.")
    parser.add_argument("--fasta", help="Input FASTA. Headers must match PDB basenames.")
    parser.add_argument("--pdb-dir", help="Directory containing <FASTA_ID>.pdb files.")
    parser.add_argument("--outdir", help="Output directory.")
    parser.add_argument("--saprot-repo", help="Path to SaProt repository.")
    parser.add_argument("--gpsite-repo", help="Path to GPSite repository.")
    parser.add_argument("--saprot-env", help="Conda environment name for SaProt.")
    parser.add_argument("--gpsite-env", help="Conda environment name for GPSite.")
    parser.add_argument("--conda-exe", help="Conda executable.")
    parser.add_argument("--saprot-model", help="Directory of SaProt_650M_AF2 model.")
    parser.add_argument("--saprot-device", help="Device for SaProt inference: cuda or cpu.")
    parser.add_argument("--foldseek", help="Foldseek binary path used by SaProt.")
    parser.add_argument("--chain", help="Single chain ID in each PDB.")
    parser.add_argument("--plddt-mask", action="store_true", default=None)
    parser.add_argument("--no-plddt-mask", action="store_false", dest="plddt_mask")
    parser.add_argument("--gpu", help="GPU id for GPSite. Omit for CPU.")
    parser.add_argument("--batch", type=int, help="GPSite batch size.")
    parser.add_argument("--prot-t5-path", help="Optional ProtT5 model path passed to GPSite.")
    parser.add_argument("--mutant-top-n", type=int, help="Select N low-score and N high-score mutants.")
    parser.add_argument("--skip-mutant-gpsite", action="store_true", default=None)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = load_config(args.config)

    fasta = Path(cfg_value(args, config, "fasta")).resolve()
    pdb_dir = Path(cfg_value(args, config, "pdb_dir")).resolve()
    outdir = Path(cfg_value(args, config, "outdir", "outputs")).resolve()
    saprot_repo = Path(cfg_value(args, config, "saprot_repo", "../SaProt-main")).resolve()
    gpsite_repo = Path(cfg_value(args, config, "gpsite_repo", "../GPSite-main")).resolve()
    saprot_env = cfg_value(args, config, "saprot_env", "saprot")
    gpsite_env = cfg_value(args, config, "gpsite_env", "GPSite")
    conda_exe = cfg_value(args, config, "conda_exe", "conda")
    saprot_model = Path(cfg_value(args, config, "saprot_model", "SaProt_650M_AF2")).resolve()
    saprot_device = cfg_value(args, config, "saprot_device", "cuda")
    foldseek = Path(cfg_value(args, config, "foldseek", saprot_repo / "bin" / "foldseek")).resolve()
    chain = cfg_value(args, config, "chain", "A")
    plddt_mask = bool(cfg_value(args, config, "plddt_mask", False))
    gpu = cfg_value(args, config, "gpu", None)
    batch = int(cfg_value(args, config, "batch", 4))
    prot_t5_path = cfg_value(args, config, "prot_t5_path", None)
    mutant_top_n = int(cfg_value(args, config, "mutant_top_n", 10))
    skip_mutant_gpsite = bool(cfg_value(args, config, "skip_mutant_gpsite", False))

    outdir.mkdir(parents=True, exist_ok=True)
    records = parse_fasta(fasta)
    ensure_pdbs(records, pdb_dir)

    saprot_csv = outdir / "saprot_mut_effect.csv"
    saprot_sorted_csv = outdir / "saprot_mut_effect.sorted.csv"
    selected_csv = outdir / "selected_mutants.csv"
    mutant_fasta = outdir / "selected_mutants.fasta"
    gpsite_original_fasta = outdir / "gpsite_original_input.fasta"
    final_csv = outdir / "final_merged.csv"
    write_fasta(records, gpsite_original_fasta)

    saprot_script = Path(__file__).resolve().parent / "saprot_mut_effect.py"
    saprot_command = [
        conda_exe,
        "run",
        "-n",
        saprot_env,
        "python",
        str(saprot_script),
        "--fasta",
        str(fasta),
        "--pdb-dir",
        str(pdb_dir),
        "--out-csv",
        str(saprot_csv),
        "--sorted-csv",
        str(saprot_sorted_csv),
        "--saprot-repo",
        str(saprot_repo),
        "--saprot-model",
        str(saprot_model),
        "--foldseek",
        str(foldseek),
        "--chain",
        chain,
        "--device",
        saprot_device,
    ]
    if plddt_mask:
        saprot_command.append("--plddt-mask")
    run_command(saprot_command)

    selected_mutants = select_mutants(
        saprot_sorted_csv, records, mutant_top_n, mutant_fasta, selected_csv
    )

    gpsite_original_root = outdir / "gpsite_original"
    gpsite_command = [
        conda_exe,
        "run",
        "-n",
        gpsite_env,
        "python",
        str(gpsite_repo / "script" / "predict.py"),
        "-i",
        str(gpsite_original_fasta),
        "-o",
        str(gpsite_original_root),
        "-b",
        str(batch),
    ]
    if gpu:
        gpsite_command.extend(["--gpu", str(gpu)])
    if prot_t5_path:
        gpsite_command.extend(["--prot-t5-path", str(Path(prot_t5_path).resolve())])
    run_command(gpsite_command, cwd=gpsite_repo)
    original_overview = gpsite_output_dir(gpsite_original_root, gpsite_original_fasta) / "pred" / "overview.txt"
    if not original_overview.exists():
        raise RuntimeError(
            "GPSite did not produce overview.txt for the original FASTA. "
            "Check the GPSite log above for input-format, GPU, or dependency errors."
        )

    mutant_overview = None
    if selected_mutants and not skip_mutant_gpsite:
        gpsite_mutant_root = outdir / "gpsite_mutants"
        mutant_command = [
            conda_exe,
            "run",
            "-n",
            gpsite_env,
            "python",
            str(gpsite_repo / "script" / "predict.py"),
            "-i",
            str(mutant_fasta),
            "-o",
            str(gpsite_mutant_root),
            "-b",
            str(batch),
        ]
        if gpu:
            mutant_command.extend(["--gpu", str(gpu)])
        if prot_t5_path:
            mutant_command.extend(["--prot-t5-path", str(Path(prot_t5_path).resolve())])
        run_command(mutant_command, cwd=gpsite_repo)
        mutant_overview = gpsite_output_dir(gpsite_mutant_root, mutant_fasta) / "pred" / "overview.txt"
        if not mutant_overview.exists():
            raise RuntimeError(
                "GPSite did not produce overview.txt for mutant FASTA. "
                "Check the GPSite log above for input-format, GPU, or dependency errors."
            )

    merge_results(
        saprot_sorted_csv=saprot_sorted_csv,
        original_overview=original_overview,
        mutant_overview=mutant_overview,
        selected_mutants=selected_mutants if not skip_mutant_gpsite else [],
        final_csv=final_csv,
    )
    print(f"\nDone. Final table: {final_csv}")


if __name__ == "__main__":
    main()
