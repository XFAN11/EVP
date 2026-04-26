#!/usr/bin/env python3
"""Score all possible single amino-acid substitutions with SaProt."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


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
    return records


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--pdb-dir", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--sorted-csv", required=True)
    parser.add_argument("--saprot-repo", required=True, help="Path to SaProt-main.")
    parser.add_argument(
        "--saprot-model",
        required=True,
        help="Path to SaProt_650M_AF2 directory. Edit this in config.json for your machine.",
    )
    parser.add_argument(
        "--foldseek",
        required=True,
        help="Path to foldseek binary. Usually SaProt-main/bin/foldseek.",
    )
    parser.add_argument("--chain", default="A", help="PDB chain ID. Single-chain PDBs usually use A.")
    parser.add_argument("--plddt-mask", action="store_true")
    parser.add_argument("--device", default="cuda", help="cuda or cpu for SaProt model inference.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    saprot_repo = Path(args.saprot_repo).resolve()
    sys.path.insert(0, str(saprot_repo))

    from model.saprot.saprot_foldseek_mutation_model import SaprotFoldseekMutationModel
    from utils.foldseek_util import get_struc_seq

    fasta_records = parse_fasta(Path(args.fasta))
    pdb_dir = Path(args.pdb_dir)

    model = SaprotFoldseekMutationModel(
        foldseek_path=None,
        config_path=str(Path(args.saprot_model).resolve()),
        load_pretrained=True,
    )
    model.eval()
    model.to(args.device)

    rows = []
    for protein_id, fasta_seq in fasta_records.items():
        pdb_path = pdb_dir / f"{protein_id}.pdb"
        parsed = get_struc_seq(
            str(Path(args.foldseek).resolve()),
            str(pdb_path),
            [args.chain],
            plddt_mask=args.plddt_mask,
        )[args.chain]
        structure_seq, _foldseek_seq, combined_seq = parsed

        if len(structure_seq) != len(fasta_seq):
            raise ValueError(
                f"Sequence length mismatch for {protein_id}: "
                f"FASTA={len(fasta_seq)}, PDB/Foldseek={len(structure_seq)}"
            )

        for position, wt_aa in enumerate(fasta_seq, start=1):
            if wt_aa not in AMINO_ACIDS:
                continue
            position_scores = model.predict_pos_mut(combined_seq, position)
            for mut_aa in AMINO_ACIDS:
                if mut_aa == wt_aa:
                    continue
                mutation = f"{wt_aa}{position}{mut_aa}"
                score = position_scores[mutation]
                rows.append(
                    {
                        "protein_id": protein_id,
                        "position": position,
                        "wt_aa": wt_aa,
                        "mut_aa": mut_aa,
                        "mutation": mutation,
                        "saprot_score": score,
                    }
                )

    out_csv = Path(args.out_csv)
    sorted_csv = Path(args.sorted_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["protein_id", "position", "wt_aa", "mut_aa", "mutation", "saprot_score"]

    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    rows.sort(key=lambda row: float(row["saprot_score"]))
    with open(sorted_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
