#!/usr/bin/env python

# This file is an entry point for the project
import argparse
import sys
from pathlib import Path

from genome_assembler import utils
from genome_assembler import graph
from genome_assembler import error_correction
from genome_assembler import progress

# Keep script name in global constant
SCRIPT_NAME = Path(__file__).name


def parse_args() -> argparse.Namespace:
    """Returns parsed command line arguments"""
    parser = argparse.ArgumentParser(
        description="A CLI genome assembler based on de Bruijn graph method written in Python."
    )
    parser.add_argument("input", type=str, help="input fasta file with reads")
    parser.add_argument("output", type=str, help="output fasta file with contigs")
    parser.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=21,
        help="size of k-mers in de Bruijn graph (default 21)",
    )
    parser.add_argument(
        "-t", "--threshold", type=int, default=2, help="threshold for data correction"
    )
    parser.add_argument(
        "-m",
        "--mismatches",
        type=int,
        default=2,
        help="mismatches for neighbour search during data correction",
    )
    return parser.parse_args()


def main(
    file: str, output: str, kmer_size: int = 21, threshold: int = 2, mismatches: int = 2
) -> int:
    K_MER_SUZE: int = kmer_size
    if any(arg < 1 for arg in (kmer_size, threshold, mismatches)):
        print(
            f"{SCRIPT_NAME}: error: kmer_size, threshold and mismatches args must all be larger than 0."
        )

    print("Parsing data...")
    data = utils.parse_fasta(file)
    k_mers = utils.draw_k_mers(data, K_MER_SUZE)
    print("Data correction... (may take long time)")
    data = error_correction.correct_dataset(
        data, K_MER_SUZE, k_mers, threshold, mismatches
    )
    data = error_correction.drop_reads_below_threshold(
        data, K_MER_SUZE, utils.draw_k_mers(data, K_MER_SUZE), 1
    )

    print("Graph construction...")
    de_bruijn_base = graph.DeBruijn.from_reads(data, K_MER_SUZE)
    de_bruijn_base = de_bruijn_base.trim_short_tips(K_MER_SUZE)
    de_bruijn_base = de_bruijn_base.trim_bubbles()
    sub_graphs = [
        graph.DeBruijn.from_graph(g) for g in de_bruijn_base.connected_components()
    ]
    for component in sub_graphs:
        component._update_tips()

    print("Creating contigs...")
    contigs = []
    graph_count = 0
    bar = progress.ProgressBar(len(sub_graphs))
    for i, g in enumerate(sub_graphs):
        graph_count += 1
        contigs.append(g.merge_nodes(K_MER_SUZE - 4, K_MER_SUZE - 2))
        bar.update(graph_count)
    bar.finish()

    print("Saving data...")
    contig_count = 0
    bar = progress.ProgressBar(len(contigs))
    with open(output, "w") as fasta:
        for i, contig in enumerate(contigs):
            contig_count += 1
            fasta.write(f">contig {i}\n")
            fasta.write(f"{contig[0]}\n")
            bar.update(contig_count)
        bar.finish()
    print("Script exited successfully!")
    return 0


if __name__ == "__main__":
    args = parse_args()
    sys.exit(
        main(args.input, args.output, args.kmer_size, args.threshold, args.mismatches)
    )
