#!/usr/bin/env python
import argparse
import sys

import utils
import graph
import error_correction
import progress

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Input fasta file with reads")
    parser.add_argument("output", type=str, help="output fasta file with contigs")
    return parser.parse_args()

def main(file: str, output: str) -> int:
    K_MER_SUZE: int = 20

    print("Parsing data")
    data = utils.parse_fasta(file)
    k_mers = utils.draw_k_mers(data, K_MER_SUZE)
    print("Data correction")
    data = error_correction.correct_dataset(data, K_MER_SUZE, k_mers, 2, 2)
    data = error_correction.drop_reads_below_threshold(
        data, K_MER_SUZE, utils.draw_k_mers(data, K_MER_SUZE), 1
    )
    
    print("Graph construction")
    de_bruijn_base = graph.DeBruijn.from_reads(data, 20)
    de_bruijn_base = de_bruijn_base.trim_short_tips(17)
    de_bruijn_base = de_bruijn_base.trim_bubbles()
    sub_graphs = [graph.DeBruijn.from_graph(g) for g in de_bruijn_base.connected_components()]
    for component in sub_graphs:
        component._update_tips()
    
    print("Creating contigs")
    contigs = []
    graph_count = 0
    bar = progress.ProgressBar(len(sub_graphs))
    for i, g in enumerate(sub_graphs):
        graph_count += 1
        contigs.append(g.merge_nodes(13, 15))
        bar.update(graph_count)
    bar.finish()

    print("Saving data")
    contig_count = 0
    bar = progress.ProgressBar(len(contigs))
    with open(output, "w") as fasta:
        for i, contig in enumerate(contigs):
            contig_count += 1
            fasta.write(f">contig {i}\n")
            fasta.write(F"{contig[0]}\n")
            bar.update(contig_count)
        bar.finish()
    return 0


if __name__ == "__main__":
    args = parse_args()
    sys.exit(main(args.input, args.output))
