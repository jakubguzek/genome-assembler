import sys

import utils
import graph
import error_correction
import greedy_scs
import progress


def main(file: str) -> int:
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
    sub_graphs = [graph.DeBruijn.from_graph(g) for g in de_bruijn_base.connected_components()]
    
    print("Creating contigs")
    contigs = []
    graph_count = 0
    bar = progress.ProgressBar(len(sub_graphs))
    for g in sub_graphs:
        graph_count += 1
        reads_from_graph = [k_minus_1_mer for k_minus_1_mer in g.nodes]
        contigs.append(greedy_scs.greedy_scs(reads_from_graph, K_MER_SUZE-22))
        bar.update(graph_count)
    bar.finish()

    print("Saving data")
    contig_count = 0
    bar = progress.ProgressBar(len(contigs))
    with open("contigs.fasta", "w") as fasta:
        for i, contig in enumerate(contigs):
            contig_count += 1
            fasta.write(f">contig {i}\n")
            fasta.write(F"{contig}\n")
            bar.update(contig_count)
        bar.finish()

    return 0


if __name__ == "__main__":
    sys.exit(main("/home/jakub/Documents/studia/technologie_w_skali_genomowej/projekty/2/training/reads/reads1.fasta"))
