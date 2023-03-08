"""File used for initial look at training data."""
import sys
from typing import Dict, List

import matplotlib.pyplot as plt

from genome_assembler import utils, error_correction

def parse_fasta(path_to_file: str) -> List[str]:
    sequences: List[str] = []
    with open(path_to_file, "r") as file:
        for line in file.readlines():
            if not line.startswith(">"):
                sequences.append(line.strip("\n"))
    return sequences


def main() -> int:
    DATA_PATH = "PATH TO DATA FILE"
    # data: List[str] = []
    # for i in range(1, 4):
    #     data.extend(parse_fasta(f"{DATA_PATH}/reads{i}.fasta"))
    data = parse_fasta(f"{DATA_PATH}/reads1.fasta")
    k_mers = utils.draw_k_mers(data, 45)
    corrected_data = error_correction.correct_dataset(data, 45, k_mers, 1, 1)
    corrected_k_mers = utils.draw_k_mers(corrected_data, 45)
    corrected_data = error_correction.drop_reads_below_threshold(corrected_data, 45, corrected_k_mers, 1)
    corrected_k_mers = utils.draw_k_mers(corrected_data, 45)
    print(bins := max(k_mers.values()))
    print(bins_corrected := max(corrected_k_mers.values()))

    _, axs = plt.subplots(1, 2, figsize = (12, 5))
    axs[0].hist(k_mers.values(), bins)
    axs[0].set_xlabel("k-mer count")
    axs[0].set_ylabel("# of distinct k-mers with that count")
    axs[0].set_yscale("log")
    axs[1].hist(corrected_k_mers.values(), bins_corrected)
    axs[1].set_xlabel("k-mer count")
    axs[1].set_yscale("log")
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
