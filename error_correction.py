import sys
from typing import List, Dict

import matplotlib.pyplot as plt

import utils


def create_histogram(k_mers: Dict[str, int]) -> None:
    bins = max(k_mers.values())
    return plt.hist(k_mers.values(), bins)


def correct_read(
    read: str,
    k: int,
    kmerhist: Dict[str, int],
    threshold: int,
    mismatches: int = 1,
    alphabet: str = "ACTG",
) -> str:
    # read_cache = read
    # print(f"Attempting correction of read {read}")
    for i in range(len(read) - (k - 1)):
        k_mer = read[i : i + k]
        if kmerhist.get(k_mer, 0) <= threshold:
            for corrected_k_mer in utils.strings_within_hamming_distance(
                k_mer, mismatches, alphabet
            ):
                if kmerhist.get(corrected_k_mer, 0) > threshold:
                    read = read[:i] + corrected_k_mer + read[i + k :]
                    # print(f"k_mer: {k_mer} corrected to {corrected_k_mer}")
                    break
    # if read_cache == read:
    #     print(f"Read {read} didn't change")
    # else: 
    #     print(f"read {read_cache} corrected to {read}")
    return read


def drop_reads_below_threshold(
    data: List[str], k: int, kmerhist: Dict[str, int], threshold: int
) -> List[str]:
    filtered: List[str] = data[:]
    for read in data:
        for k_mer, _ in utils.chop(read, k):
            if kmerhist.get(k_mer, 0) <= threshold:
                filtered.remove(read)
                break
    return filtered


def correct_dataset(
    data: List[str],
    k: int,
    kmer_dist: Dict[str, int],
    threshold: int = 1,
    mismatches: int = 1,
    alphabet: str = "ACTG",
) -> List[str]:
    return [
        correct_read(read, k, kmer_dist, threshold, mismatches, alphabet)
        for read in data
    ]


def main() -> int:
    # counts: Dict[int, int] = {}
    # with open("./test_data", "r") as file:
    #     data = [read.strip("\n") for read in file.readlines()]
    # k_mers = utils.draw_k_mers(data, 20)
    # corrected_data = []
    # for read in data:
    #     corrected_data.append(correct(read, 20, k_mers, 1))
    # corrected_k_mers = utils.draw_k_mers(corrected_data, 20)
    # bins = max(corrected_k_mers.values())
    # fig, axs = plt.subplots(1, 2, figsize = (12, 5), sharey=True)
    # axs[0].hist(k_mers.values(), bins)
    # axs[0].set_xlabel("k-mer count")
    # axs[0].set_ylabel("# of distinct k-mers with that count")
    # axs[1].hist(corrected_k_mers.values(), bins)
    # axs[1].set_xlabel("k-mer count")
    # plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
