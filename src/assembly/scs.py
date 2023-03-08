import itertools
from typing import Tuple, List

def overlap(a, b, min_length=3):
    """Return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0."""
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous matcho


def maximal_overlap(strings: List[str], k: int, max_length: int) -> Tuple[str, str, int]:
    string_1, string_2 = None, None
    best = 0
    for a, b in itertools.permutations(strings, 2):
        overlap_length = overlap(a, b, min_length=k)
        if overlap_length >= best:
            string_1, string_2 = a, b
            best = overlap_length
        if best >= max_length:
            break
    return string_1, string_2, best


def greedy_scs(reads, k, max_length):
    """Greedy shortest-common-superstring merge.
    Repeat until no edges (overlaps of length >= k)
    remain."""
    read_a, read_b, olen = maximal_overlap(reads, k, max_length)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[-(len(read_b) - olen) :])
        read_a, read_b, olen = maximal_overlap(reads, k, max_length)
    return "".join(reads[0])
