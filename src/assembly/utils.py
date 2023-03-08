"""Set of utilities that may come in handy during the developement of the project."""
import random
import sys
import itertools
from typing import List, Generator, Tuple, Dict, Protocol


class Node(Protocol):
    k_minus_1_mer: str
    _indegree: int
    _outdegree: int

    def __hash__(self) -> int:
        ...


    @property
    def indegree(self) -> int:
        """Number of input egdes."""
        return self._indegree

    @property
    def outdegree(self) -> int:
        """Number of output edges"""
        return self._outdegree

    def is_balanced(self) -> bool:
        ...

    def is_semi_balanced(self) -> bool:
        ...


class NoOverlapError(ValueError):
    """Raised when two strings required to be overlapping are not."""

    pass


def random_sequence(length: int, symbols: str) -> str:
    """Returns random sequence of given length, consisting of given set of symbols.

    length - length of in characters of sequence to be generated.
    symbols - an alphabet over which the sequence will be generated.
    """
    return "".join([random.choice(symbols) for _ in range(length)])


def random_SN_error(sequence: str) -> str:
    position = random.randint(0, len(sequence))
    return (
        sequence[:position]
        + random.choice("ACTG".replace(sequence[position], ""))
        + sequence[position + 1 :]
    )


def reads_from_sequence(
    sequence: str, length: int, coverage: int, error: float = 0
) -> List[str]:
    """Returns a list of random reads of given length sampled at random from given sequence.

    sequence - a sequence that will be sampled.
    length - length of individual reads that will be generated.
    coverage - a number representing how many times the sequence will be covered with reads.
        Its purely stochastic, as reads are generated fully randomly, and this paramter
        is used only to calculate how many reads to generate.
    """
    required_reads: int = int(coverage * len(sequence) / length)
    reads: List[str] = []
    for _ in range(required_reads):
        start: int = random.randint(0, len(sequence) - length)
        end: int = start + length
        read = sequence[start:end]
        if random.random() < error:
            read = random_SN_error(read)
        reads.append(read)
    return reads


def are_overlapped(seq1: str, seq2: str, min_overlap: int = 1) -> bool:
    """Returns True if suffix of seq1 is equal to prefix of seq2 -- if seq1 and seq2 are ovelapped.

    seq1 - string whose suffix will be checked against the prefix of seq2.
    seq2 - string whose prefix will be checked against the suffix of seq1.
    min_overlap - minimal length by which seq1 and seq2 are required to be overlapped."""
    suffix_of_1_in_2_position: int = seq2.rfind(seq1[-min_overlap:])
    prefix_length: int = suffix_of_1_in_2_position + min_overlap
    return seq2[:prefix_length] == seq1[-prefix_length:]


def overlapped_by(seq1: str, seq2: str, min_overlap: int) -> int:
    """Returns a length by which seq1 and seq2 are overlapped."""
    suffix_of_1_in_2_position: int = seq2.rfind(seq1[-min_overlap:])
    prefix_length: int = suffix_of_1_in_2_position + min_overlap
    print(seq2[:prefix_length])
    print(seq1[-prefix_length:])
    if seq2[:prefix_length] == seq1[-prefix_length:]:
        return prefix_length
    return -1


def concat_overlapped(string1: str, string2: str) -> str:
    """Returns overlap sensitive concatenation of string1 and string2.

    strin1, string2 - strings to be concatenated.

    Raises NoOverlapError if strings are not overlapping
    """
    if (overlap := overlapped_by(string1, string2, min_overlap=1)) > 0:
        return string1 + string2[overlap:]
    elif (overlap := overlapped_by(string2, string1, min_overlap=1)) > 0:
        return string2 + string1[overlap:]
    else:
        print(overlap)
        raise NoOverlapError


def chop(string: str, k: int) -> Generator[Tuple[str, int], None, None]:
    for i in range(len(string) - (k - 1)):
        yield string[i : i + k], i


def hamming_distance(x: str, y: str) -> int:
    if len(x) != len(y):
        raise ValueError(
            "Hamming distance is undefined for sequences of uneven length!"
        )
    return sum(xi != yi for xi, yi in zip(x, y))


def strings_within_hamming_distance(
    string: str, distance: int, alphabet: str = "ACTG"
) -> Generator[str, None, None]:
    for positions in itertools.combinations(range(len(string)), distance):
        for substitutions in itertools.product(alphabet, repeat=distance):
            neighbour: list[str] = list(string)
            for i, char in zip(positions, substitutions):
                if string[i] != char:
                    neighbour[i] = char
            yield "".join(neighbour)


def neighbours1mm(k_mer: str, alphabet: str = "ACTG") -> List[str]:
    """Generate all neighbours at Hamming distance 1 frmo k_mer"""
    neighbours: List[str] = []
    for i in range(len(k_mer) - 1, -1, -1):
        old_char = k_mer[i]
        for char in alphabet:
            if char == old_char: continue
            neighbours.append(k_mer[:i] + char + k_mer[i + 1:])
    return neighbours


def draw_k_mers(reads: List[str], k) -> Dict[str, int]:
    k_mers: Dict[str, int] = {}
    for read in reads:
        for k_mer, _ in chop(read, k):
            k_mers[k_mer] = k_mers.get(k_mer, 0) + 1
    return k_mers


def dfs(G: Dict[Node, List[Node]], node: Node):
    discovered: List[Node] = []

    def _travel(node: Node):
        discovered.append(node)
        for neighbour in G[node]:
            if neighbour not in discovered:
                _travel(neighbour)

    _travel(node)
    return discovered

def parse_fasta(path_to_file: str) -> List[str]:
    sequences: List[str] = []
    with open(path_to_file, "r") as file:
        for line in file.readlines():
            if not line.startswith(">"):
                sequences.append(line.strip("\n"))
    return sequences


def generate_all_k_mers(length: int, alphabet: str = "ACTG") -> List[str]:
    k_mers: List[str] = []
    for combination in itertools.product(alphabet, repeat=length):
        k_mers.append("".join(combination))
    return k_mers


def main() -> int:
    """MAIN FUNCTIOM IN THIS FILE IS USED ONLY FOR DEBUGGING."""
    # Generate random sequence
    random_seq = random_sequence(55, "ACTG")
    # Generate all hexamers
    all_6_mers = generate_all_k_mers(6)
    # Print how many poosible hexamers there are
    print(len(all_6_mers))
    # Print created random sequence
    print(random_seq)
    # Print hexmaers from random sequence, and then how many of them there is
    print(draw_k_mers([random_seq], 6))
    print(len(draw_k_mers([random_seq], 6)))
    return 0


# Run main only if this file is run directly and not if it's imported by other module
if __name__ == "__main__":
    sys.exit(main())
