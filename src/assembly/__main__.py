from sys import exit
from . import assembly

if __name__ == "__main__":
    args = parse_args()
    exit(
        main(args.input, args.output, args.kmer_size, args.threshold, args.mismatches)
    )
