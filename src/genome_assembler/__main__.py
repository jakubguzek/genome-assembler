from sys import exit
from genome_assembler.assembly import main, parse_args

if __name__ == "__main__":
    args = parse_args()
    exit(
        main(args.input, args.output, args.kmer_size, args.threshold, args.mismatches)
    )
