
import argparse
from ui.cleaner import clean_fasta_stream

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--min-len", type=int, default=200)
    p.add_argument("--max-n", type=float, default=0.05)
    p.add_argument("--no-dedup", action="store_true")
    args = p.parse_args()

    counters = clean_fasta_stream(
        args.inp, args.out_fasta, args.summary,
        min_len=args.min_len, max_N_fraction=args.max_n, dedup=not args.no_dedup
    )
    print("Done. Summary:", counters)

if __name__ == "__main__":
    main()
