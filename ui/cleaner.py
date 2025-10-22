from Bio import SeqIO
import hashlib
import sqlite3
import csv
import os
import gzip
from typing import Optional, Dict, TextIO

ALLOWED_CHARS = set("ACGTN")

def gc_pct(seq: str) -> float:
    """Return GC percentage of a sequence."""
    s = seq.upper().replace('U', 'T')
    if not s:
        return 0.0
    g = s.count('G')
    c = s.count('C')
    return round(100.0 * (g + c) / len(s), 2)

def normalize_seq(seq: str) -> str:
    """Uppercase and replace U with T."""
    return seq.upper().replace('U', 'T')

def open_fasta_maybe_gz(path: str) -> TextIO:
    """Open compressed or regular FASTA."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def clean_fasta(
    input_path: str,
    out_fasta: str,
    summary_csv: str,
    min_len: int = 200,
    max_N_fraction: float = 0.05,
    dedup: bool = True,
    dedup_db_path: Optional[str] = None
) -> Dict[str, int]:
    """
    Clean a FASTA file and write outputs.

    - Removes sequences that are too short, contain too many Ns, or invalid characters.
    - Deduplicates based on MD5 of normalized sequence.
    - Returns a summary counter dict.
    """

    # Dedup setup
    conn = None
    seen_hashes = set()

    if dedup:
        if dedup_db_path:
            conn = sqlite3.connect(dedup_db_path)
            conn.execute("CREATE TABLE IF NOT EXISTS seq_hash(h TEXT PRIMARY KEY)")
            conn.commit()

    counters = {
        "total": 0,
        "kept": 0,
        "too_short": 0,
        "high_Ns": 0,
        "non_allowed": 0,
        "duplicate": 0
    }

    def record_status(rec_id, seq, md5, status):
        writer.writerow([
            rec_id,
            len(seq),
            gc_pct(seq),
            seq.count("N"),
            md5,
            status
        ])

    with open(out_fasta, "w") as outf, open(summary_csv, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(["id", "length", "gc_pct", "n_count", "md5", "status"])

        with open_fasta_maybe_gz(input_path) as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                counters["total"] += 1
                seq = normalize_seq(str(rec.seq))
                md5 = hashlib.md5(seq.encode()).hexdigest()

                # Deduplication
                if dedup:
                    duplicate = False
                    if conn:
                        try:
                            conn.execute("INSERT INTO seq_hash(h) VALUES(?)", (md5,))
                            conn.commit()
                        except sqlite3.IntegrityError:
                            duplicate = True
                    else:
                        if md5 in seen_hashes:
                            duplicate = True
                        else:
                            seen_hashes.add(md5)

                    if duplicate:
                        record_status(rec.id, seq, md5, "duplicate")
                        counters["duplicate"] += 1
                        continue

                # Length filter
                if len(seq) < min_len:
                    record_status(rec.id, seq, md5, "too_short")
                    counters["too_short"] += 1
                    continue

                # N content filter
                if (seq.count("N") / len(seq)) > max_N_fraction:
                    record_status(rec.id, seq, md5, "high_Ns")
                    counters["high_Ns"] += 1
                    continue

                # Allowed character filter
                if not set(seq).issubset(ALLOWED_CHARS):
                    record_status(rec.id, seq, md5, "non_allowed")
                    counters["non_allowed"] += 1
                    continue

                 # Passed all checks
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord

                    # Create a new SeqRecord with normalized sequence and same record ID
                clean_rec = SeqRecord(Seq(seq), id=rec.id, description="")
                SeqIO.write(clean_rec, outf, "fasta")
                record_status(rec.id, seq, md5, "kept")
                counters["kept"] += 1

    # Cleanup
    if conn:
        conn.close()
        if dedup_db_path is None:
            try:
                os.remove(dedup_db_path)
            except:
                pass

    return counters
