#getting necessary imports for the cleaning function
from Bio import SeqIO
import streamlit as st
import hashlib
import sqlite3
import csv
import os
import gzip
from typing import Optional, Dict, TextIO

#only allowed characters in sequences given by user
ALLOWED_CHARS = set("ACGTN")

#function to calculate GC percentage
def gc_pct(seq: str) -> float:
    """Return GC percentage of a sequence."""
    s = seq.upper().replace('U', 'T')
    if not s:
        return 0.0
    g = s.count('G')
    c = s.count('C')
    return round(100.0 * (g + c) / len(s), 2)

#function to normalize sequences into DNA form
def normalize_seq(seq: str) -> str:
    """Convert RNA to DNA and uppercase."""
    return seq.upper().replace('U', 'T')

#function to open fasta files, whether gzipped(genebank format) or not.

def open_fasta_maybe_gz(path: str) -> TextIO:
    """Open compressed or regular FASTA."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

#main cleaning function
def clean_fasta(
    input_path: str,
    out_fasta: str,
    summary_csv: str,
    min_len: int = 200,
    max_N_fraction: float = 0.05,
    dedup: bool = True,
    dedup_db_path: Optional[str] = None,
    progress_callback=None  # Added callback parameter
) -> Dict[str, int]:
    # Setup connections and deduplication data structures as before
    conn = None
    seen_hashes = set()
    if dedup:
        if dedup_db_path:
            conn = sqlite3.connect(dedup_db_path)
            conn.execute("CREATE TABLE IF NOT EXISTS seq_hash(h TEXT PRIMARY KEY)")
            conn.commit()
            #setting the counters variables to zero
    total = kept = too_short = high_Ns = non_allowed = duplicate = 0

    #function to record status of each sequence in the summary csv
    def record_status(rec_id, seq, md5, status):
        writer.writerow([
            rec_id,
            len(seq),
            gc_pct(seq),
            seq.count("N"),
            md5,
            status
        ])
        # Flush to ensure data is written immediately reducing risk of data loss
    with open(out_fasta, "w") as outf, open(summary_csv, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(["id", "length", "gc_pct", "n_count", "md5", "status"])
        #processing each sequence record in the input fasta file
        with open_fasta_maybe_gz(input_path) as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                total += 1
                seq = normalize_seq(str(rec.seq))
                md5 = hashlib.md5(seq.encode()).hexdigest()

                # Deduplication check
                is_duplicate = False
                if dedup:
                    if conn:
                        try:
                            conn.execute("INSERT INTO seq_hash(h) VALUES(?)", (md5,))
                            conn.commit()
                        except sqlite3.IntegrityError:
                            is_duplicate = True
                    else:
                        if md5 in seen_hashes:
                            is_duplicate = True
                        else:
                            seen_hashes.add(md5)

                if is_duplicate:
                    record_status(rec.id, seq, md5, "duplicate")
                    duplicate += 1
                    if progress_callback:
                        progress_callback({
                            "total": total,
                            "kept": kept,
                            "too_short": too_short,
                            "high_Ns": high_Ns,
                            "non_allowed": non_allowed,
                            "duplicate": duplicate
                        })
                    continue

                # Length filter
                if len(seq) < min_len:
                    record_status(rec.id, seq, md5, "too_short")
                    too_short += 1
                    if progress_callback:
                        progress_callback({
                            "total": total,
                            "kept": kept,
                            "too_short": too_short,
                            "high_Ns": high_Ns,
                            "non_allowed": non_allowed,
                            "duplicate": duplicate
                        })
                    continue
                
                # N content filter
                if (seq.count("N") / len(seq)) > max_N_fraction:
                    record_status(rec.id, seq, md5, "high_Ns")
                    high_Ns += 1
                    if progress_callback:
                        progress_callback({
                            "total": total,
                            "kept": kept,
                            "too_short": too_short,
                            "high_Ns": high_Ns,
                            "non_allowed": non_allowed,
                            "duplicate": duplicate
                        })
                    continue
                
                # Allowed character filter
                if not set(seq).issubset(ALLOWED_CHARS):
                    record_status(rec.id, seq, md5, "non_allowed")
                    non_allowed += 1
                    if progress_callback:
                        progress_callback({
                            "total": total,
                            "kept": kept,
                            "too_short": too_short,
                            "high_Ns": high_Ns,
                            "non_allowed": non_allowed,
                            "duplicate": duplicate
                        })
                    continue

                # Passed all filters; write output
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                clean_rec = SeqRecord(Seq(seq), id=rec.id, description="")
                SeqIO.write(clean_rec, outf, "fasta")
                record_status(rec.id, seq, md5, "kept")
                kept += 1

                # Update live metrics after each record
                if progress_callback:
                    progress_callback({
                        "total": total,
                        "kept": kept,
                        "too_short": too_short,
                        "high_Ns": high_Ns,
                        "non_allowed": non_allowed,
                        "duplicate": duplicate
                    })
    # Cleanup
    if conn:
        conn.close()
    if dedup_db_path is None:
        try:
            os.remove(dedup_db_path)
        except:
            pass
    #returning the counters dictionary to show status in the main app ui
    return {
        "total": total,
        "kept": kept,
        "too_short": too_short,
        "high_Ns": high_Ns,
        "non_allowed": non_allowed,
        "duplicate": duplicate
    }

