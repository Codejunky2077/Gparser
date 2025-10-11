# backend/cleaner.py
from Bio import SeqIO
import hashlib
import sqlite3
import csv
import tempfile
import os
from typing import Optional, Tuple

ALLOWED_CHARS = set("ACGTN") 

def gc_pct(seq: str) -> float:
    s = seq.upper().replace('U','T')
    if len(s) == 0:
        return 0.0
    g = s.count('G')
    c = s.count('C')
    return 100.0 * (g + c) / len(s)

def normalize_seq(seq: str) -> str:
    return seq.upper().replace('U','T')

def clean_fasta(
    input_path: str,
    out_fasta: str,
    summary_csv: str,
    min_len: int = 200,
    max_N_fraction: float = 0.05,
    dedup: bool = True,
    dedup_db_path: Optional[str] = None,
) -> dict:
    """
    Streams through input FASTA, writes cleaned FASTA and a summary CSV.
    Returns counters summary.
    """
    # prepare dedup db
    if dedup:
        if dedup_db_path is None:
            dedup_db_path = tempfile.NamedTemporaryFile(delete=False, suffix='.db').name
        conn = sqlite3.connect(dedup_db_path)
        conn.execute("CREATE TABLE IF NOT EXISTS seq_hash(h TEXT PRIMARY KEY)")
        conn.commit()
    else:
        conn = None

    counters = {
        "total": 0, "kept": 0,
        "too_short": 0, "high_Ns": 0, "non_allowed": 0, "duplicate": 0
    }

    out_handle = open(out_fasta, "w")
    csvf = open(summary_csv, "w", newline="")
    writer = csv.writer(csvf)
    writer.writerow(["id","length","gc_pct","n_count","md5","status"])

    try:
        for rec in SeqIO.parse(input_path, "fasta"):
            counters["total"] += 1
            seq = normalize_seq(str(rec.seq))
            md5 = hashlib.md5(seq.encode()).hexdigest()
            # dedup check
            if conn is not None:
                try:
                    conn.execute("INSERT INTO seq_hash(h) VALUES(?)", (md5,))
                    conn.commit()
                except Exception as e:
                    # likely IntegrityError -> duplicate
                    writer.writerow([rec.id, len(seq), round(gc_pct(seq),2), seq.count('N'), md5, "duplicate"])
                    counters["duplicate"] += 1
                    continue

            # length filter
            if len(seq) < min_len:
                writer.writerow([rec.id, len(seq), round(gc_pct(seq),2), seq.count('N'), md5, "too_short"])
                counters["too_short"] += 1
                continue

            # N's filter
            if seq.count('N') / max(1, len(seq)) > max_N_fraction:
                writer.writerow([rec.id, len(seq), round(gc_pct(seq),2), seq.count('N'), md5, "high_Ns"])
                counters["high_Ns"] += 1
                continue

            # allowed characters (AGCTN only)
            if any(ch not in ALLOWED_CHARS for ch in set(seq)):
                writer.writerow([rec.id, len(seq), round(gc_pct(seq),2), seq.count('N'), md5, "non_allowed"])
                counters["non_allowed"] += 1
                continue

            # passed all checks -> write
            SeqIO.write(rec, out_handle, "fasta")
            writer.writerow([rec.id, len(seq), round(gc_pct(seq),2), seq.count('N'), md5, "kept"])
            counters["kept"] += 1

    finally:
        out_handle.close()
        csvf.close()
        if conn:
            conn.close()
        # try remove temp db if it was created automatically
        if dedup and dedup_db_path and os.path.exists(dedup_db_path):
            try:
                os.remove(dedup_db_path)
            except:
                pass

    return counters
