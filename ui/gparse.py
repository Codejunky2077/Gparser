# ui/app.py
import streamlit as st
import tempfile, os
from cleaner import clean_fasta_stream

st.set_page_config(page_title="Gparser🧬", layout="wide")
st.title("Gparser🧬 -FASTA Cleaner")

uploaded = st.file_uploader("Upload FASTA (fa / fasta / gz)", type=["fa","fasta","gz","fa.gz","txt"])
min_len = st.number_input("Minimum sequence length", value=200, step=10)
max_n = st.slider("Max fraction of Ns allowed", min_value=0.0, max_value=0.5, value=0.05, step=0.01)
dedup = st.checkbox("Enable deduplication", value=True)

if uploaded is not None:
    tmp_in = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    tmp_in.write(uploaded.getbuffer())
    tmp_in.flush()

    st.write("File saved for processing:", tmp_in.name)
    run = st.button("Start cleaning")

    if run:
        tmp_out = tempfile.NamedTemporaryFile(delete=False, suffix=".clean.fasta")
        tmp_csv = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
        with st.spinner("Cleaning... (this may take time for large files)"):
            counters = clean_fasta_stream(
                input_path=tmp_in.name,
                out_fasta=tmp_out.name,
                summary_csv=tmp_csv.name,
                min_len=int(min_len),
                max_N_fraction=float(max_n),
                dedup=dedup
            )
        st.success("Cleaning finished")
        st.write(counters)
        with open(tmp_out.name, "rb") as f:
            st.download_button("Download cleaned FASTA", data=f, file_name="cleaned.fasta")
        with open(tmp_csv.name, "rb") as f:
            st.download_button("Download summary CSV", data=f, file_name="summary.csv")

        # optional cleanup
        # os.remove(tmp_in.name)
