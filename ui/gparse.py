# ui/app.py
import streamlit as st
import tempfile
from PIL import Image
from cleaner import clean_fasta
#page tab setup
st.set_page_config(page_title="Gparser🧬",page_icon=Image.open("dna.png"),layout="centered")
#page UI cosmetic config
st.markdown(f"""
            <img src="https://github.com/Codejunky2077/Gparser/blob/022317ba401641adf54177ab5a586835f689afd9/ui/dna.png?raw=true" width="100" style="margin-left:150px;margin-right:10px;margin-bottom:45px; vertical-align: middle;" />
            <span style="font-style:Fjalla One;font-size: 80px; font-weight: 300;">Gparser</span>
            <br><br>
            """,unsafe_allow_html=True)





#page ui procedure config
uploaded = st.file_uploader("Upload FASTA (fa / fasta / gz)", type=["fa","fasta","gz","fa.gz","txt"],help="Upload your GENOME file here for cleaning.")
min_len = st.number_input("Minimum sequence length", value=100, step=10,help="Sequences shorter than this will be removed")
max_n = st.slider("Max fraction of Ns allowed", min_value=0.0, max_value=0.5, value=0.05, step=0.01,help="Sequences with a higher fraction of 'N' bases will be removed")
dedup = st.checkbox("Enable deduplication", value=True,help="Remove duplicate sequences from the FASTA")

if uploaded is not None:
    tmp_in = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    tmp_in.write(uploaded.getbuffer())
    tmp_in.flush()

    run = st.button("Start cleaning")

    if run:
        tmp_out = tempfile.NamedTemporaryFile(delete=False, suffix=".clean.fasta")
        tmp_csv = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
        with st.spinner("Cleaning... please wait..."):
            counters = clean_fasta(
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


