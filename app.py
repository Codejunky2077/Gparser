# ui/app.py
import streamlit as st
import tempfile
from PIL import Image
from cleaner_init_ import clean_fasta
#page tab setup
st.set_page_config(page_title="Gparser",page_icon=Image.open("dna.png"),layout="centered")
#page UI cosmetic config
st.markdown(f"""
            <img src="https://github.com/Codejunky2077/Gparser/blob/022317ba401641adf54177ab5a586835f689afd9/ui/dna.png?raw=true" width="100" style="margin-left:150px;margin-right:10px;margin-bottom:45px; vertical-align: middle;" />
            <span style="font-style:Fjalla One;font-size: 80px; font-weight: 300;">Gparser</span>
            <br>
            <span style="font-size:20px; font-weight: 200;text-align:centre;">A simple and effective genome parser saves your time and sanity brilliantly.</span>
            <br><br>
            """,unsafe_allow_html=True)





#page ui procedure config
uploaded = st.file_uploader("Upload File for parsing", type=["fa","fasta","gz","fa.gz","txt"],help="Upload your GENOME file here for cleaning.")
min_len = st.number_input("Minimum sequence length", value=100, step=10,help="Sequences shorter than this will be removed")
max_n = st.slider("Max fraction of N allowed", min_value=0.0, max_value=0.5, value=0.05, step=0.01,help="Sequences with a higher fraction of 'N' bases will be removed")
dedup = st.checkbox("Enable deduplication", value=True,help="Remove duplicate sequences from the FASTA")

if uploaded is not None:
    tmp_in = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    tmp_in.write(uploaded.getbuffer())
    tmp_in.flush()

    run = st.button("Start cleaning 🧹")

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
        st.success("Cleaning finished✅")


        #saves file path and counters
        st.session_state["counters"]=counters
        st.session_state["out_fasta"]=tmp_out.name
        st.session_state["summary_csv"]=tmp_csv.name

if "counters" in st.session_state:
            counters = st.session_state["counters"]

            col1, col2, col3, col4, col5, col6 = st.columns(6)
            col1.metric("Total", counters["total"])
            col2.metric("Kept", counters["kept"])
            col3.metric("Too Short", counters["too_short"])
            col4.metric("High Ns", counters["high_Ns"])
            col5.metric("Invalid codons", counters["non_allowed"])
            col6.metric("Duplicates", counters["duplicate"])

            # Download buttons 
            with open(st.session_state["out_fasta"], "rb") as f:
                st.download_button("Download cleaned file", data=f, file_name="cleaned.fasta")

            with open(st.session_state["summary_csv"], "rb") as f:
                st.download_button("Download summary CSV", data=f, file_name="summary.csv")





        