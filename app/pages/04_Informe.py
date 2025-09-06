import sys, os
_APP_DIR = os.path.dirname(os.path.dirname(__file__))
if _APP_DIR not in sys.path: sys.path.insert(0, _APP_DIR)
from state import get_state

import streamlit as st
from pathlib import Path

S = get_state()
st.set_page_config(page_title="Informe", layout="wide")
st.title("Informe (Executive Summary)")
st.write("Downloadables & cloud export")

st.markdown("---")
st.markdown("### Export to Google Drive")
drive_folder_link = st.text_input("Google Drive Folder Link (share this folder with your service account email)", key="drive_folder_link_input")

if st.button("☁️ Export outputs to Drive", key="export_to_drive_btn"):
    try:
        from etl.gdrive_export import folder_id_from_link, upload_file
        folder_id = folder_id_from_link(drive_folder_link or "")
        if not folder_id:
            st.error("Invalid folder link. It should contain /folders/<FOLDER_ID> or ?id=<FOLDER_ID>.")
        else:
            uploaded = []
            infil = S.get("infiltration_path")
            aoi = S.get("aoi_geojson")
            if infil and Path(infil).exists():
                _, link = upload_file(Path(infil), folder_id, mimetype="image/tiff")
                uploaded.append(("InfiltrationScore.tif", link))
            if aoi and Path(aoi).exists():
                _, link = upload_file(Path(aoi), folder_id, mimetype="application/geo+json")
                uploaded.append(("AOI.geojson", link))
            if uploaded:
                st.success("Uploaded to Google Drive:")
                for name, link in uploaded:
                    st.write(f"• **{name}** → {link}")
            else:
                st.warning("No output files found to upload yet. Run the analysis first.")
    except Exception as e:
        st.error(f"Drive export failed: {e}\n\nIn Streamlit Cloud, set `gdrive_service_account` in Secrets (service account JSON).")

