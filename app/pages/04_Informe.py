from app.state import get_state
import streamlit as st

st.set_page_config(page_title="Informe", layout="wide")
S = get_state()

st.title("Informe (Executive Summary)")

c1, c2, c3, c4 = st.columns(4)
with c1:
    st.metric("Water Score", "72/100")
with c2:
    st.write("**Key Risks**\n\n• Flooding 12%\n\n• Low AWC patches")
with c3:
    st.write("**Opportunities**\n\n• 125 ha Priority A\n\n• Green infra")
with c4:
    st.write("**Site Rating**\n\n⭐️⭐️⭐️⭐️")

st.markdown("### Key Findings")
st.write("- Above-average infiltration potential\n- 12% of AOI in flood-prone areas\n- Recommend nature-based infiltration in Priority A zones")

st.markdown("### Downloads")
st.download_button("📄 PDF Report", b"", file_name="informe.pdf", disabled=True, key="download_pdf")
st.download_button("📝 Markdown", b"", file_name="informe.md", disabled=True, key="download_md")
st.download_button("🗺️ GeoPackage (All Data)", b"", file_name="outputs.gpkg", disabled=True, key="download_gpkg")

with st.expander("Methodology & Data Sources"):
    st.write("• AWC normalized (0–1)\n• Slope derived from DEM (5m)\n• Infiltration = 0.35*AWC + 0.35*(1 - slope_norm) + 0.30*LC\n• Data: PNOA DEM, user AWC (Drive), CORINE/CLC")
