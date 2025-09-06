from app.state import get_state
import streamlit as st
import pandas as pd

st.set_page_config(page_title="Dashboard", layout="wide")
S = get_state()

st.title("Dashboard")

st.subheader("Comparative Analysis: Site vs AOI vs Region")
left, right = st.columns([3,2])

with left:
    st.write("Infiltration Score (0–100)")
    df = pd.DataFrame({
        "Area": ["Site","AOI","Region"],
        "Infiltration": [72,68,54]
    })
    st.bar_chart(df, x="Area", y="Infiltration", height=300)

with right:
    st.write("Land Cover Composition (%)")
    st.dataframe(pd.DataFrame({
        "Cover": ["Forest","Agriculture","Urban"],
        "Site": [45,30,25],
        "AOI": [38,42,20],
        "Region": [40,35,25]
    }))

st.markdown("---")
st.subheader("Key Metrics")
metrics = pd.DataFrame({
    "Metric":["AWC avg (m³/m³)","Slope avg (°)","Floodplain %","Natura 2000 %"],
    "Site":[0.18,2.3,12,8],
    "AOI":[0.16,2.8,10,7],
    "Region":[0.15,3.1,9,6]
})
st.dataframe(metrics, use_container_width=True)
st.download_button("📥 Export Data (CSV)", metrics.to_csv(index=False).encode("utf-8"), "dashboard_metrics.csv")
