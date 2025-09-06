# sustai-geo-app (MVP)

Streamlit MVP to assess environmental factors around a site. Focus on **Water** module using **AWC (Available Water Capacity)** + **DEM**.

## Quick start
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app/Home.py
```

## How it works
- Paste public **Google Drive** links for **AWC (GeoTIFF)** and **DEM** in the sidebar.
- Set coordinates and radius, click **Run Analysis**.
- View **Mapa**, **Zonificación**, **Dashboard**, **Informe**.

## Data vault
- Large files are **not committed**. They download into `data/raw/`.
- Outputs go to `data/processed/` and `data/reports/`.

## Roadmap
- MVP: Water module w/ AWC + DEM, base maps, basic KPIs.
- v1.1+: OGC API MAPA (flood, Natura), connectivity, benchmarks.
