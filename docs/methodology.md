# Methodology (MVP)
- InfiltrationScore = 0.35*AWC_norm + 0.35*(1 - Slope_norm) + 0.30*LC_norm
- AWC: Available Water Capacity (m³/m³), provided via Google Drive link
- DEM: derive slope with Horn gradient; normalize 0-1 by percentile clamping
- Land cover (optional): if provided, convert to permeability index (0-1)
