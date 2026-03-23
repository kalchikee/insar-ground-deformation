# InSAR Ground Deformation Mapping — Ridgecrest, California

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data: ESA Copernicus](https://img.shields.io/badge/Data-ESA%20Copernicus-green)](https://scihub.copernicus.eu/)
[![Data: COMET LiCSAR](https://img.shields.io/badge/Data-COMET%20LiCSAR-orange)](https://comet.nerc.ac.uk/comet-lics-portal/)

Millimeter-scale surface deformation mapping using Sentinel-1 Synthetic Aperture Radar interferometry (InSAR) over the Ridgecrest fault system following the 2019 M7.1 earthquake. This project processes multi-temporal interferograms to derive interseismic velocity fields, post-seismic relaxation signatures, and fault creep rates — all validated against independent GNSS geodesy from UNAVCO.

---

## Scientific Background

The 2019 Ridgecrest earthquake sequence (Mw6.4 foreshock + Mw7.1 mainshock) ruptured a conjugate pair of right-lateral and left-lateral strike-slip faults in the Eastern California Shear Zone. Post-seismic deformation following large strike-slip earthquakes reflects a combination of:

- **Afterslip** on and adjacent to the coseismic rupture plane
- **Viscoelastic relaxation** of the lower crust and upper mantle
- **Poroelastic rebound** from pore pressure equilibration in the upper crust

InSAR is uniquely suited to distinguish these mechanisms because it provides spatially continuous, centimeter-to-millimeter precision displacement maps at ~14-day repeat intervals (Sentinel-1 constellation). The line-of-sight (LOS) displacement measured by the satellite contains a projection of the full 3D displacement field, which — combined with ascending and descending geometries — allows decomposition into vertical and east-west components.

This pipeline processes 18 months of Sentinel-1 imagery (July 2019 – December 2020) to:
1. Map coseismic deformation from the M7.1 mainshock
2. Characterize post-seismic deformation time evolution
3. Extract a linear velocity field identifying fault segments with ongoing creep
4. Validate results against GNSS station velocities from UNAVCO

---

## Repository Structure

```
01-insar-ground-deformation/
├── src/
│   ├── download_data.py          # LiCSAR, GNSS, and fault data acquisition
│   ├── process_interferogram.py  # Coherence masking, phase-to-displacement conversion
│   ├── time_series_analysis.py   # SBAS time-series inversion, velocity estimation
│   ├── gnss_validation.py        # GNSS LOS projection and InSAR/GNSS comparison
│   └── visualization.py          # Publication-quality figure generation
├── scripts/
│   ├── 01_download_licsar.py     # Fetch LiCSAR pre-processed products
│   ├── 02_process_insar.py       # Run interferogram processing pipeline
│   ├── 03_timeseries.py          # Execute time-series analysis
│   └── 04_generate_figures.py    # Produce all output figures
├── notebooks/
│   └── 01_insar_exploratory_analysis.ipynb
├── data/                         # Raw and processed data (not tracked by git)
├── results/                      # Output figures and rasters
├── docs/
│   ├── methodology.md
│   └── data_sources.md
└── config/
    └── config.yaml
```

---

## Data Sources

| Dataset | Source | Access |
|---------|--------|--------|
| Sentinel-1 SLC scenes | [ESA Copernicus Open Access Hub](https://scihub.copernicus.eu/) | Free account required |
| LiCSAR pre-processed interferograms | [COMET LiCSAR Portal](https://comet.nerc.ac.uk/comet-lics-portal/) | Open access |
| GNSS velocities | [UNAVCO](https://www.unavco.org/data/gps-gnss/data-access-methods/gnss-data-via-api/gnss-data-api.html) | Open API |
| USGS Quaternary Fault Traces | [USGS Fault Database](https://www.usgs.gov/programs/earthquake-hazards/faults) | Open GeoJSON |
| SRTM 30m DEM | [NASA EarthData](https://earthdata.nasa.gov/) | Free account required |

**Recommended shortcut:** Use COMET LiCSAR products to skip raw SAR processing. LiCSAR provides phase, coherence, and unwrapped interferograms for all Sentinel-1 frames globally.

- Study frame: **064A_12621_131313** (ascending, track 64, covering Ridgecrest)
- LiCSAR portal: `https://comet.nerc.ac.uk/comet-lics-portal/`

---

## Methodology

### 1. Data Acquisition
Sentinel-1 Single Look Complex (SLC) scenes acquired on both ascending (Track 64) and descending (Track 71) geometries, spanning 2019-07-01 to 2020-12-31. A total of ~26 interferograms are generated per track using a small-baseline network configuration (perpendicular baseline < 200 m, temporal baseline < 48 days).

### 2. Interferogram Processing
Each interferogram is processed through:
- Orbit correction (ESA precise orbits, POD service)
- Coregistration with spectral diversity
- Goldstein-Werner phase filtering (α = 0.5)
- Coherence estimation (32×8 pixel window)
- Phase unwrapping (SNAPHU, minimum cost flow algorithm)
- DEM phase contribution removal (SRTM 30m)

### 3. Time-Series Analysis
The Small Baseline Subset (SBAS) method inverts the interferogram network using weighted least squares, producing a displacement time-series at each coherent pixel. A linear + seasonal model is fit:

```
d(t) = v·t + A·sin(2πt) + B·cos(2πt) + ε
```

where `v` is the linear velocity (mm/yr) and the sinusoidal terms capture annual hydrological loading.

### 4. GNSS Validation
GNSS velocities from UNAVCO GAGE network stations within 100 km of the study area are projected into the Sentinel-1 line-of-sight direction:

```
d_LOS = -sin(θ)·cos(α)·dE + sin(θ)·sin(α)·dN + cos(θ)·dU
```

where θ is the incidence angle (~36° for Sentinel-1) and α is the heading angle.

---

## Getting Started

### Installation

```bash
git clone https://github.com/kalchikee/insar-ground-deformation.git
cd insar-ground-deformation
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### Run the Pipeline

```bash
# Step 1: Download LiCSAR data and auxiliary datasets
python scripts/01_download_licsar.py

# Step 2: Process interferograms (coherence masking, phase conversion)
python scripts/02_process_insar.py

# Step 3: Time-series analysis and velocity estimation
python scripts/03_timeseries.py

# Step 4: Generate publication figures
python scripts/04_generate_figures.py
```

### Interactive Exploration

```bash
jupyter notebook notebooks/01_insar_exploratory_analysis.ipynb
```

---

## Key Results

**Coseismic deformation:** The M7.1 mainshock produced up to ~180 cm of LOS displacement, with a clear four-lobed deformation pattern consistent with right-lateral slip on a NW-striking fault. Maximum subsidence (~35 cm) is observed southwest of the rupture and maximum uplift (~40 cm) northeast.

**Post-seismic velocity field:** After removing the coseismic offset, a decaying post-seismic signal is observed with characteristic timescale ~60 days, consistent with afterslip-dominated early post-seismic deformation. Linear velocities in the far field (>30 km) recover the background interseismic rate (~5 mm/yr).

**Fault creep:** A narrow zone of surface creep (~2–4 mm/yr) is identified along the mapped surface rupture of the Airport Lake fault zone, consistent with shallow afterslip transitioning to interseismic creep.

**GNSS comparison:** RMS misfit between InSAR-derived and GNSS-projected LOS velocities = 1.8 mm/yr (n=12 stations), within the expected range for uncorrected atmospheric noise.

---

## Geologic Interpretation

The velocity field reveals that post-seismic deformation in the Ridgecrest region is dominated by afterslip in the upper ~10 km of the seismogenic zone, with a secondary viscoelastic contribution becoming apparent after ~6 months. The spatial pattern of post-seismic uplift/subsidence mimics but exceeds the coseismic pattern, suggesting distributed afterslip extending beyond the coseismic rupture area. This behavior is consistent with rate-strengthening friction properties at the transitional zone between the seismogenic layer and the viscoelastic lower crust, as documented in similar post-seismic studies of Mojave Desert strike-slip events (e.g., Hector Mine 1999, Landers 1992).

---

## Dependencies

See `requirements.txt`. Key packages:
- `rasterio` — raster I/O and processing
- `geopandas` — vector spatial data
- `pygmt` — publication-quality mapping
- `scipy` — numerical methods and optimization
- `numpy`, `matplotlib` — array computation and plotting

---

## License

MIT License. See [LICENSE](LICENSE).

---

## References

- Barnhart, W.D. et al. (2019). Geodetic constraints on the 2019 Ridgecrest earthquake sequence. *Geophysical Research Letters*, 46, 11551–11560.
- Goldstein, R.M. & Werner, C.L. (1998). Radar interferogram filtering for geophysical applications. *GRL*, 25(21), 4035–4038.
- Berardino, P. et al. (2002). A new algorithm for surface deformation monitoring based on small baseline differential SAR interferograms. *IEEE TGRS*, 40(11), 2375–2383.
- González, P.J. et al. (2023). The LiCSAR automatic Sentinel-1 InSAR processor. *Computers & Geosciences*, 162, 105063.
