"""
time_series_analysis.py
------------------------
Small Baseline Subset (SBAS) time-series inversion for multi-temporal InSAR.

Stacks multiple interferograms to produce a displacement time-series at each
coherent pixel, then fits a linear + seasonal model to extract:
  - Linear velocity (mm/year)
  - Annual amplitude and phase
  - Temporal coherence

Reference: Berardino et al. (2002), Casu et al. (2006)
"""

import logging
from pathlib import Path

import numpy as np
import rasterio
from scipy import linalg
from scipy.optimize import curve_fit

log = logging.getLogger(__name__)

MM_PER_YEAR = 365.25


def build_design_matrix(
    acquisition_dates: list[str],
    reference_date: str | None = None,
) -> np.ndarray:
    """
    Build the SBAS design matrix for time-series inversion.

    For n acquisition dates and m interferograms, the design matrix G (m×n)
    encodes which date pairs contribute to each interferogram.

    Parameters
    ----------
    acquisition_dates : list[str]  ISO 8601 dates of all acquisitions
    reference_date : str  Reference epoch (default: first date)

    Returns
    -------
    G : np.ndarray, shape (n-1, n-1)  SBAS design matrix (incremental formulation)
    t : np.ndarray  Decimal year for each acquisition
    """
    import pandas as pd

    dates = pd.to_datetime(acquisition_dates)
    if reference_date:
        ref = pd.Timestamp(reference_date)
    else:
        ref = dates[0]

    # Convert to decimal years relative to reference
    t = (dates - ref).days / MM_PER_YEAR

    return t.values


def fit_linear_model(
    t: np.ndarray,
    displacement_stack: np.ndarray,
    n_jobs: int = -1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Fit linear + seasonal displacement model at each pixel.

    Model: d(t) = v*t + A*sin(2π*t) + B*cos(2π*t) + offset

    Parameters
    ----------
    t : np.ndarray, shape (N,)  Decimal years
    displacement_stack : np.ndarray, shape (N, rows, cols)  Displacement time-series

    Returns
    -------
    velocity : np.ndarray, shape (rows, cols)  Linear velocity (mm/yr)
    amplitude : np.ndarray  Annual amplitude (mm)
    phase_deg : np.ndarray  Annual phase (degrees)
    temporal_coherence : np.ndarray  RMS residual (proxy for temporal coherence)
    """
    N, rows, cols = displacement_stack.shape

    # Design matrix: [t, sin(2πt), cos(2πt), 1]
    A = np.column_stack([
        t,
        np.sin(2 * np.pi * t),
        np.cos(2 * np.pi * t),
        np.ones(N),
    ])

    # Reshape for vectorized least squares
    D = displacement_stack.reshape(N, rows * cols)

    # Replace NaN columns with zeros for lstsq (will be masked in output)
    valid = ~np.any(np.isnan(D), axis=0)
    D_clean = np.where(np.isnan(D), 0, D)

    # Least squares solution: x = (AᵀA)⁻¹AᵀD
    x, residuals, rank, sv = linalg.lstsq(A, D_clean, check_finite=False)

    velocity_flat = x[0]  # mm/yr
    sin_amp = x[1]
    cos_amp = x[2]

    # Annual amplitude and phase
    amplitude_flat = np.sqrt(sin_amp**2 + cos_amp**2)
    phase_flat = np.degrees(np.arctan2(sin_amp, cos_amp))

    # Temporal coherence (inverse of normalized RMS residual)
    fitted = A @ x
    rms_flat = np.sqrt(np.nanmean((D - fitted)**2, axis=0))
    # Normalize so 0 = incoherent, 1 = perfectly coherent
    max_rms = np.nanpercentile(rms_flat[valid], 95)
    temporal_coherence_flat = np.clip(1.0 - rms_flat / (max_rms + 1e-10), 0, 1)

    # Mask invalid pixels
    velocity_flat[~valid] = np.nan
    amplitude_flat[~valid] = np.nan
    phase_flat[~valid] = np.nan
    temporal_coherence_flat[~valid] = np.nan

    velocity = velocity_flat.reshape(rows, cols)
    amplitude = amplitude_flat.reshape(rows, cols)
    phase_deg = phase_flat.reshape(rows, cols)
    temporal_coherence = temporal_coherence_flat.reshape(rows, cols)

    log.info(
        f"Velocity field: min={np.nanmin(velocity):.1f}, "
        f"max={np.nanmax(velocity):.1f} mm/yr, "
        f"mean={np.nanmean(velocity):.2f} mm/yr"
    )

    return velocity, amplitude, phase_deg, temporal_coherence


def save_velocity_raster(
    velocity: np.ndarray,
    reference_raster: Path,
    output_path: Path,
) -> None:
    """Save velocity field as GeoTIFF using projection from reference raster."""
    with rasterio.open(reference_raster) as src:
        profile = src.profile.copy()

    profile.update(dtype=rasterio.float32, nodata=np.nan, count=1)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(velocity.astype(np.float32), 1)

    log.info(f"Velocity raster → {output_path}")


def run_sbas(
    displacement_tifs: list[Path],
    acquisition_dates: list[str],
    output_dir: Path = Path("results"),
) -> dict:
    """
    Full SBAS time-series pipeline.

    Parameters
    ----------
    displacement_tifs : list[Path]  Processed displacement TIFs (chronological)
    acquisition_dates : list[str]   Corresponding acquisition dates
    output_dir : Path

    Returns
    -------
    dict with output file paths
    """
    if len(displacement_tifs) < 3:
        raise ValueError("SBAS requires at least 3 interferograms")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Running SBAS on {len(displacement_tifs)} interferograms...")

    # Load all displacement rasters into 3D stack
    stack = []
    ref_path = displacement_tifs[0]
    for tif in displacement_tifs:
        with rasterio.open(tif) as src:
            stack.append(src.read(1).astype(np.float32))

    displacement_stack = np.array(stack)  # shape: (N, rows, cols)
    t = build_design_matrix(acquisition_dates)

    # Fit model
    velocity, amplitude, phase_deg, coherence = fit_linear_model(
        t, displacement_stack
    )

    # Save outputs
    outputs = {}
    for name, data in [
        ("velocity_mm_yr", velocity),
        ("annual_amplitude_mm", amplitude),
        ("annual_phase_deg", phase_deg),
        ("temporal_coherence", coherence),
    ]:
        out_path = output_dir / f"{name}.tif"
        save_velocity_raster(data, ref_path, out_path)
        outputs[name] = out_path

    log.info("SBAS complete. Outputs:")
    for name, path in outputs.items():
        log.info(f"  {name}: {path}")

    return outputs
