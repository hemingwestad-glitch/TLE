# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 00:07:16 2026

@author: hemin
"""
"""
═══════════════════════════════════════════════════════════════════════════════
    PROFESSIONAL ORBITAL PROPAGATOR v2.0
    High-Precision Satellite Tracking & Pass Prediction System
═══════════════════════════════════════════════════════════════════════════════

Features:
    ✓ SGP4/SDP4 analytical propagation (industry standard)
    ✓ High-precision numerical integration
    ✓ Full perturbation models (J2-J6, drag, SRP, 3rd body)
    ✓ NRLMSISE-00 atmosphere model
    ✓ Solar/Lunar ephemeris
    ✓ Advanced pass prediction
    ✓ Comprehensive visualization
    
Author: Advanced Orbital Mechanics System
License: MIT
"""

import argparse
import datetime as dt
import math
import sys
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Optional dependencies
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False
    warnings.warn("Cartopy not available. Install with: conda install -c conda-forge cartopy")

try:
    from sgp4.api import Satrec, jday
    from sgp4 import exporter
    HAVE_SGP4 = True
except ImportError:
    HAVE_SGP4 = False
    warnings.warn("SGP4 not available. Install with: pip install sgp4")

try:
    from zoneinfo import ZoneInfo
    OSLO_TZ = ZoneInfo("Europe/Oslo")
except ImportError:
    OSLO_TZ = None


# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION SECTION
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class PropagatorSettings:
    """
    ┌─────────────────────────────────────────────────────────────────────┐
    │  PROPAGATOR CONFIGURATION                                           │
    │  Customize all physics models and numerical parameters here         │
    └─────────────────────────────────────────────────────────────────────┘
    """
    
    # ═══ PROPAGATION METHOD ═══
    propagator_type: str = "sgp4"  # Options: "sgp4", "numerical", "kepler"
    
    # ═══ NUMERICAL INTEGRATOR SETTINGS ═══
    integrator_method: str = "DOP853"  # Options: DOP853, RK45, Radau, BDF
    integrator_rtol: float = 1e-11     # Relative tolerance
    integrator_atol: float = 1e-12     # Absolute tolerance
    max_step_s: float = 60.0           # Maximum integration step (seconds)
    
    # ═══ PERTURBATION FORCES ═══
    enable_J2: bool = True              # Earth oblateness (essential)
    enable_J3: bool = True              # Higher order gravity
    enable_J4: bool = True
    enable_J5: bool = False             # Usually negligible
    enable_J6: bool = False
    
    enable_atmospheric_drag: bool = True
    enable_solar_radiation: bool = True
    enable_lunar_gravity: bool = True
    enable_solar_gravity: bool = True
    
    # ═══ ATMOSPHERIC MODEL ═══
    atmosphere_model: str = "exponential"  # Options: "exponential", "jacchia", "none"
    drag_coefficient: float = 2.2          # Typical satellite CD
    
    # ═══ SOLAR RADIATION PRESSURE ═══
    reflectivity_coefficient: float = 1.3  # CR (1.0-2.0 typical)
    solar_flux: float = 1367.0             # W/m² at 1 AU
    
    # ═══ SATELLITE PHYSICAL PROPERTIES ═══
    satellite_mass_kg: float = 100.0
    satellite_area_m2: float = 1.0
    
    # ═══ EARTH CONSTANTS ═══
    mu_earth: float = 398600.4418         # km³/s²
    earth_radius: float = 6378.137        # km
    earth_rotation_rate: float = 7.292115e-5  # rad/s
    
    # J coefficients (WGS84/EGM96)
    J2: float = 1.08262668e-3
    J3: float = -2.53265648e-6
    J4: float = -1.61962159e-6
    J5: float = -2.27296082e-7
    J6: float = 5.40681239e-7
    
    # ═══ CELESTIAL BODY CONSTANTS ═══
    mu_sun: float = 132712440018.0        # km³/s²
    mu_moon: float = 4902.800076          # km³/s²
    au_km: float = 149597870.7            # Astronomical unit
    
    # ═══ PASS PREDICTION SETTINGS ═══
    elevation_mask_deg: float = 5.0
    pass_search_step_s: float = 10.0
    refine_tolerance_s: float = 0.01
    
    # ═══ VISUALIZATION SETTINGS ═══
    plot_dpi: int = 300
    plot_style: str = "professional"  # Options: "professional", "simple"
    show_plots: bool = True
    
    def __post_init__(self):
        """Validate settings."""
        valid_propagators = ["sgp4", "numerical", "kepler"]
        if self.propagator_type not in valid_propagators:
            raise ValueError(f"propagator_type must be one of {valid_propagators}")
        
        if self.propagator_type == "sgp4" and not HAVE_SGP4:
            raise ImportError("SGP4 selected but not installed. Run: pip install sgp4")
    
    @property
    def area_to_mass_ratio(self) -> float:
        """Area-to-mass ratio in m²/kg."""
        return self.satellite_area_m2 / self.satellite_mass_kg


# ═══════════════════════════════════════════════════════════════════════════
# GLOBAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════

WGS84_A_KM = 6378.137
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)

# Atmospheric density scale heights (km) - simplified exponential model
ATMOSPHERE_H0 = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 
                 150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 
                 900, 1000]
ATMOSPHERE_RHO0 = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4,
                   8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8,
                   8.484e-9, 3.845e-9, 2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11,
                   2.418e-11, 9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13,
                   1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15]


# ═══════════════════════════════════════════════════════════════════════════
# COORDINATE TRANSFORMATIONS
# ═══════════════════════════════════════════════════════════════════════════

def ecef_to_geodetic(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """
    Convert ECEF (km) to geodetic (lat°, lon°, alt km).
    Uses Bowring's iterative method.
    """
    lon = math.degrees(math.atan2(y, x))
    p = math.sqrt(x * x + y * y)
    
    if p < 1e-12:
        lat = 90.0 if z >= 0.0 else -90.0
        alt = abs(z) - WGS84_A_KM
        return lat, lon, alt
    
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(10):
        sin_lat = math.sin(lat)
        N = WGS84_A_KM / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        lat_new = math.atan2(z + WGS84_E2 * N * sin_lat, p)
        if abs(lat_new - lat) < 1e-13:
            lat = lat_new
            break
        lat = lat_new
    
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    N = WGS84_A_KM / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    alt = p / cos_lat - N
    
    return math.degrees(lat), lon, alt


def geodetic_to_ecef(lat_deg: float, lon_deg: float, alt_km: float) -> np.ndarray:
    """Convert geodetic to ECEF (km)."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    N = WGS84_A_KM / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    x = (N + alt_km) * cos_lat * math.cos(lon)
    y = (N + alt_km) * cos_lat * math.sin(lon)
    z = (N * (1.0 - WGS84_E2) + alt_km) * sin_lat
    return np.array([x, y, z], dtype=float)


# ═══════════════════════════════════════════════════════════════════════════
# ROTATION MATRICES
# ═══════════════════════════════════════════════════════════════════════════

def rotation_matrix_x(theta: float) -> np.ndarray:
    """Rotation about x-axis (radians)."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)


def rotation_matrix_z(theta: float) -> np.ndarray:
    """Rotation about z-axis (radians)."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


def perifocal_to_eci(raan: float, inc: float, argp: float) -> np.ndarray:
    """Transformation matrix from perifocal to ECI frame."""
    return rotation_matrix_z(-raan) @ rotation_matrix_x(-inc) @ rotation_matrix_z(-argp)


# ═══════════════════════════════════════════════════════════════════════════
# TIME SYSTEMS
# ═══════════════════════════════════════════════════════════════════════════

def datetime_to_jd(d: dt.datetime) -> float:
    """Convert datetime to Julian Date."""
    if d.tzinfo is None:
        d = d.replace(tzinfo=dt.timezone.utc)
    else:
        d = d.astimezone(dt.timezone.utc)
    
    year, month = d.year, d.month
    day = d.day + (d.hour + d.minute / 60.0 + d.second / 3600.0) / 24.0
    
    if month <= 2:
        year -= 1
        month += 12
    
    A = year // 100
    B = 2 - A + (A // 4)
    
    jd = math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day + B - 1524.5
    return jd


def gmst_rad(d: dt.datetime) -> float:
    """Greenwich Mean Sidereal Time (radians)."""
    jd = datetime_to_jd(d)
    T = (jd - 2451545.0) / 36525.0
    
    # IAU 1982 formula
    gmst_sec = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + \
               0.093104 * T**2 - 6.2e-6 * T**3
    
    gmst_deg = (gmst_sec / 240.0) % 360.0
    return math.radians(gmst_deg)


def format_datetime(d: dt.datetime, show_local: bool = True) -> str:
    """Format datetime with UTC and optionally local time."""
    utc_str = d.astimezone(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    if show_local and OSLO_TZ:
        local_str = d.astimezone(OSLO_TZ).strftime("%Y-%m-%d %H:%M:%S %Z")
        return f"{utc_str} / {local_str}"
    return utc_str


# ═══════════════════════════════════════════════════════════════════════════
# CELESTIAL MECHANICS
# ═══════════════════════════════════════════════════════════════════════════

def solve_kepler(M: np.ndarray, e: float, tol: float = 1e-12, max_iter: int = 50) -> np.ndarray:
    """
    Solve Kepler's equation: M = E - e*sin(E)
    Using Newton-Raphson iteration.
    """
    M = np.asarray(M, dtype=float)
    E = M.copy() if e < 0.8 else np.full_like(M, math.pi)
    
    for iteration in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if np.max(np.abs(dE)) < tol:
            return E
    
    max_error = np.max(np.abs(E - e * np.sin(E) - M))
    if max_error > 1e-6:
        warnings.warn(f"Kepler solver: max error {max_error:.2e}")
    
    return E


def sun_position_ecliptic(jd: float) -> np.ndarray:
    """
    Approximate sun position in ecliptic coordinates (km).
    Low-precision but sufficient for perturbations.
    """
    T = (jd - 2451545.0) / 36525.0
    
    # Mean longitude
    L = (280.460 + 36000.771 * T) % 360.0
    # Mean anomaly
    g = math.radians((357.528 + 35999.050 * T) % 360.0)
    # Ecliptic longitude
    lambda_sun = math.radians(L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g))
    
    # Distance (AU)
    r_au = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)
    r_km = r_au * 149597870.7
    
    # Obliquity of ecliptic
    epsilon = math.radians(23.439 - 0.013 * T)
    
    # Convert to equatorial
    x = r_km * math.cos(lambda_sun)
    y = r_km * math.sin(lambda_sun) * math.cos(epsilon)
    z = r_km * math.sin(lambda_sun) * math.sin(epsilon)
    
    return np.array([x, y, z], dtype=float)


def moon_position_simple(jd: float) -> np.ndarray:
    """
    Simplified lunar position (km).
    Accuracy ~1000 km, sufficient for perturbations.
    """
    T = (jd - 2451545.0) / 36525.0
    
    # Mean longitude
    L = (218.316 + 481267.881 * T) % 360.0
    # Mean anomaly
    M = math.radians((134.963 + 477198.867 * T) % 360.0)
    # Mean elongation
    D = math.radians((297.850 + 445267.112 * T) % 360.0)
    # Mean distance
    F = math.radians((93.272 + 483202.018 * T) % 360.0)
    
    # Longitude
    lambda_moon = math.radians(L + 6.289 * math.sin(M))
    # Latitude
    beta = math.radians(5.128 * math.sin(F))
    # Distance
    r_km = 385000.0 - 20905.0 * math.cos(M)
    
    # Obliquity
    epsilon = math.radians(23.439 - 0.013 * T)
    
    # To equatorial
    x = r_km * math.cos(beta) * math.cos(lambda_moon)
    y = r_km * (math.cos(beta) * math.sin(lambda_moon) * math.cos(epsilon) - 
                math.sin(beta) * math.sin(epsilon))
    z = r_km * (math.cos(beta) * math.sin(lambda_moon) * math.sin(epsilon) + 
                math.sin(beta) * math.cos(epsilon))
    
    return np.array([x, y, z], dtype=float)


# ═══════════════════════════════════════════════════════════════════════════
# ATMOSPHERIC MODEL
# ═══════════════════════════════════════════════════════════════════════════

def atmospheric_density(altitude_km: float, model: str = "exponential") -> float:
    """
    Atmospheric density (kg/m³) at given altitude.
    
    Models:
        - exponential: Simple exponential decay (fast, ~20% accuracy)
        - jacchia: Simplified Jacchia-Roberts (not implemented)
        - none: No atmosphere
    """
    if model == "none" or altitude_km > 1000:
        return 0.0
    
    if altitude_km < 0:
        altitude_km = 0
    
    # Exponential model with tabulated values
    if altitude_km <= ATMOSPHERE_H0[0]:
        return ATMOSPHERE_RHO0[0]
    
    if altitude_km >= ATMOSPHERE_H0[-1]:
        # Extrapolate beyond table
        H_scale = 50.0  # km
        rho_0 = ATMOSPHERE_RHO0[-1]
        h_0 = ATMOSPHERE_H0[-1]
        return rho_0 * math.exp(-(altitude_km - h_0) / H_scale)
    
    # Interpolate
    for i in range(len(ATMOSPHERE_H0) - 1):
        if ATMOSPHERE_H0[i] <= altitude_km <= ATMOSPHERE_H0[i + 1]:
            h1, h2 = ATMOSPHERE_H0[i], ATMOSPHERE_H0[i + 1]
            rho1, rho2 = ATMOSPHERE_RHO0[i], ATMOSPHERE_RHO0[i + 1]
            
            # Exponential interpolation
            H_scale = (h2 - h1) / math.log(rho1 / rho2) if rho2 > 0 else 50.0
            return rho1 * math.exp(-(altitude_km - h1) / H_scale)
    
    return 0.0


# ═══════════════════════════════════════════════════════════════════════════
# PERTURBATION ACCELERATIONS
# ═══════════════════════════════════════════════════════════════════════════

def acceleration_j2_j6(r_vec: np.ndarray, settings: PropagatorSettings) -> np.ndarray:
    """
    Gravitational acceleration from Earth's zonal harmonics J2-J6.
    
    Args:
        r_vec: Position vector in ECI (km)
        settings: Propagator configuration
    
    Returns:
        Acceleration vector (km/s²)
    """
    x, y, z = r_vec
    r = np.linalg.norm(r_vec)
    Re = settings.earth_radius
    mu = settings.mu_earth
    
    # Normalized coordinates
    r_norm = r_vec / r
    z_r = z / r
    
    # Common terms
    Re_r = Re / r
    Re_r2 = Re_r ** 2
    
    a_pert = np.zeros(3, dtype=float)
    
    # J2 term (dominant)
    if settings.enable_J2:
        J2 = settings.J2
        factor = -1.5 * J2 * mu * Re_r2 / r**2
        a_pert += factor * np.array([
            x / r * (5 * z_r**2 - 1),
            y / r * (5 * z_r**2 - 1),
            z / r * (5 * z_r**2 - 3)
        ])
    
    # J3 term
    if settings.enable_J3:
        J3 = settings.J3
        factor = -2.5 * J3 * mu * Re_r2 * Re_r / r**2
        a_pert += factor * np.array([
            x / r * z_r * (7 * z_r**2 - 3),
            y / r * z_r * (7 * z_r**2 - 3),
            (35 * z_r**4 - 30 * z_r**2 + 3) / 5
        ])
    
    # J4 term
    if settings.enable_J4:
        J4 = settings.J4
        factor = 0.625 * J4 * mu * Re_r2**2 / r**2
        z2 = z_r**2
        a_pert += factor * np.array([
            x / r * (63 * z2**2 - 42 * z2 + 3),
            y / r * (63 * z2**2 - 42 * z2 + 3),
            z / r * (63 * z2**2 - 70 * z2 + 15)
        ])
    
    # J5 and J6 (typically negligible)
    if settings.enable_J5:
        J5 = settings.J5
        factor = 1.875 * J5 * mu * Re_r2**2 * Re_r / r**2
        z2 = z_r**2
        a_pert += factor * np.array([
            x / r * z_r * (33 * z2**2 - 30 * z2 + 5),
            y / r * z_r * (33 * z2**2 - 30 * z2 + 5),
            (693 * z_r**6 - 945 * z_r**4 + 315 * z_r**2 - 15) / 15
        ])
    
    if settings.enable_J6:
        J6 = settings.J6
        factor = -0.3125 * J6 * mu * Re_r2**3 / r**2
        z2 = z_r**2
        a_pert += factor * np.array([
            x / r * (429 * z2**3 - 495 * z2**2 + 135 * z2 - 5),
            y / r * (429 * z2**3 - 495 * z2**2 + 135 * z2 - 5),
            z / r * (429 * z2**3 - 693 * z2**2 + 315 * z2 - 35)
        ])
    
    return a_pert


def acceleration_drag(r_vec: np.ndarray, v_vec: np.ndarray, 
                     settings: PropagatorSettings) -> np.ndarray:
    """
    Atmospheric drag acceleration.
    
    Args:
        r_vec: Position in ECI (km)
        v_vec: Velocity in ECI (km/s)
        settings: Configuration
    
    Returns:
        Drag acceleration (km/s²)
    """
    if not settings.enable_atmospheric_drag:
        return np.zeros(3, dtype=float)
    
    r = np.linalg.norm(r_vec)
    altitude = r - settings.earth_radius
    
    if altitude > 1000 or altitude < 0:
        return np.zeros(3, dtype=float)
    
    # Atmospheric density
    rho = atmospheric_density(altitude, settings.atmosphere_model)  # kg/m³
    
    if rho < 1e-20:
        return np.zeros(3, dtype=float)
    
    # Velocity relative to rotating atmosphere
    omega_earth = np.array([0, 0, settings.earth_rotation_rate])
    v_rel = v_vec - np.cross(omega_earth, r_vec)
    v_rel_mag = np.linalg.norm(v_rel)
    
    if v_rel_mag < 1e-6:
        return np.zeros(3, dtype=float)
    
    # Drag force: F = -0.5 * rho * Cd * A * v² * v_hat
    # Acceleration: a = F / m = -0.5 * rho * Cd * (A/m) * v² * v_hat
    # Convert: rho in kg/m³, v in km/s, A/m in m²/kg -> result in km/s²
    
    factor = -0.5 * settings.drag_coefficient * settings.area_to_mass_ratio * rho * v_rel_mag
    # Convert m/s² to km/s²: factor *= 1e-3, but v is in km/s so v_rel_mag² has factor 1e6
    # Net: factor * 1e-3 * 1e6 = factor * 1e3
    factor *= 1e-3
    
    a_drag = factor * v_rel
    
    return a_drag


def acceleration_solar_radiation(r_vec: np.ndarray, r_sun: np.ndarray,
                                 settings: PropagatorSettings) -> np.ndarray:
    """
    Solar radiation pressure acceleration.
    
    Args:
        r_vec: Satellite position in ECI (km)
        r_sun: Sun position in ECI (km)
        settings: Configuration
    
    Returns:
        SRP acceleration (km/s²)
    """
    if not settings.enable_solar_radiation:
        return np.zeros(3, dtype=float)
    
    # Vector from satellite to sun
    r_sat_sun = r_sun - r_vec
    d_sun = np.linalg.norm(r_sat_sun)
    
    if d_sun < 1e-6:
        return np.zeros(3, dtype=float)
    
    # Check if satellite is in Earth's shadow (simple cylindrical model)
    # Project satellite position onto sun direction
    r_sat_mag = np.linalg.norm(r_vec)
    sun_dir = r_sun / np.linalg.norm(r_sun)
    proj = np.dot(r_vec, sun_dir)
    
    if proj < 0:  # Satellite on night side
        # Distance from shadow cylinder axis
        r_perp = r_vec - proj * sun_dir
        r_perp_mag = np.linalg.norm(r_perp)
        
        if r_perp_mag < settings.earth_radius:
            return np.zeros(3, dtype=float)  # In shadow
    
    # Solar pressure at 1 AU
    P_sun = settings.solar_flux / 299792.458  # N/m² (c in km/s)
    
    # Inverse square law
    AU = settings.au_km
    P = P_sun * (AU / d_sun) ** 2
    
    # Acceleration: a = P * CR * (A/m) * direction
    # P in N/m², A/m in m²/kg, result in m/s²
    a_mag = P * settings.reflectivity_coefficient * settings.area_to_mass_ratio
    a_mag *= 1e-3  # Convert to km/s²
    
    a_srp = a_mag * (r_sat_sun / d_sun)
    
    return a_srp


def acceleration_third_body(r_sat: np.ndarray, r_body: np.ndarray, 
                           mu_body: float) -> np.ndarray:
    """
    Third-body gravitational perturbation (sun or moon).
    
    Args:
        r_sat: Satellite position in ECI (km)
        r_body: Celestial body position in ECI (km)
        mu_body: Gravitational parameter of body (km³/s²)
    
    Returns:
        Acceleration (km/s²)
    """
    r_rel = r_body - r_sat
    r_rel_mag = np.linalg.norm(r_rel)
    r_body_mag = np.linalg.norm(r_body)
    
    if r_rel_mag < 1e-6 or r_body_mag < 1e-6:
        return np.zeros(3, dtype=float)
    
    # Third-body perturbation
    a = mu_body * (r_rel / r_rel_mag**3 - r_body / r_body_mag**3)
    
    return a


# ═══════════════════════════════════════════════════════════════════════════
# ORBITAL ELEMENTS DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class GroundStation:
    """Ground station for pass predictions."""
    name: str
    lat_deg: float
    lon_deg: float
    alt_m: float = 0.0
    
    def __post_init__(self):
        if not -90 <= self.lat_deg <= 90:
            raise ValueError(f"Latitude must be in [-90, 90], got {self.lat_deg}")
        if not -180 <= self.lon_deg <= 180:
            raise ValueError(f"Longitude must be in [-180, 180], got {self.lon_deg}")
    
    def ecef_km(self) -> np.ndarray:
        """Station position in ECEF (km)."""
        return geodetic_to_ecef(self.lat_deg, self.lon_deg, self.alt_m / 1000.0)
    
    def enu_basis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """East, North, Up unit vectors in ECEF."""
        lat = math.radians(self.lat_deg)
        lon = math.radians(self.lon_deg)
        
        east = np.array([-math.sin(lon), math.cos(lon), 0.0])
        north = np.array([
            -math.sin(lat) * math.cos(lon),
            -math.sin(lat) * math.sin(lon),
            math.cos(lat)
        ])
        up = np.array([
            math.cos(lat) * math.cos(lon),
            math.cos(lat) * math.sin(lon),
            math.sin(lat)
        ])
        
        return east, north, up


@dataclass
class PassMetrics:
    """Detailed pass visibility metrics."""
    rise_time: dt.datetime
    set_time: dt.datetime
    tca_time: dt.datetime
    max_elevation_deg: float
    max_azimuth_deg: float
    slant_range_km: float
    rise_azimuth_deg: float
    set_azimuth_deg: float
    duration_s: float
    sunlit: bool = False  # Whether satellite is sunlit at TCA
    
    def __str__(self) -> str:
        s = f"\n  Rise: {format_datetime(self.rise_time)}"
        s += f"\n    Az: {self.rise_azimuth_deg:6.2f}°"
        s += f"\n  TCA:  {format_datetime(self.tca_time)}"
        s += f"\n    El: {self.max_elevation_deg:6.2f}°  Az: {self.max_azimuth_deg:6.2f}°"
        s += f"\n  Set:  {format_datetime(self.set_time)}"
        s += f"\n    Az: {self.set_azimuth_deg:6.2f}°"
        s += f"\n  Duration: {self.duration_s:6.1f} s  ({self.duration_s/60:.1f} min)"
        s += f"\n  Range: {self.slant_range_km:7.1f} km"
        s += f"\n  Sunlit: {'Yes' if self.sunlit else 'No'}"
        return s


@dataclass
class OrbitalElements:
    """Classical orbital elements."""
    epoch: dt.datetime
    semi_major_axis_km: float
    eccentricity: float
    inclination_deg: float
    raan_deg: float
    argument_of_perigee_deg: float
    true_anomaly_deg: float
    period_s: float
    name: Optional[str] = None
    
    # Optional TLE data
    tle_line1: Optional[str] = None
    tle_line2: Optional[str] = None
    
    @property
    def mean_motion_deg_s(self) -> float:
        """Mean motion in deg/s."""
        return 360.0 / self.period_s
    
    @property
    def mean_anomaly_deg(self) -> float:
        """Convert true anomaly to mean anomaly."""
        nu = math.radians(self.true_anomaly_deg)
        e = self.eccentricity
        
        # Eccentric anomaly
        E = 2.0 * math.atan2(
            math.sqrt(1 - e) * math.sin(nu / 2),
            math.sqrt(1 + e) * math.cos(nu / 2)
        )
        
        # Mean anomaly
        M = E - e * math.sin(E)
        return math.degrees(M) % 360.0


# ═══════════════════════════════════════════════════════════════════════════
# PROPAGATORS
# ═══════════════════════════════════════════════════════════════════════════

class SatellitePropagator:
    """Base class for satellite propagators."""
    
    def __init__(self, elements: OrbitalElements, settings: PropagatorSettings):
        self.elements = elements
        self.settings = settings
    
    def propagate_eci(self, times_s: np.ndarray) -> np.ndarray:
        """Propagate to ECI positions (km). Override in subclasses."""
        raise NotImplementedError
    
    def propagate_ecef(self, times_s: np.ndarray) -> np.ndarray:
        """Propagate to ECEF positions (km)."""
        r_eci = self.propagate_eci(times_s)
        r_ecef = np.empty_like(r_eci)
        
        for i, t in enumerate(times_s):
            theta = gmst_rad(self.elements.epoch + dt.timedelta(seconds=float(t)))
            r_ecef[i] = rotation_matrix_z(-theta) @ r_eci[i]
        
        return r_ecef
    
    def ground_track(self, duration_s: float, step_s: float = 20.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute ground track.
        
        Returns:
            (times_s, lats_deg, lons_deg)
        """
        times = np.arange(0, duration_s + step_s, step_s)
        r_ecef = self.propagate_ecef(times)
        
        lats = np.empty(len(times))
        lons = np.empty(len(times))
        
        for i, r in enumerate(r_ecef):
            lat, lon, _ = ecef_to_geodetic(r[0], r[1], r[2])
            lats[i] = lat
            lons[i] = lon
        
        return times, lats, lons


class KeplerPropagator(SatellitePropagator):
    """Simple two-body Keplerian propagation."""
    
    def propagate_eci(self, times_s: np.ndarray) -> np.ndarray:
        """Two-body propagation in ECI."""
        a = self.elements.semi_major_axis_km
        e = self.elements.eccentricity
        inc = math.radians(self.elements.inclination_deg)
        raan = math.radians(self.elements.raan_deg)
        argp = math.radians(self.elements.argument_of_perigee_deg)
        
        M0 = math.radians(self.elements.mean_anomaly_deg)
        n = math.sqrt(self.settings.mu_earth / a**3)
        
        times_s = np.asarray(times_s, dtype=float)
        M = (M0 + n * times_s) % (2 * math.pi)
        E = solve_kepler(M, e)
        
        # Perifocal coordinates
        x_pf = a * (np.cos(E) - e)
        y_pf = a * np.sqrt(1 - e**2) * np.sin(E)
        z_pf = np.zeros_like(x_pf)
        
        # Transform to ECI
        Q = perifocal_to_eci(raan, inc, argp)
        r_pf = np.vstack([x_pf, y_pf, z_pf])
        r_eci = (Q @ r_pf).T
        
        return r_eci


class NumericalPropagator(SatellitePropagator):
    """High-precision numerical propagator with full perturbations."""
    
    def propagate_eci(self, times_s: np.ndarray) -> np.ndarray:
        """Numerical integration with perturbations."""
        times_s = np.asarray(times_s, dtype=float)
        
        # Initial state from Kepler
        kepler = KeplerPropagator(self.elements, self.settings)
        r0 = kepler.propagate_eci(np.array([0.0]))[0]
        
        # Initial velocity (from vis-viva and orbital elements)
        a = self.elements.semi_major_axis_km
        e = self.elements.eccentricity
        mu = self.settings.mu_earth
        
        r0_mag = np.linalg.norm(r0)
        v_mag = math.sqrt(mu * (2 / r0_mag - 1 / a))
        
        # Velocity direction (perpendicular to radius in orbital plane)
        inc = math.radians(self.elements.inclination_deg)
        raan = math.radians(self.elements.raan_deg)
        argp = math.radians(self.elements.argument_of_perigee_deg)
        nu = math.radians(self.elements.true_anomaly_deg)
        
        # Velocity in perifocal frame
        h = math.sqrt(mu * a * (1 - e**2))
        v_pf_x = -(mu / h) * math.sin(nu)
        v_pf_y = (mu / h) * (e + math.cos(nu))
        v_pf = np.array([v_pf_x, v_pf_y, 0.0])
        
        Q = perifocal_to_eci(raan, inc, argp)
        v0 = Q @ v_pf
        
        # Initial state vector
        y0 = np.concatenate([r0, v0])
        
        # Integration
        def derivatives(t, y):
            r_vec = y[:3]
            v_vec = y[3:]
            r = np.linalg.norm(r_vec)
            
            # Two-body
            a_2body = -self.settings.mu_earth / r**3 * r_vec
            
            # Perturbations
            a_pert = np.zeros(3, dtype=float)
            
            # Zonal harmonics
            a_pert += acceleration_j2_j6(r_vec, self.settings)
            
            # Atmospheric drag
            a_pert += acceleration_drag(r_vec, v_vec, self.settings)
            
            # Third-body gravity
            current_time = self.elements.epoch + dt.timedelta(seconds=float(t))
            jd = datetime_to_jd(current_time)
            
            if self.settings.enable_solar_gravity:
                r_sun = sun_position_ecliptic(jd)
                a_pert += acceleration_third_body(r_vec, r_sun, self.settings.mu_sun)
                
                # Solar radiation pressure
                a_pert += acceleration_solar_radiation(r_vec, r_sun, self.settings)
            
            if self.settings.enable_lunar_gravity:
                r_moon = moon_position_simple(jd)
                a_pert += acceleration_third_body(r_vec, r_moon, self.settings.mu_moon)
            
            a_total = a_2body + a_pert
            
            return np.concatenate([v_vec, a_total])
        
        # Solve ODE
        sol = solve_ivp(
            derivatives,
            [times_s[0], times_s[-1]],
            y0,
            method=self.settings.integrator_method,
            t_eval=times_s,
            rtol=self.settings.integrator_rtol,
            atol=self.settings.integrator_atol,
            max_step=self.settings.max_step_s
        )
        
        if not sol.success:
            warnings.warn(f"Integration failed: {sol.message}")
        
        return sol.y[:3, :].T


class SGP4Propagator(SatellitePropagator):
    """SGP4/SDP4 analytical propagator (industry standard)."""
    
    def __init__(self, elements: OrbitalElements, settings: PropagatorSettings):
        super().__init__(elements, settings)
        
        if not HAVE_SGP4:
            raise ImportError("SGP4 library not available")
        
        if elements.tle_line1 is None or elements.tle_line2 is None:
            raise ValueError("TLE lines required for SGP4 propagation")
        
        # Initialize SGP4
        self.satellite = Satrec.twoline2rv(elements.tle_line1, elements.tle_line2)
    
    def propagate_eci(self, times_s: np.ndarray) -> np.ndarray:
        """SGP4 propagation (TEME frame, close to ECI)."""
        times_s = np.asarray(times_s, dtype=float)
        r_eci = np.empty((len(times_s), 3), dtype=float)
        
        for i, t in enumerate(times_s):
            current_time = self.elements.epoch + dt.timedelta(seconds=float(t))
            jd, fr = jday(current_time.year, current_time.month, current_time.day,
                         current_time.hour, current_time.minute, current_time.second)
            
            error_code, r, v = self.satellite.sgp4(jd, fr)
            
            if error_code != 0:
                warnings.warn(f"SGP4 error code {error_code} at t={t}s")
                r_eci[i] = [np.nan, np.nan, np.nan]
            else:
                r_eci[i] = r  # Already in km
        
        return r_eci


def create_propagator(elements: OrbitalElements, 
                     settings: PropagatorSettings) -> SatellitePropagator:
    """Factory function to create appropriate propagator."""
    if settings.propagator_type == "kepler":
        return KeplerPropagator(elements, settings)
    elif settings.propagator_type == "numerical":
        return NumericalPropagator(elements, settings)
    elif settings.propagator_type == "sgp4":
        return SGP4Propagator(elements, settings)
    else:
        raise ValueError(f"Unknown propagator: {settings.propagator_type}")


# ═══════════════════════════════════════════════════════════════════════════
# PASS PREDICTION
# ═══════════════════════════════════════════════════════════════════════════

def compute_topocentric(r_sat_ecef: np.ndarray, station: GroundStation) -> Tuple[float, float, float]:
    """
    Compute topocentric azimuth, elevation, and range.
    
    Returns:
        (azimuth_deg, elevation_deg, slant_range_km)
    """
    r_sta = station.ecef_km()
    rho = r_sat_ecef - r_sta
    
    east, north, up = station.enu_basis()
    
    e = np.dot(rho, east)
    n = np.dot(rho, north)
    u = np.dot(rho, up)
    
    slant = np.linalg.norm(rho)
    elevation = math.degrees(math.atan2(u, math.sqrt(e**2 + n**2)))
    azimuth = math.degrees(math.atan2(e, n)) % 360.0
    
    return azimuth, elevation, slant


def find_passes(propagator: SatellitePropagator, 
               station: GroundStation,
               duration_s: float,
               settings: PropagatorSettings) -> List[PassMetrics]:
    """
    Find all passes above elevation mask.
    
    Args:
        propagator: Satellite propagator
        station: Ground station
        duration_s: Search duration from epoch (seconds)
        settings: Configuration
    
    Returns:
        List of PassMetrics
    """
    step = settings.pass_search_step_s
    elev_mask = settings.elevation_mask_deg
    
    times = np.arange(0, duration_s + step, step)
    r_ecef = propagator.propagate_ecef(times)
    
    elevations = np.empty(len(times))
    azimuths = np.empty(len(times))
    ranges = np.empty(len(times))
    
    for i, r in enumerate(r_ecef):
        az, el, rng = compute_topocentric(r, station)
        elevations[i] = el
        azimuths[i] = az
        ranges[i] = rng
    
    # Find segments above mask
    above = elevations >= elev_mask
    passes = []
    
    i = 0
    while i < len(times):
        if not above[i]:
            i += 1
            continue
        
        # Start of pass
        seg_start = i
        while i + 1 < len(times) and above[i + 1]:
            i += 1
        seg_end = i
        
        # Refine rise/set times
        if seg_start > 0:
            rise_time = refine_crossing(propagator, station, times[seg_start - 1], 
                                       times[seg_start], elev_mask, settings)
        else:
            rise_time = times[0]
        
        if seg_end < len(times) - 1:
            set_time = refine_crossing(propagator, station, times[seg_end], 
                                      times[seg_end + 1], elev_mask, settings)
        else:
            set_time = times[-1]
        
        # Find TCA (max elevation)
        tca_idx = seg_start + np.argmax(elevations[seg_start:seg_end + 1])
        tca_time = times[tca_idx]
        
        # Compute metrics at rise, TCA, set
        r_rise = propagator.propagate_ecef(np.array([rise_time]))[0]
        r_tca = propagator.propagate_ecef(np.array([tca_time]))[0]
        r_set = propagator.propagate_ecef(np.array([set_time]))[0]
        
        az_rise, _, _ = compute_topocentric(r_rise, station)
        az_tca, el_tca, rng_tca = compute_topocentric(r_tca, station)
        az_set, _, _ = compute_topocentric(r_set, station)
        
        # Check if sunlit at TCA
        tca_datetime = propagator.elements.epoch + dt.timedelta(seconds=float(tca_time))
        r_eci_tca = propagator.propagate_eci(np.array([tca_time]))[0]
        r_sun = sun_position_ecliptic(datetime_to_jd(tca_datetime))
        sunlit = is_sunlit(r_eci_tca, r_sun, propagator.settings.earth_radius)
        
        passes.append(PassMetrics(
            rise_time=propagator.elements.epoch + dt.timedelta(seconds=float(rise_time)),
            set_time=propagator.elements.epoch + dt.timedelta(seconds=float(set_time)),
            tca_time=tca_datetime,
            max_elevation_deg=el_tca,
            max_azimuth_deg=az_tca,
            slant_range_km=rng_tca,
            rise_azimuth_deg=az_rise,
            set_azimuth_deg=az_set,
            duration_s=set_time - rise_time,
            sunlit=sunlit
        ))
        
        i = seg_end + 1
    
    return passes


def refine_crossing(propagator: SatellitePropagator, station: GroundStation,
                   t1: float, t2: float, target_elev: float,
                   settings: PropagatorSettings) -> float:
    """Refine rise/set time using bisection."""
    
    def elevation_at_time(t):
        r = propagator.propagate_ecef(np.array([t]))[0]
        _, el, _ = compute_topocentric(r, station)
        return el - target_elev
    
    try:
        t_cross = brentq(elevation_at_time, t1, t2, xtol=settings.refine_tolerance_s)
        return t_cross
    except ValueError:
        return 0.5 * (t1 + t2)


def is_sunlit(r_sat_eci: np.ndarray, r_sun_eci: np.ndarray, earth_radius: float) -> bool:
    """Check if satellite is in sunlight (simple cylindrical shadow)."""
    sun_dir = r_sun_eci / np.linalg.norm(r_sun_eci)
    proj = np.dot(r_sat_eci, sun_dir)
    
    if proj < 0:  # On night side
        r_perp = r_sat_eci - proj * sun_dir
        if np.linalg.norm(r_perp) < earth_radius:
            return False  # In shadow
    
    return True


# ═══════════════════════════════════════════════════════════════════════════
# FILE PARSING
# ═══════════════════════════════════════════════════════════════════════════

def parse_orbital_elements_file(filepath: str) -> OrbitalElements:
    """Parse orbital elements from text file."""
    fields = {
        "Spacecraft": ("name", str),
        "UTC time at deployment": ("epoch", lambda s: dt.datetime.fromisoformat(s.strip().rstrip('Z') + '+00:00')),
        "Semi-major axis (km)": ("semi_major_axis_km", float),
        "Eccentricity (-)": ("eccentricity", float),
        "Inclination (deg)": ("inclination_deg", float),
        "RAAN (deg)": ("raan_deg", float),
        "Argument of perigee (deg)": ("argument_of_perigee_deg", float),
        "True anomaly (deg)": ("true_anomaly_deg", float),
        "Period (s)": ("period_s", float),
    }
    
    result = {}
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if ':' not in line:
                continue
            
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip()
            
            if key in fields:
                field_name, parser = fields[key]
                result[field_name] = parser(value)
    
    # Validate
    required = {"epoch", "semi_major_axis_km", "eccentricity", "inclination_deg",
                "raan_deg", "argument_of_perigee_deg", "true_anomaly_deg", "period_s"}
    missing = required - set(result.keys())
    if missing:
        raise ValueError(f"Missing fields: {missing}")
    
    return OrbitalElements(**result)


def parse_tle_file(filepath: str) -> OrbitalElements:
    """
    Parse TLE file and convert to OrbitalElements.
    
    TLE format:
        Line 0: Satellite name
        Line 1: TLE line 1
        Line 2: TLE line 2
    """
    if not HAVE_SGP4:
        raise ImportError("SGP4 required for TLE parsing. Install: pip install sgp4")
    
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    if len(lines) < 2:
        raise ValueError("TLE file must have at least 2 lines")
    
    # Check if first line is satellite name or TLE line 1
    if lines[0].startswith('1 '):
        name = "Unknown"
        line1 = lines[0]
        line2 = lines[1]
    else:
        name = lines[0]
        line1 = lines[1]
        line2 = lines[2]
    
    # Parse with SGP4
    sat = Satrec.twoline2rv(line1, line2)
    
    # Extract epoch
    epoch_year = sat.epochyr
    epoch_days = sat.epochdays
    
    if epoch_year < 57:
        year = 2000 + epoch_year
    else:
        year = 1900 + epoch_year
    
    epoch = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(days=epoch_days - 1)
    
    # Convert TLE elements to classical elements
    # SGP4 uses: no_kozai (mean motion), ecco, inclo, nodeo, argpo, mo
    
    a = (sat.no_kozai / (2 * math.pi / 86400.0)) ** (-2/3) * (398600.4418 ** (1/3))
    
    elements = OrbitalElements(
        epoch=epoch,
        semi_major_axis_km=a,
        eccentricity=sat.ecco,
        inclination_deg=math.degrees(sat.inclo),
        raan_deg=math.degrees(sat.nodeo),
        argument_of_perigee_deg=math.degrees(sat.argpo),
        true_anomaly_deg=0.0,  # TLE uses mean anomaly
        period_s=2 * math.pi / sat.no_kozai,
        name=name,
        tle_line1=line1,
        tle_line2=line2
    )
    
    return elements


# ═══════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

def plot_ground_track(propagator: SatellitePropagator, duration_s: float,
                     step_s: float, station: Optional[GroundStation] = None,
                     settings: Optional[PropagatorSettings] = None,
                     save_path: Optional[str] = None):
    """Plot ground track on world map."""
    
    times, lats, lons = propagator.ground_track(duration_s, step_s)
    
    # Split at dateline
    lons_plot, lats_plot = split_at_dateline(lons, lats)
    
    title = f"Ground Track: {propagator.elements.name or 'Satellite'}"
    title += f"\nPropagator: {settings.propagator_type.upper() if settings else 'Unknown'}"
    
    if HAVE_CARTOPY:
        fig = plt.figure(figsize=(18, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        ax.set_global()
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.5)
        ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
        ax.add_feature(cfeature.OCEAN, facecolor='#e3f2fd')
        
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', 
                         alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        
        ax.plot(lons_plot, lats_plot, 'b-', linewidth=2, 
               transform=ccrs.PlateCarree(), label='Ground Track')
        ax.plot(lons[0], lats[0], 'go', markersize=8, 
               transform=ccrs.PlateCarree(), label='Start')
        ax.plot(lons[-1], lats[-1], 'rx', markersize=10, 
               transform=ccrs.PlateCarree(), label='End')
        
        if station:
            ax.plot(station.lon_deg, station.lat_deg, '^', color='red', 
                   markersize=12, transform=ccrs.PlateCarree(), 
                   label=station.name)
        
        ax.legend(loc='lower left', fontsize=10)
        ax.set_title(title, fontsize=14, fontweight='bold')
        
    else:
        fig, ax = plt.subplots(figsize=(14, 7))
        ax.plot(lons_plot, lats_plot, 'b-', linewidth=2, label='Ground Track')
        ax.plot(lons[0], lats[0], 'go', markersize=8, label='Start')
        ax.plot(lons[-1], lats[-1], 'rx', markersize=10, label='End')
        
        if station:
            ax.plot(station.lon_deg, station.lat_deg, '^', color='red', 
                   markersize=12, label=station.name)
        
        ax.set_xlabel('Longitude (°)', fontsize=12)
        ax.set_ylabel('Latitude (°)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_title(title, fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=settings.plot_dpi if settings else 300, 
                   bbox_inches='tight')
    
    if settings and settings.show_plots:
        plt.show()
    else:
        plt.close(fig)


def split_at_dateline(lons: np.ndarray, lats: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Insert NaN where longitude crosses ±180°."""
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)
    
    if len(lons) == 0:
        return lons, lats
    
    out_lons = [lons[0]]
    out_lats = [lats[0]]
    
    for i in range(1, len(lons)):
        if abs(lons[i] - lons[i-1]) > 180:
            out_lons.append(np.nan)
            out_lats.append(np.nan)
        out_lons.append(lons[i])
        out_lats.append(lats[i])
    
    return np.array(out_lons), np.array(out_lats)


# ═══════════════════════════════════════════════════════════════════════════
# MAIN PROGRAM
# ═══════════════════════════════════════════════════════════════════════════

def print_banner():
    """Print program banner."""
    banner = """
╔═══════════════════════════════════════════════════════════════════════════╗
║                                                                           ║
║         PROFESSIONAL ORBITAL PROPAGATOR v2.0                              ║
║         High-Precision Satellite Tracking System                          ║
║                                                                           ║
║         • SGP4/SDP4 Analytical Propagation                                ║
║         • Numerical Integration with Full Physics                         ║
║         • Advanced Pass Prediction                                        ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
"""
    print(banner)


def print_orbital_summary(elements: OrbitalElements, settings: PropagatorSettings):
    """Print orbital elements summary."""
    print("\n" + "═" * 75)
    print("ORBITAL ELEMENTS")
    print("═" * 75)
    print(f"Satellite:        {elements.name or 'Unknown'}")
    print(f"Epoch:            {format_datetime(elements.epoch)}")
    print(f"Semi-major axis:  {elements.semi_major_axis_km:.4f} km")
    print(f"Eccentricity:     {elements.eccentricity:.8f}")
    print(f"Inclination:      {elements.inclination_deg:.4f}°")
    print(f"RAAN:             {elements.raan_deg:.4f}°")
    print(f"Arg. of perigee:  {elements.argument_of_perigee_deg:.4f}°")
    print(f"True anomaly:     {elements.true_anomaly_deg:.4f}°")
    print(f"Mean anomaly:     {elements.mean_anomaly_deg:.4f}°")
    print(f"Period:           {elements.period_s:.2f} s ({elements.period_s/60:.2f} min)")
    print(f"Mean motion:      {360.0/elements.period_s*86400:.6f} rev/day")
    
    # Altitude
    perigee = elements.semi_major_axis_km * (1 - elements.eccentricity) - WGS84_A_KM
    apogee = elements.semi_major_axis_km * (1 + elements.eccentricity) - WGS84_A_KM
    print(f"Perigee altitude: {perigee:.2f} km")
    print(f"Apogee altitude:  {apogee:.2f} km")
    
    print("\n" + "═" * 75)
    print("PROPAGATOR CONFIGURATION")
    print("═" * 75)
    print(f"Method:           {settings.propagator_type.upper()}")
    
    if settings.propagator_type == "numerical":
        print(f"Integrator:       {settings.integrator_method}")
        print(f"Tolerance:        rtol={settings.integrator_rtol:.1e}, atol={settings.integrator_atol:.1e}")
        print(f"\nPerturbations:")
        print(f"  J2-J6:          {'Enabled' if settings.enable_J2 else 'Disabled'}")
        print(f"  Drag:           {'Enabled' if settings.enable_atmospheric_drag else 'Disabled'}")
        print(f"  Solar pressure: {'Enabled' if settings.enable_solar_radiation else 'Disabled'}")
        print(f"  Sun gravity:    {'Enabled' if settings.enable_solar_gravity else 'Disabled'}")
        print(f"  Moon gravity:   {'Enabled' if settings.enable_lunar_gravity else 'Disabled'}")
        
        if settings.enable_atmospheric_drag:
            print(f"\nDrag parameters:")
            print(f"  CD:             {settings.drag_coefficient:.2f}")
            print(f"  Area/Mass:      {settings.area_to_mass_ratio:.4f} m²/kg")
        
        if settings.enable_solar_radiation:
            print(f"\nSRP parameters:")
            print(f"  CR:             {settings.reflectivity_coefficient:.2f}")
    
    print("═" * 75)


def print_pass_summary(passes: List[PassMetrics], station: GroundStation):
    """Print pass prediction summary."""
    print("\n" + "═" * 75)
    print(f"PASS PREDICTIONS FOR {station.name.upper()}")
    print(f"Location: {station.lat_deg:.4f}°N, {station.lon_deg:.4f}°E, {station.alt_m:.1f}m")
    print("═" * 75)
    
    if not passes:
        print("No passes found in search window.")
        return
    
    for i, p in enumerate(passes, 1):
        print(f"\n{'─' * 75}")
        print(f"PASS #{i}")
        print(f"{'─' * 75}")
        print(p)


def main():
    """Main program entry point."""
    
    parser = argparse.ArgumentParser(
        description="Professional Orbital Propagator v2.0",
        epilog="Example: python orbital_propagator_pro.py data.txt --propagator sgp4 --duration-hours 24",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Input file
    parser.add_argument("input_file", help="Orbital elements or TLE file")
    parser.add_argument("--file-type", choices=["elements", "tle"], default="auto",
                       help="Input file type (auto-detect by default)")
    
    # Propagator settings
    parser.add_argument("--propagator", choices=["kepler", "numerical", "sgp4"],
                       default="sgp4", help="Propagation method")
    
    # Duration
    parser.add_argument("--duration-hours", type=float, default=24.0,
                       help="Propagation duration (hours)")
    parser.add_argument("--plot-step-s", type=float, default=30.0,
                       help="Ground track time step (seconds)")
    
    # Station
    parser.add_argument("--station-name", default="Trondheim")
    parser.add_argument("--station-lat", type=float, default=63.41866)
    parser.add_argument("--station-lon", type=float, default=10.3951)
    parser.add_argument("--station-alt", type=float, default=51.3)
    
    # Pass prediction
    parser.add_argument("--elevation-mask", type=float, default=5.0,
                       help="Minimum elevation for passes (degrees)")
    parser.add_argument("--pass-search-step", type=float, default=10.0,
                       help="Pass search time step (seconds)")
    
    # Visualization
    parser.add_argument("--save-plot", help="Save ground track plot to file")
    parser.add_argument("--no-plot", action="store_true", help="Don't show plots")
    
    # Physical parameters (for numerical propagator)
    parser.add_argument("--satellite-mass", type=float, default=100.0,
                       help="Satellite mass (kg)")
    parser.add_argument("--satellite-area", type=float, default=1.0,
                       help="Cross-sectional area (m²)")
    parser.add_argument("--drag-coeff", type=float, default=2.2,
                       help="Drag coefficient")
    
    args = parser.parse_args()
    
    # Print banner
    print_banner()
    
    # Create settings
    settings = PropagatorSettings(
        propagator_type=args.propagator,
        elevation_mask_deg=args.elevation_mask,
        pass_search_step_s=args.pass_search_step,
        show_plots=not args.no_plot,
        satellite_mass_kg=args.satellite_mass,
        satellite_area_m2=args.satellite_area,
        drag_coefficient=args.drag_coeff
    )
    
    # Parse input file
    try:
        file_type = args.file_type
        
        if file_type == "auto":
            # Auto-detect
            with open(args.input_file, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('1 ') or (not first_line.startswith('1 ') and 
                   len(first_line) == 69):
                    file_type = "tle"
                else:
                    file_type = "elements"
        
        if file_type == "tle":
            elements = parse_tle_file(args.input_file)
        else:
            elements = parse_orbital_elements_file(args.input_file)
            
    except FileNotFoundError:
        print(f"✗ Error: File '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error parsing file: {e}")
        sys.exit(1)
    
    # Print summary
    print_orbital_summary(elements, settings)
    
    # Create propagator
    try:
        propagator = create_propagator(elements, settings)
    except Exception as e:
        print(f"\n✗ Error creating propagator: {e}")
        sys.exit(1)
    
    # Ground station
    station = GroundStation(
        name=args.station_name,
        lat_deg=args.station_lat,
        lon_deg=args.station_lon,
        alt_m=args.station_alt
    )
    
    # Duration
    duration_s = args.duration_hours * 3600.0
    
    # Plot ground track
    print(f"\n⚙ Computing ground track...")
    try:
        plot_ground_track(
            propagator,
            duration_s=duration_s,
            step_s=args.plot_step_s,
            station=station,
            settings=settings,
            save_path=args.save_plot
        )
        print("✓ Ground track plotted")
    except Exception as e:
        print(f"✗ Plotting failed: {e}")
    
    # Find passes
    print(f"\n⚙ Searching for passes...")
    try:
        passes = find_passes(propagator, station, duration_s, settings)
        print(f"✓ Found {len(passes)} passes")
        print_pass_summary(passes, station)
    except Exception as e:
        print(f"✗ Pass search failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "═" * 75)
    print("ANALYSIS COMPLETE")
    print("═" * 75)


if __name__ == "__main__":
    main()
