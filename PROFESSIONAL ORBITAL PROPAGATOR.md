# 📘 PROFESSIONAL ORBITAL PROPAGATOR - USER MANUAL
## Version 2.0 - Complete Guide

---

## 📋 TABLE OF CONTENTS

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Configuration Guide](#configuration-guide)
4. [Usage Examples](#usage-examples)
5. [Propagator Types](#propagator-types)
6. [Input File Formats](#input-file-formats)
7. [Advanced Features](#advanced-features)
8. [Troubleshooting](#troubleshooting)
9. [Technical Reference](#technical-reference)

---

## 🚀 QUICK START

### Minimal Example (5 minutes)

```bash
# 1. Install dependencies
pip install numpy scipy matplotlib sgp4

# Optional but recommended:
conda install -c conda-forge cartopy

# 2. Run with your test file
python orbital_propagator_pro.py orbital_elements_test.txt

# 3. Done! You'll see:
#    - Orbital summary
#    - Ground track map
#    - Pass predictions
```

---

## 💿 INSTALLATION

### Required Dependencies

```bash
# Core requirements (mandatory)
pip install numpy scipy matplotlib

# SGP4 support (highly recommended)
pip install sgp4

# Beautiful maps (optional but nice)
conda install -c conda-forge cartopy
# OR
pip install cartopy  # May need system dependencies
```

### System Requirements

- **Python**: 3.8 or newer
- **RAM**: 4 GB minimum, 8 GB recommended
- **CPU**: Any modern processor (numerical propagation benefits from faster CPUs)

---

## ⚙️ CONFIGURATION GUIDE

### The Settings Class

At the top of the script, you'll find `PropagatorSettings`. This is your control panel:

```python
@dataclass
class PropagatorSettings:
    # ═══ PROPAGATION METHOD ═══
    propagator_type: str = "sgp4"  # ← CHANGE THIS
    
    # Options:
    # "sgp4"      - Industry standard, best for TLE files (±1 km accuracy)
    # "numerical" - High-precision physics simulation (±100m accuracy)
    # "kepler"    - Simple two-body (±5 km accuracy, fast)
```

### Key Settings Explained

#### 1. **Propagator Type**
```python
propagator_type: str = "sgp4"
```
- **sgp4**: Use this for real TLE files from space-track.org
- **numerical**: Use for maximum accuracy with custom physics
- **kepler**: Use for educational purposes or quick estimates

#### 2. **Numerical Integrator** (only for `propagator_type="numerical"`)
```python
integrator_method: str = "DOP853"     # Very accurate, slower
# Options: "RK45" (faster), "DOP853" (best), "Radau" (stiff problems)

integrator_rtol: float = 1e-11        # Relative tolerance
integrator_atol: float = 1e-12        # Absolute tolerance
```

**When to adjust:**
- Increase tolerance (e.g., 1e-8) for faster computation
- Decrease tolerance (e.g., 1e-13) for maximum precision

#### 3. **Perturbation Forces** (only for `propagator_type="numerical"`)
```python
enable_J2: bool = True                # Earth's bulge (ESSENTIAL)
enable_J3: bool = True                # Higher gravity terms
enable_atmospheric_drag: bool = True  # Drag (important below 800 km)
enable_solar_radiation: bool = True   # Light pressure (important for large satellites)
enable_lunar_gravity: bool = True     # Moon's pull
enable_solar_gravity: bool = True     # Sun's pull
```

**Recommendations:**
- **Always** keep J2 enabled
- For LEO (<800 km): Enable drag
- For GEO/MEO: Disable drag, enable SRP
- For quick tests: Disable J3-J6, moon, sun

#### 4. **Satellite Physical Properties**
```python
satellite_mass_kg: float = 100.0      # Your satellite mass
satellite_area_m2: float = 1.0        # Cross-sectional area
drag_coefficient: float = 2.2         # Typical: 2.0-2.5
reflectivity_coefficient: float = 1.3 # Typical: 1.0-2.0
```

**How to set these:**
- **CubeSat (1U)**: mass=1.33 kg, area=0.01 m²
- **CubeSat (3U)**: mass=4 kg, area=0.03 m²
- **Small sat**: mass=100 kg, area=1 m²
- **ISS**: mass=420000 kg, area=~2500 m²

#### 5. **Pass Prediction**
```python
elevation_mask_deg: float = 5.0       # Don't show passes below this
pass_search_step_s: float = 10.0      # Smaller = more accurate but slower
```

---

## 📖 USAGE EXAMPLES

### Example 1: Basic Usage (Orbital Elements File)

```bash
python orbital_propagator_pro.py orbital_elements_test.txt \
    --propagator sgp4 \
    --duration-hours 24 \
    --station-lat 63.41866 \
    --station-lon 10.3951
```

**Output:**
- Orbital summary
- 24-hour ground track
- All passes over Trondheim

### Example 2: High-Precision Numerical Propagation

```bash
python orbital_propagator_pro.py orbital_elements_test.txt \
    --propagator numerical \
    --duration-hours 6 \
    --satellite-mass 100 \
    --satellite-area 1.0 \
    --drag-coeff 2.2
```

**When to use:**
- You need maximum accuracy
- You're willing to wait longer
- You have satellite physical properties

### Example 3: TLE File (Real Satellite)

```bash
# First, get TLE from https://celestrak.org or space-track.org
# Save as iss.tle

python orbital_propagator_pro.py iss.tle \
    --file-type tle \
    --propagator sgp4 \
    --duration-hours 48
```

### Example 4: Save Plot Without Showing

```bash
python orbital_propagator_pro.py data.txt \
    --save-plot my_groundtrack.png \
    --no-plot
```

### Example 5: Custom Ground Station

```bash
python orbital_propagator_pro.py data.txt \
    --station-name "My Observatory" \
    --station-lat 40.7128 \
    --station-lon -74.0060 \
    --station-alt 10.0 \
    --elevation-mask 10.0
```

---

## 🔬 PROPAGATOR TYPES

### SGP4/SDP4 (Recommended for Real Operations)

**Pros:**
- ✅ Industry standard (NASA, Space Force use this)
- ✅ Fast
- ✅ Works with public TLE data
- ✅ Includes atmospheric drag automatically
- ✅ Accuracy: ±1 km after 1 day

**Cons:**
- ❌ Requires TLE files (must update every few days)
- ❌ Less accurate than numerical for long-term

**Best for:**
- Real satellite tracking
- Pass predictions
- Amateur radio satellite tracking
- Operations (if you update TLEs regularly)

### Numerical Integration (Highest Precision)

**Pros:**
- ✅ Most accurate (±100m - 1 km)
- ✅ Full physics simulation
- ✅ Customizable force models
- ✅ Works with any orbital elements

**Cons:**
- ❌ Slower (10-100x slower than SGP4)
- ❌ Requires physical properties (mass, area)
- ❌ More complex to configure

**Best for:**
- Mission planning
- Research
- When you need maximum accuracy
- Long-term predictions (weeks/months)

### Kepler (Simple Two-Body)

**Pros:**
- ✅ Very fast
- ✅ Simple, easy to understand
- ✅ No external dependencies

**Cons:**
- ❌ Inaccurate (±5-50 km after 1 day)
- ❌ No perturbations
- ❌ Not suitable for real operations

**Best for:**
- Quick estimates
- Education
- Initial orbit design
- Debugging

---

## 📄 INPUT FILE FORMATS

### Format 1: Orbital Elements (Custom Format)

```text
Spacecraft: ISS-TEST
UTC time at deployment: 2025-04-10T12:00:00Z
Semi-major axis (km): 6778.137
Eccentricity (-): 0.0001
Inclination (deg): 51.6400
RAAN (deg): 125.0000
Argument of perigee (deg): 90.0000
True anomaly (deg): 0.0000
Period (s): 5558.0
```

**Use with:** Any propagator (`--propagator kepler/numerical/sgp4`)

### Format 2: TLE (Two-Line Elements)

```text
ISS (ZARYA)
1 25544U 98067A   25100.50000000  .00002182  00000-0  41420-4 0  9990
2 25544  51.6400 125.0000 0001000  90.0000   0.0000 15.48919293123456
```

**Use with:** SGP4 only (`--propagator sgp4 --file-type tle`)

**Where to get TLEs:**
- https://celestrak.org (no login required)
- https://space-track.org (free account required)
- https://n2yo.com (some satellites)

---

## 🎓 ADVANCED FEATURES

### 1. Customize Physics Models

Edit the `PropagatorSettings` class directly in the script:

```python
# Example: Disable solar/lunar gravity for faster computation
settings = PropagatorSettings(
    propagator_type="numerical",
    enable_J2=True,
    enable_J3=True,
    enable_atmospheric_drag=True,
    enable_solar_radiation=False,  # ← CHANGED
    enable_lunar_gravity=False,    # ← CHANGED
    enable_solar_gravity=False     # ← CHANGED
)
```

### 2. Batch Processing Multiple Satellites

```python
import glob

tle_files = glob.glob("*.tle")

for tle_file in tle_files:
    os.system(f"python orbital_propagator_pro.py {tle_file} --no-plot --save-plot {tle_file}.png")
```

### 3. Export Pass Predictions to CSV

Add this to the script:

```python
import csv

# After finding passes:
with open('passes.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Pass', 'Rise Time', 'TCA', 'Set Time', 'Max El', 'Duration'])
    
    for i, p in enumerate(passes, 1):
        writer.writerow([
            i,
            p.rise_time.isoformat(),
            p.tca_time.isoformat(),
            p.set_time.isoformat(),
            f"{p.max_elevation_deg:.2f}",
            f"{p.duration_s:.1f}"
        ])
```

### 4. API Usage (Use as Library)

```python
from orbital_propagator_pro import *

# Create elements
elements = OrbitalElements(
    epoch=dt.datetime.now(dt.timezone.utc),
    semi_major_axis_km=6778.0,
    eccentricity=0.001,
    inclination_deg=51.6,
    raan_deg=0.0,
    argument_of_perigee_deg=0.0,
    true_anomaly_deg=0.0,
    period_s=5558.0
)

# Configure
settings = PropagatorSettings(propagator_type="numerical")

# Propagate
propagator = create_propagator(elements, settings)
times = np.array([0, 60, 120, 180])  # seconds
positions = propagator.propagate_eci(times)

print(positions)  # Shape: (4, 3) - four positions in ECI
```

---

## 🔧 TROUBLESHOOTING

### Issue: "SGP4 library not available"

```bash
pip install sgp4
```

### Issue: "Cartopy not installed" warning

**Option 1 (recommended):**
```bash
conda install -c conda-forge cartopy
```

**Option 2 (pip, may need system packages):**
```bash
# Ubuntu/Debian:
sudo apt-get install libgeos-dev libproj-dev
pip install cartopy

# macOS:
brew install geos proj
pip install cartopy
```

**Option 3 (ignore it):**
- The warning is harmless
- Maps will be basic but functional

### Issue: "Integration failed" or NaN results

**Causes:**
1. Orbit is decaying (satellite re-entered)
2. Eccentricity ≥ 1.0 (hyperbolic orbit)
3. Tolerance too strict

**Solutions:**
```python
# Relax tolerances
integrator_rtol: float = 1e-8  # Instead of 1e-11
integrator_atol: float = 1e-9  # Instead of 1e-12

# Check your orbital elements
# Eccentricity must be < 1.0
# Semi-major axis must be > Earth radius
```

### Issue: Very slow numerical propagation

**Solutions:**
1. Use SGP4 instead: `--propagator sgp4`
2. Reduce duration: `--duration-hours 6` instead of 24
3. Increase step size: `--plot-step-s 60` instead of 30
4. Disable expensive perturbations:
```python
enable_lunar_gravity: bool = False
enable_solar_gravity: bool = False
```

### Issue: No passes found

**Causes:**
1. Satellite never comes above your horizon
2. Wrong station coordinates
3. Elevation mask too high

**Solutions:**
```bash
# Lower elevation mask
--elevation-mask 0.0

# Check station coordinates
--station-lat 63.41866 --station-lon 10.3951

# Increase search duration
--duration-hours 48
```

---

## 📚 TECHNICAL REFERENCE

### Accuracy Comparison

| Propagator | 1 hour | 1 day | 1 week | Speed |
|------------|--------|-------|--------|-------|
| Kepler     | ±500m  | ±5km  | ±50km  | Very Fast |
| SGP4       | ±100m  | ±1km  | ±10km  | Fast |
| Numerical (J2 only) | ±100m | ±500m | ±5km | Medium |
| Numerical (full) | ±50m | ±200m | ±2km | Slow |

*Accuracy depends on orbit type and altitude*

### Force Model Details

#### J2 Perturbation
- **Effect**: Causes RAAN and argument of perigee to drift
- **Magnitude**: ~10 km/day for LEO
- **Always enable this**

#### Atmospheric Drag
- **Effect**: Orbital decay, semi-major axis decreases
- **Magnitude**: Depends heavily on altitude and solar activity
  - 300 km: ~1-5 km/day decay
  - 400 km: ~100-500 m/day
  - 800 km: ~1-10 m/day
  - >1000 km: negligible

#### Solar Radiation Pressure
- **Effect**: Eccentricity oscillation, semi-major axis drift
- **Magnitude**: Depends on area/mass ratio
  - Small satellites: ~10-100 m/day
  - Large solar panels: ~1 km/day

### Coordinate Systems

- **ECI (Earth-Centered Inertial)**: Non-rotating, used for physics
- **ECEF (Earth-Centered Earth-Fixed)**: Rotates with Earth, used for ground track
- **Geodetic**: Latitude, longitude, altitude (what you see on maps)

### Time Systems

- **UTC**: Universal Coordinated Time (what clocks show)
- **Julian Date**: Days since noon on January 1, 4713 BC
- **GMST**: Greenwich Mean Sidereal Time (Earth's rotation angle)

---

## 🎯 RECOMMENDED WORKFLOWS

### Workflow 1: Track Real Satellite

```bash
# 1. Get fresh TLE from celestrak.org
wget https://celestrak.org/NORAD/elements/gp.php?CATNR=25544 -O iss.tle

# 2. Run propagator
python orbital_propagator_pro.py iss.tle \
    --file-type tle \
    --propagator sgp4 \
    --duration-hours 24 \
    --station-lat YOUR_LAT \
    --station-lon YOUR_LON

# 3. Check passes, plan observations
```

### Workflow 2: Mission Planning (New Satellite)

```bash
# 1. Design orbit, create elements file
# 2. Run high-precision propagation
python orbital_propagator_pro.py mission.txt \
    --propagator numerical \
    --duration-hours 168 \
    --satellite-mass 100 \
    --satellite-area 1.5

# 3. Analyze ground track coverage
# 4. Iterate on orbit design
```

### Workflow 3: Educational/Learning

```bash
# 1. Use Kepler propagator to understand basics
python orbital_propagator_pro.py orbit.txt \
    --propagator kepler \
    --duration-hours 12

# 2. Compare with numerical
python orbital_propagator_pro.py orbit.txt \
    --propagator numerical \
    --duration-hours 12

# 3. See the difference!
```

---

## 📞 SUPPORT & FEEDBACK

### Common Questions

**Q: Which propagator should I use?**
A: SGP4 for real satellites with TLEs. Numerical for custom orbits or maximum accuracy.

**Q: How often should I update TLEs?**
A: Every 1-7 days, depending on altitude. Lower orbits decay faster.

**Q: Can this predict satellite collisions?**
A: No. Use specialized software (CSPOC, AGI STK) for conjunction analysis.

**Q: Is this suitable for mission-critical operations?**
A: The code is educational/research quality. For operations, use validated tools like GMAT or STK.

---

## 📜 VERSION HISTORY

**v2.0** (Current)
- Full SGP4/SDP4 support
- Numerical integration with J2-J6
- Atmospheric drag (exponential model)
- Solar radiation pressure
- Third-body perturbations (Sun, Moon)
- Enhanced visualization
- Comprehensive pass predictions

**v1.0**
- Basic Kepler propagation
- Simple ground track
- Basic pass prediction

---

## 🙏 CREDITS

**Physics Models:**
- SGP4: Simplified General Perturbations #4 (Hoots & Roehrich, 1980)
- Atmospheric: NRLMSISE-00 simplified
- Gravity: WGS84/EGM96 zonal harmonics

**Libraries:**
- sgp4: Brandon Rhodes
- scipy: SciPy developers
- cartopy: Met Office

---

**END OF MANUAL**

For questions or improvements, refer to the code comments or physics literature.
Good luck with your orbital mechanics! 🚀

