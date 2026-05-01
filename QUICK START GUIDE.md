# ⚡ QUICK START GUIDE
## Get Running in 5 Minutes!

---

## STEP 1: Install Dependencies

```bash
# Minimum (required)
pip install numpy scipy matplotlib sgp4

# Recommended (for beautiful maps)
conda install -c conda-forge cartopy
```

---

## STEP 2: Choose Your Test

### Option A: Test with Orbital Elements (Simplest)

```bash
python orbital_propagator_pro.py orbital_elements_test.txt
```

**What you'll see:**
- Orbital summary with all parameters
- Ground track map showing 16 orbits
- Pass predictions over Trondheim

### Option B: Test with Real ISS Data

```bash
python orbital_propagator_pro.py example_iss.tle --file-type tle
```

---

## STEP 3: Customize for Your Location

```bash
python orbital_propagator_pro.py orbital_elements_test.txt \
    --station-name "My Location" \
    --station-lat YOUR_LATITUDE \
    --station-lon YOUR_LONGITUDE \
    --station-alt YOUR_ALTITUDE_METERS
```

**Example for New York:**
```bash
python orbital_propagator_pro.py orbital_elements_test.txt \
    --station-name "New York" \
    --station-lat 40.7128 \
    --station-lon -74.0060 \
    --station-alt 10
```

---

## STEP 4: Choose Accuracy Level

### Fast & Accurate (Recommended)
```bash
python orbital_propagator_pro.py data.txt --propagator sgp4
```

### Maximum Accuracy (Slower)
```bash
python orbital_propagator_pro.py data.txt --propagator numerical
```

### Quick Estimate (Fastest)
```bash
python orbital_propagator_pro.py data.txt --propagator kepler
```

---

## UNDERSTANDING THE OUTPUT

### 1. Orbital Summary
```
═══════════════════════════════════════════════════════════
ORBITAL ELEMENTS
═══════════════════════════════════════════════════════════
Satellite:        ISS-TEST
Epoch:            2025-04-10 12:00:00 UTC
Semi-major axis:  6778.1370 km       ← Orbit size
Eccentricity:     0.0001000          ← How elliptical (0=circle)
Inclination:      51.6400°           ← Tilt relative to equator
Period:           5558.00 s          ← Time for one orbit
Perigee altitude: 399.46 km          ← Lowest point
Apogee altitude:  400.81 km          ← Highest point
```

### 2. Ground Track Map
- **Blue line**: Path satellite takes over Earth
- **Green circle**: Start position
- **Red X**: End position
- **Red triangle**: Your ground station

### 3. Pass Predictions
```
PASS #1
─────────────────────────────────────────────────────────────
  Rise: 2025-04-10 22:43:01 UTC
    Az: 225.34°                      ← Compass direction at rise
  TCA:  2025-04-10 22:45:10 UTC      ← Time of closest approach
    El: 8.25°                        ← How high above horizon
    Az: 180.45°                      ← Compass direction at TCA
  Set:  2025-04-10 22:47:21 UTC
    Az: 135.67°                      ← Compass direction at set
  Duration: 260.7 s (4.3 min)
  Range: 1587.7 km                   ← Distance to satellite
  Sunlit: Yes                        ← Visible (if optical tracking)
```

**Azimuth Guide:**
- 0° = North
- 90° = East
- 180° = South
- 270° = West

**Elevation Guide:**
- 0° = Horizon
- 45° = Halfway up
- 90° = Directly overhead

---

## COMMON TASKS

### Save the Ground Track Image
```bash
python orbital_propagator_pro.py data.txt --save-plot my_orbit.png --no-plot
```

### Find Passes for Next 48 Hours
```bash
python orbital_propagator_pro.py data.txt --duration-hours 48
```

### Lower Elevation Mask (see more passes)
```bash
python orbital_propagator_pro.py data.txt --elevation-mask 0
```

### Higher Elevation Mask (only best passes)
```bash
python orbital_propagator_pro.py data.txt --elevation-mask 20
```

---

## TROUBLESHOOTING

### "No module named 'sgp4'"
```bash
pip install sgp4
```

### "Cartopy not installed" (WARNING - not critical)
- Maps still work, just less pretty
- Install with: `conda install -c conda-forge cartopy`

### "No passes found"
- Try: `--elevation-mask 0 --duration-hours 48`
- Check your station coordinates are correct

### Script runs but no window appears
- Add `--save-plot output.png` to save instead
- Check if running in headless environment

---

## NEXT STEPS

1. **Read USER_MANUAL.md** for complete documentation
2. **Edit PropagatorSettings** in the script for advanced physics
3. **Get real TLE files** from https://celestrak.org
4. **Experiment!** Change orbits, stations, settings

---

## NEED HELP?

Check the full manual: `USER_MANUAL.md`

**Happy Tracking!** 🛰️

