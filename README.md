# lst_downscale

Land surface temperature (LST) from Sentinel-3 is measured on 1 km pixels, too coarse to see a single
field. This project sharpens it to 10 m with Sentinel-2 as the guide, and checks how much of the added
detail is real against Landsat's 100 m thermal band. Study area: 40 x 40 km south of Lahore, around
the fields of [crop_project](https://github.com/mburhankhan101-glitch/crop_project). 

```powershell
python lst_data.py 2026-03-27 2026-05-06 2026-05-14 2026-06-07 2026-06-23   # download and cache
python run.py 2026-03-27 2026-05-06 2026-05-14 2026-06-07 2026-06-23        # sharpen and score (~30 min)
python run.py --figures 2026-05-06                                          # redraw one figure
```

## Data

All from Microsoft Planetary Computer, on three aligned grids in UTM zone 43N (`lst_data.py`):

| Source | Pixel | Role |
|---|---|---|
| Sentinel-3 SLSTR LST (L2) | 1 km | what is sharpened |
| Sentinel-2 L2A, 10 bands | 10 m (20 m bands resampled) | the guide: what is in each 10 m pixel |
| Landsat 8/9 Collection 2 surface temperature | 100 m (delivered at 30 m) | the check |

Five dates in spring 2026 when all three passed on the same day, with Landsat on track 149, the one
that covers the whole area: 27 March, 6 May, 14 May, 7 June and 23 June. Things that had to be
handled: the Sentinel-3 flag files declare 0 as missing although 0 means "no flag raised"; the
Bayesian cloud flag marked 85% of a clear day as cloudy over hot, bright ground, so only the basic
cloud tests are used; and the area sits where two Landsat scenes of a pass meet, so they are stitched.

## Methods (`sharpen.py`)

- **resample**: enlarge the 1 km cells smoothly. No new information: the line to beat.
- **TsHARP**: one straight line, LST = a + b x NDVI, fitted on the 1 km cells, applied at 10 m.
- **Data mining sharpener (DMS)**: an ensemble of regression trees (a random forest) on all ten bands
  plus NDVI, NDWI and NDMI, fitted per scene on the 80% most uniform 1 km cells, applied at 10 m.

TsHARP and DMS end with **residual correction**: inside each 1 km cell the 10 m values are shifted so
their average equals the measured 1 km value. Two kinds:

- **constant**: one shift per 1 km cell. Exact, but neighbouring cells get different shifts, which
  draws straight 1 km edges into the 10 m map.
- **smooth**: the shifts are treated as values at the cell centres and blended bilinearly; the small
  error this leaves in each cell's mean is smoothed and added again, 8 rounds in all.

## Validation (`run.py`)

- **Test A**: Landsat averaged to 1 km, sharpened back, compared with the Landsat it came from. The
  truth is known, so this measures the method alone.
- **Test B**: real Sentinel-3 sharpened and compared with Landsat of the same day: adds two sensors,
  minutes apart, different viewing angles and algorithms.

Both are scored at 100 m (Landsat's real resolution): RMSE and bias, and the **detail scores**, which
look only at the pattern inside each 1 km cell (each 100 m value minus its 1 km mean). The detail is
what sharpening is meant to add: copying the 1 km value would score r = 0. **Nothing can check the
10 m detail itself.** Two more numbers at 10 m: the **edge ratio** (temperature jumps between
neighbouring pixels across 1 km borders divided by the jumps inside cells; about 1 means no visible
grid) and the **conservation error** (mean gap between each 1 km measurement and the mean of its
10 m pixels).

## Results

Mean over the five dates ([per date](runs/results.csv), [summary](runs/summary.csv)):

| Method | Test A RMSE (deg C) | Test A detail r | Test B detail r | Test B detail RMSE (deg C) | Edge ratio (B) | Conservation error (B, deg C) |
|---|---|---|---|---|---|---|
| resample | 1.85 | 0.30 | 0.22 | 1.84 | 1.03 | 0.16 |
| TsHARP | 1.92 | 0.69 | 0.69 | 1.45 | 1.66 | 0.00 |
| DMS trees, constant correction | 1.40 | 0.77 | 0.73 | **1.30** | 1.75 | 0.00 |
| **DMS trees, smooth correction** | **1.37** | **0.78** | 0.72 | 1.32 | **1.03** | 0.01 |

- **The tree sharpener is best on every date and both tests.** It explains most of the temperature
  pattern inside a 1 km cell (r 0.73-0.77) and cuts the error of that pattern by about 30% against
  resampling (1.30 against 1.84 deg C on real Sentinel-3).
- **Smooth correction removes the 1 km grid without losing accuracy.** With constant correction the
  jumps across 1 km borders are 1.5-1.75 times the jumps inside cells; smooth correction brings that
  to 1.03, the same as the real landscape, while each 1 km mean stays within 0.01 deg C of the
  measurement. The 100 m scores barely move (1.40 to 1.37 deg C in Test A, 1.30 to 1.32 in Test B):
  the change is to how the 10 m map looks, not to its accuracy. Smooth is now the default for the maps.

  ![Constant vs smooth residual correction](runs/correction_explained.png)

  Top left, the idea: a 1 km cell measured at 36 deg C next to one at 33 deg C, and a model that says
  34 in both. Constant correction shifts every pixel of a cell by the same amount (+2, then -1), so two
  pixels 10 m apart end up 3 deg C apart at the border: a cliff. Smooth correction blends the shift
  from one cell centre to the next: a ramp. Bottom, the same on real data (6 May), at the border where
  the two corrections differ most (`python figures.py 2026-05-06`): with constant correction the 10 m
  map changes abruptly on the dashed 1 km line, and the column averages drop by 3.5 deg C right at the
  border; with smooth correction they don't. Top right, the edge ratio: how much bigger the jumps
  between neighbouring pixels are across 1 km borders than inside cells (1.0 = no seams).
- **TsHARP finds the pattern (r 0.69) but overshoots it**: in Test A its error (1.92) is no better than
  doing nothing (1.85). One NDVI line exaggerates the contrast between green and bare fields.
- **Moisture matters more than greenness.** In the trees, the shortwave infrared band B12 and the
  moisture index NDMI are the most important inputs on four of five dates; NDVI comes well after.
  Evaporation cools, and SWIR sees water in soil and leaves; that is why a model with SWIR beats a
  model with NDVI alone.
- **Sentinel-3 and Landsat disagree before any sharpening** ([sensors](runs/sensors.csv)). At 1 km they
  correlate well (r 0.86-0.95), but the offset swings from +0.7 to -4.4 deg C with the time of the
  Sentinel-3 pass: within 5 minutes of Landsat (05:39, 05:41 UTC) it is +0.7 and +0.2 deg C; 25 and
  40 minutes earlier (05:11, 04:56) it is -2.7 and -4.4 deg C, as the ground is still heating up in the
  morning. 14 May (-3.1 deg C at 05:33) does not fit that pattern. This offset dominates Test B's overall
  RMSE (3.2 deg C for every method); sharpening cannot fix it.

![6 May 2026](runs/lst_2026-05-06.png)

The 10 m maps are saved as GeoTIFFs in `runs/` (not in git: 30-60 MB each).

## The 2022 heatwave: 27 May 2022

```powershell
python find_dates.py lahore 2022-03-20 2022-05-30       # rank days by how closely the satellites passed
python lst_data.py 2022-05-27
python run.py --tag heatwave2022 2022-05-27
```

Spring 2022 brought a long heatwave and dry spell to Pakistan. `find_dates.py` searched 20 March to
30 May 2022 over Lahore and Islamabad for the days when Sentinel-3, Sentinel-2 and Landsat passed
closest together ([Lahore](runs/dates_lahore_2022-03-20_2022-05-30.csv),
[Islamabad](runs/dates_islamabad_2022-03-20_2022-05-30.csv)). The best by far: **27 May 2022 over
Lahore**, with Sentinel-3B at 05:35, Sentinel-2A at 05:36 and Landsat 8 at 05:36 UTC, and all three
100% clear over the area. Landsat measured a median ground temperature of 45 deg C, above 50 deg C on
the hottest 2%. Islamabad had no comparable day: its only day with all three had Sentinel-3 32
minutes off.

| 27 May 2022 | Test A RMSE (deg C) | Test A detail r | Test B detail r | Test B detail RMSE (deg C) |
|---|---|---|---|---|
| resample | 2.12 | 0.27 | 0.24 | 2.09 |
| TsHARP | 2.15 | 0.74 | 0.74 | 1.50 |
| DMS trees, smooth | **1.32** | **0.82** | **0.78** | **1.33** |

([results](runs/results_heatwave2022.csv), [summary](runs/summary_heatwave2022.csv))

- The tree sharpener does at least as well as on the 2026 dates: it explains 82% (Test A) and 78% (Test
  B) of the pattern inside each 1 km cell, and cuts the error of that pattern by 36% against resampling
  on real Sentinel-3.
- **Simultaneous passes do not remove the gap between the sensors.** One minute apart, Sentinel-3 was
  still 1.5 deg C cooler than Landsat at 1 km ([sensors](runs/sensors_heatwave2022.csv)). So the time
  between passes explains part of the offset seen in 2026, not all of it: the two sensors and their
  temperature algorithms also differ, here by about 1.5 deg C on a very hot day.

![27 May 2022](runs/lst_2022-05-27.png)

## Limitations

- The 10 m detail is inferred from Sentinel-2, not measured; it is checked only at 100 m.
- One area, six dates (five in 2026, one in 2022), one model per scene. 14 May and 23 June had partly cloudy Sentinel-2
  images and are scored on a quarter of the area.
- The time and angle offset between Sentinel-3 and Landsat is reported, not corrected.

## Roadmap

- [x] Data for five same-day dates; resample, TsHARP and tree sharpener; Tests A and B
- [x] Smooth the residual correction to remove the 1 km block edges (edge ratio 1.75 to 1.03)
- [ ] A U-Net trained on Landsat (1 km input, 100 m answer key) and tested on unseen dates
- [ ] Correct the Sentinel-3 to Landsat offset (pass time, viewing angle)
- [ ] Use it: crop temperature at the head and tail of a canal, as evidence of water shortage
