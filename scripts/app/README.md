# Grid Weights Generator (Shiny App)

This Shiny web application allows users to **intersect NetCDF grid data with HRU (Hydrological Response Unit) polygons** to generate area-weighted grid-cell weights. It is especially useful for hydrological or land-surface model preprocessing, where spatial aggregation of gridded climate or forcing data is required.

---

## 📦 Features

* Upload NetCDF file (`.nc`) with 2D latitude/longitude grids.
* Upload HRU shapefile (`.zip`) with polygon features.
* Select variables and dimensions interactively from NetCDF file.
* Choose the attribute in the shapefile to use as HRU ID.
* Generate:

  * Intersected grid polygons
  * HRU boundaries
  * Grid centroids
  * Area-based weights for each HRU–grid pair
* Download all results in a single ZIP file.
* Optional: visualize the grid and HRU overlay on an interactive Leaflet map.
* Automatically checks for and installs missing R packages on first run.

---

## 🖥️ How to Run Locally

### 🔹 Simple Method (No Git Required)

If you're not familiar with Git or command line:

1. Go to:
   👉 [https://github.com/rarabzad/RDRS/tree/main/scripts/app](https://github.com/rarabzad/RDRS/tree/main/scripts/app)

2. Click the green **`Code`** button (top right), then select **`Download ZIP`**

3. After downloading:

   * Extract the ZIP file
   * Navigate into the folder:
     `RDRS-main/scripts/app`

4. Install [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) (if you haven’t already)

5. Open **RStudio**, click **File → Open File...**, and select the `app.R` file inside that folder

6. Click **Run App** (top-right of the RStudio window)

> 💡 The app will automatically install missing R packages the first time you run it. Please be patient — it may take a few minutes depending on your internet speed.

---

### 🔹 Git Method

```bash
git clone https://github.com/rarabzad/RDRS.git
cd RDRS/scripts/app
```

Then open `app.R` in RStudio or run:

```r
shiny::runApp()
```

---

## 📁 Required Inputs

### 1. NetCDF File (`.nc`)

* Must contain:

  * **2D latitude and longitude variables** (e.g., `lat[rlat, rlon]`, `lon[rlat, rlon]`)
  * **Time-dependent variable(s)** such as temperature or precipitation (e.g., `tas[time, rlat, rlon]`)
  * **Dimension names** for the spatial grid axes (e.g., `rlat`, `rlon`)

> 📝 You’ll be prompted to select:
>
> * Two variables representing **longitude and latitude**
> * Two **dimensions** representing the grid (e.g., `rlat`, `rlon`)

### 2. HRU Shapefile (`.zip`)

* Upload a zipped shapefile that includes:

  * `.shp`, `.shx`, `.dbf`, `.prj` (minimum required)
* Files must be **at the root of the ZIP** (not inside a folder)
* The shapefile should contain **polygon features** representing HRUs

> 📝 You’ll be prompted to select one **attribute field** (e.g., `SubId`) to use as the HRU ID

---

## 🎛️ Output Files (in Downloaded ZIP)

| File              | Description                             |
| ----------------- | --------------------------------------- |
| `grid_cells.shp`  | Shapefile of intersected NetCDF grid    |
| `hru_cells.shp`   | Processed HRU polygons used in analysis |
| `centroids.shp`   | Grid cell centroids                     |
| `weights.txt`     | Area-weighted table of HRU–grid pairs   |
| `grid_cells.json` | GeoJSON version of grid polygons        |
| `plot.pdf`        | Optional visualization plot (PDF)       |

---

## 🗺️ Interface Overview

| Component         | Description                                   |
| ----------------- | --------------------------------------------- |
| NetCDF Upload     | Upload `.nc` file with spatial variables      |
| Shapefile Upload  | Upload `.zip` with shapefile components       |
| Select Variables  | Choose `lon` / `lat` variables from NetCDF    |
| Select Dimensions | Choose grid dimensions (`rlat`, `rlon`, etc.) |
| Select HRU Field  | Choose column in shapefile for HRU ID         |
| Show Map          | Toggle interactive Leaflet map                |
| Download Results  | Download ZIP with all outputs                 |
| Log Output        | View real-time logs and error messages        |

---

## ⚠️ Notes

* Max upload size: **100 MB**
* Long processing time (2–5 mins) expected for large inputs
* Geometry validation is applied before intersection
* Ensure spatial alignment between shapefile and NetCDF grid (e.g., both in geographic coordinates like WGS84)

---

## 🙋 Need Help?

If you're having issues or want to suggest improvements:

* [Open an Issue](https://github.com/rarabzad/RDRS/issues)
