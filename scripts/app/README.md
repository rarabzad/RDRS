# Grid Weights Generator (Shiny App)

This Shiny web application allows users to **intersect NetCDF grid data with HRU (Hydrological Response Unit) polygons** to generate area-weighted grid-cell weights. It is especially useful for hydrological or land-surface model preprocessing, where spatial aggregation of gridded climate or forcing data is required.

---

## 📦 Features

- Upload NetCDF file (`.nc`) with 2D latitude/longitude grids.
- Upload HRU shapefile (`.zip`) with polygon features.
- Select variables and dimensions interactively from NetCDF file.
- Choose the attribute in the shapefile to use as HRU ID.
- Generate:
  - Intersected grid polygons
  - HRU boundaries
  - Grid centroids
  - Area-based weights for each HRU–grid pair
- Download all results in a single ZIP file.
- Optional: visualize the grid and HRU overlay on an interactive Leaflet map.
- Automatically checks for and installs missing R packages on first run.

---

## 🖥️ How to Run Locally

### 1. Prerequisites

Make sure you have:
- R installed (≥ 4.1 recommended)
- RStudio (optional but recommended)
- Internet connection (for downloading required packages and external script)

### 2. Clone or Download This Repository

```bash
git clone https://github.com/<your-username>/<your-repo-name>.git
cd <your-repo-name>
````

Or download the ZIP and extract it.

### 3. Launch the App

Open R or RStudio and set the working directory to the project folder.

Then run:

```r
shiny::runApp()
```

> 💡 On first run, the app will automatically install all required packages. This might take a few minutes depending on your internet speed.

---

## 📁 Required Inputs

### 1. NetCDF File (`.nc`)

* Must contain:

  * **2D latitude and longitude variables** (e.g., `lat[rlat, rlon]`, `lon[rlat, rlon]`)
  * **Time-dependent variable(s)** such as temperature or precipitation (e.g., `tas[time, rlat, rlon]`)
  * **Dimension names** for the spatial grid axes (e.g., `rlat`, `rlon`)

> 📝 You’ll be prompted to select:
>
> * Two variables representing **longitude and latitude** (e.g., `lon`, `lat`)
> * Two **dimensions** that represent the spatial grid (e.g., `rlat`, `rlon`)

### 2. HRU Shapefile (`.zip`)

* Upload a zipped shapefile that contains:

  * `.shp`, `.shx`, `.dbf`, and `.prj` files (minimum)
* Make sure all files are inside the ZIP root (not in a folder)
* The shapefile must contain polygon features, each representing an HRU

> 📝 You’ll be prompted to select:
>
> * One column from the shapefile to serve as **HRU ID** (e.g., `SubId`)

---

## 🎛️ Output Files (downloadable as ZIP)

* `grid_cells.shp` — Shapefile of NetCDF grid polygons
* `hru_cells.shp` — Original HRU shapefile (used in processing)
* `centroids.shp` — Grid cell centroids
* `weights.txt` — Text file with HRU-to-grid cell weights
* `grid_cells.json` — GeoJSON of grid polygons
* `plot.pdf` — Optional plot saved as PDF (if plotting enabled)

---

## 🗺️ Interface Overview

| Component         | Description                                        |
| ----------------- | -------------------------------------------------- |
| NetCDF Upload     | Upload `.nc` file with spatial variables           |
| Shapefile Upload  | Upload `.zip` with shapefile components            |
| Select Variables  | Choose lon/lat variables from NetCDF               |
| Select Dimensions | Choose grid dimension names (e.g., `rlat`, `rlon`) |
| Select HRU Field  | Choose column in shapefile as HRU ID               |
| Show Map          | Toggle interactive Leaflet map                     |
| Download Results  | Get ZIP archive with all output layers             |
| Log Output        | View processing messages and errors                |

---

## ⚠️ Notes

* App supports up to 100 MB file upload.
* Long processing time is expected (2–5 mins) for large shapefiles or grids.
* Geometry validation is internally applied for intersections.
* For best results, ensure the shapefile and NetCDF grid are in **compatible coordinate systems** (e.g., both in geographic lat-lon or reprojected appropriately).

---

## 📄 License

MIT License. See `LICENSE.md` for details.

---

## 🙋 Need Help?

If you're having trouble running the app locally or have a dataset question, feel free to open an issue or start a discussion on the GitHub repository.

Let me know if you'd like:
- Me to prepare this as an actual file you can copy-paste.
- A sample NetCDF file and shapefile to use as demo inputs.
- Instructions for deploying it with Docker or any platform.
