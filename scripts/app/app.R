required_packages <- c(
  "shiny", "leaflet", "sf", "DT", "shinyWidgets", "zip", "shinyjs",
  "ncdf4", "geosphere", "dplyr", "sp", "lwgeom", "rmapshaper"
)

install_log <- ""

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_log <- paste0(install_log, Sys.time(), " — Installing missing package: ", pkg, "\n")
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org")
      install_log <- paste0(install_log, Sys.time(), " — ✅ Successfully installed: ", pkg, "\n")
    }, error = function(e) {
      install_log <- paste0(install_log, Sys.time(), " — ❌ Failed to install ", pkg, ": ", e$message, "\n")
    })
  }
}

if (install_log == "") {
  install_log <- paste0(Sys.time(), " — ✅ All required packages are already installed.\n")
}

library(shiny)
library(leaflet)
library(sf)
library(DT)
library(shinyWidgets)
library(zip)
library(shinyjs)
library(ncdf4)
library(geosphere)
library(dplyr)
library(sp)
library(lwgeom)
library(rmapshaper)

setwd(if (!is.null(sys.frame(1)$ofile)) dirname(sys.frame(1)$ofile) else ".")
source("https://raw.githubusercontent.com/rarabzad/RDRS/refs/heads/main/scripts/app/grids_weights_generator.R")
options(shiny.maxRequestSize = 50 * 1024^2)  # 100 MB
initial_log <- install_log

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Grid Weights Generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("ncfile", "Upload NetCDF File (.nc)", accept = ".nc"),
      fileInput("shpfile", "Upload HRU Shapefile (.zip)", accept = ".zip"),
      uiOutput("var_select"),
      uiOutput("dim_select"),
      uiOutput("hru_select"),
      checkboxInput("show_map", "Show Map", value = TRUE),
      actionButton("generate", "Generate Weights", icon = icon("play")),
      div(id = "waiting_msg", style = "color: blue; font-style: italic; margin-top: 10px;"),
      br(),
      downloadButton("download_zip", "Download Results")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 12,
          div(
            style = "height: 250px; overflow-y: auto; border: 1px solid #ccc; padding: 10px; margin-bottom: 10px;",
            strong("Log Output"),
            verbatimTextOutput("log")
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          conditionalPanel(
            condition = "input.show_map == true",
            leafletOutput("map", height = "500px")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  log_data <- reactiveVal(initial_log)
  shp_path_r <- reactiveVal(NULL)  # Store path of extracted .shp file
  zip_path <- reactiveVal(NULL)
  result_data <- reactiveVal(NULL)
  nc_data <- reactiveVal(NULL)
  shp_path_r <- reactiveVal(NULL)
  
  log_text <- function(msg) {
    if (is.null(msg) || msg == "") return()  # skip empty messages
    isolate({
      old <- log_data()
      new <- paste0(old, Sys.time(), " — ", msg, "\n")
      log_data(new)
    })
  }
  
  update_log <- function(msg) {
    isolate({
      current_log <- log_data()
      new_log <- paste0(current_log, Sys.time(), " — ", msg, "\n")
      log_data(new_log)
    })
  }
  
  nc_data <- reactiveVal(NULL)
  shp_data <- reactiveVal(NULL)
  result_data <- reactiveVal(NULL)
  
  observeEvent(input$ncfile, {
    req(input$ncfile)
    nc <- ncdf4::nc_open(input$ncfile$datapath)
    nc_data(nc)
    
    vars <- names(nc$var)
    dims <- names(nc$dim)
    
    output$var_select <- renderUI({
      req(nc_data())
      vars <- names(nc_data()$var)
      
      selected_vars <- input$varnames
      # Make sure selected_vars are in vars and unique
      selected_vars <- selected_vars[selected_vars %in% vars]
      # Put selected_vars first, then rest
      ordered_vars <- c(selected_vars, setdiff(vars, selected_vars))
      
      selectizeInput("varnames", "Select Variable Names (lon, lat)", 
                     choices = ordered_vars,
                     selected = selected_vars,
                     multiple = TRUE,
                     options = list(maxItems = 2))
    })
    
    output$dim_select <- renderUI({
      req(nc_data())
      dims <- names(nc_data()$dim)
      
      selected_dims <- input$dimnames
      selected_dims <- selected_dims[selected_dims %in% dims]
      ordered_dims <- c(selected_dims, setdiff(dims, selected_dims))
      
      selectizeInput("dimnames", "Select Dimension Names (rlat, rlon)", 
                     choices = ordered_dims,
                     selected = selected_dims,
                     multiple = TRUE,
                     options = list(maxItems = 2))
    })
    update_log("✅ NetCDF file loaded.")
  })
  
  observeEvent(input$shpfile, {
    req(input$shpfile)
    log_text("📦 Extracting uploaded shapefile ZIP...")
    
    unzip_dir <- tempfile(pattern = "shpzip_")
    dir.create(unzip_dir)
    unzip(input$shpfile$datapath, exdir = unzip_dir)
    
    shp_files <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) {
      log_text("❌ Error: No .shp file found in the ZIP archive.")
      return()
    }
    
    shp_path_r(shp_files[1])  # Store the shapefile path here
    
    shp_data_tmp <- tryCatch(sf::st_read(shp_path_r(), quiet = TRUE), error = function(e) NULL)
    if (is.null(shp_data_tmp)) {
      log_text("❌ Error: Failed to read shapefile.")
      return()
    }
    
    shp_data(shp_data_tmp)  # assign to reactiveVal
    
    
    # populate HRU_ID select input dynamically
    output$hru_select <- renderUI({
      req(shp_data())
      selectInput("hru_id", "Select HRU ID Field", choices = names(shp_data()), selected = names(shp_data())[1])
    })
    if (!is.null(shp_data())) {
      log_text(sprintf("✅ Shapefile loaded with %d features.", nrow(shp_data())))
    } else {
      log_text("⚠️ Shapefile loaded, but no features found.")
    }
  })
  
  observeEvent(input$generate, {
    req(input$ncfile, shp_path_r(), input$hru_id)
    shinyjs::html("waiting_msg", "⏳ Please wait ...")
    
    outdir <- tempfile()
    dir.create(outdir)
    owd <- setwd(outdir)
    on.exit(setwd(owd), add = TRUE)
    
    tryCatch({
      withProgress(message = "Processing grid weights...", value = 0, {
        incProgress(0.1, detail = "Running grids_weights_generator()...")
        res <- grids_weights_generator(
          ncfile   = input$ncfile$datapath,
          hrufile  = shp_path_r(),
          varnames = input$varnames,
          dimnames = input$dimnames,
          HRU_ID   = input$hru_id,
          plot     = TRUE
        )
        incProgress(0.4, detail = "Generating leaflet map and writing shapefiles...")
        result_data(res)
        
        output$map <- renderLeaflet({
          req(input$show_map)
          req(result_data())
          leaflet() %>%
            addTiles() %>%
            addPolygons(data = res$grid_sf, color = "black", weight = 2, fillOpacity = 0.3) %>%
            addPolygons(data = res$hru_sf, color = "red", weight = 0.5, fillOpacity = 0.3) %>%
            addCircleMarkers(data = res$centroids, color = "blue", radius = 2)
        })
        
        # Write outputs
        grid_cells_shp_path <- file.path(outdir, "grid_cells.shp")
        grid_cells_json_path <- file.path(outdir, "grid_cells.json")
        hru_cells_shp_path <- file.path(outdir, "hru_cells.shp")
        centroids_shp_path <- file.path(outdir, "centroids.shp")
        weights_txt_path <- file.path(outdir, "weights.txt")
        
        sf::st_write(res$grid_sf, grid_cells_shp_path, delete_layer = TRUE, quiet = TRUE)
        log_text(paste("✅ Shapefile written at:", grid_cells_shp_path))
        incProgress(0.1)
        
        sf::st_write(res$grid_sf, grid_cells_json_path, driver="GeoJSON", delete_layer = TRUE, quiet = TRUE)
        log_text(paste("✅ GeoJSON written at:", grid_cells_json_path))
        incProgress(0.05)
        
        sf::st_write(res$hru_sf,  hru_cells_shp_path,  delete_layer = TRUE, quiet = TRUE)
        sf::st_write(res$centroids, centroids_shp_path, delete_layer = TRUE, quiet = TRUE)
        log_text(paste("✅ Centroids shapefile written at:", centroids_shp_path))
        incProgress(0.1)
        
        writeLines(res$weights_txt, weights_txt_path)
        log_text(paste("✅ Weights text file written at:", weights_txt_path))
        incProgress(0.05)
        
        # Generate plot PDF if plot_path function is available
        if (!is.null(res$plot_path) && is.function(res$plot_path)) {
          pdf_path <- file.path(outdir, "plot.pdf")
          res$plot_path(pdf_path)  # Call function to save plot to PDF file
          log_text(paste("✅ Plot PDF created at:", pdf_path))
          incProgress(0.05)
        }
        
        # Zip shapefile components and other outputs
        shapefile_bases <- c("grid_cells", "hru_cells", "centroids")
        shapefile_extensions <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")
        
        shapefile_files <- unlist(lapply(shapefile_bases, function(base) {
          file.path(outdir, paste0(base, shapefile_extensions))
        }))
        shapefile_files <- shapefile_files[file.exists(shapefile_files)]
        
        extra_files <- c(file.path(outdir, "weights.txt"), file.path(outdir, "plot.pdf"), file.path(outdir, "grid_cells.json"))
        extra_files <- extra_files[file.exists(extra_files)]
        
        files_to_zip <- c(shapefile_files, extra_files)
        
        # Normalize and relative paths for zip
        outdir_norm <- normalizePath(outdir, winslash = "/")
        files_norm <- normalizePath(files_to_zip, winslash = "/")
        files_rel <- sub(paste0("^", outdir_norm, "/"), "", files_norm)
        
        zipfile <- file.path(tempdir(), "weights_output.zip")
        
        incProgress(0.1, detail = "Creating ZIP archive...")
        zip::zip(zipfile, files = files_rel)
        incProgress(0.1)
        
        zip_path(zipfile)
        
        log_text("✅ Weight calculation complete and ZIP archive created.")
      })
    }, error = function(e) {
      log_text(paste0("❌ Error during processing: ", e$message))
    })
    shinyjs::html("waiting_msg", "")  # clear the message
  })
  
  
  output$log <- renderText({
    log_data()
  })
  
  output$download_zip <- downloadHandler(
    filename = function() {
      "weights_output.zip"
    },
    content = function(file) {
      req(zip_path())  # ensure the zip path exists
      file.copy(zip_path(), file)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui = ui, server = server)
