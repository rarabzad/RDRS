FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    libv8-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    libjq-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean

# Install R packages
RUN R -e "install.packages(c('shiny', 'leaflet', 'sf', 'DT', 'shinyWidgets', 'zip', 'shinyjs', 'ncdf4', 'geosphere', 'dplyr', 'sp', 'lwgeom', 'rmapshaper'), repos='https://cloud.r-project.org')"

# Copy the app
COPY . /srv/shiny-server/

# Make sure permissions are correct
RUN chown -R shiny:shiny /srv/shiny-server

# Expose default Shiny port
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
