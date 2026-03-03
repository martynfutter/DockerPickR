# =============================================================
#  ParetoPick-R — Dockerfile
#  Base image: rocker/shiny (R 4.4.2 — closest stable to 4.5.2)
#  NOTE: Change to rocker/shiny:4.5.2 once it is published on
#        Docker Hub (usually lags a few weeks behind R releases).
# =============================================================
FROM rocker/shiny:4.4.2

# -----------------------------------------------------------------
# 1. System libraries required by R packages
#    sf / spdep / tmap  -> GDAL, GEOS, PROJ, units
#    svglite            -> fontconfig, harfbuzz, fribidi
#    curl / openssl     -> HTTPS / network support
#    chromium           -> replaces PhantomJS for webshot2 map exports
# -----------------------------------------------------------------


RUN apt-get update && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libglpk-dev \       
    chromium \
    chromium-driver \
    git \
    && rm -rf /var/lib/apt/lists/*

# Tell webshot2 where to find the browser
ENV CHROMOTE_CHROME=/usr/bin/chromium

# -----------------------------------------------------------------
# 2. Restore R packages from renv.lock
#    Copying only the lockfile first lets Docker cache this expensive
#    layer separately from app-code changes.
# -----------------------------------------------------------------
WORKDIR /srv/shiny-server

COPY renv.lock renv.lock

RUN R -e "install.packages('renv', repos='https://cloud.r-project.org'); \
          install.packages('remotes', repos='https://cloud.r-project.org'); \
          renv::restore(lockfile='renv.lock', prompt=FALSE)"

# -----------------------------------------------------------------
# 3. Install webshot2 to replace webshot + PhantomJS
#    *** You must also update global.R (see README note below) ***
# -----------------------------------------------------------------
RUN R -e "install.packages('webshot2', repos='https://cloud.r-project.org')"

# -----------------------------------------------------------------
# 4. Copy app source and Shiny Server config
# -----------------------------------------------------------------
COPY app/               /srv/shiny-server/app/
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

# Create the directories global.R writes to at runtime
RUN mkdir -p \
    /srv/shiny-server/data \
    /srv/shiny-server/input \
    /srv/shiny-server/output

# -----------------------------------------------------------------
# 5. Permissions — shiny user needs write access to data dirs
# -----------------------------------------------------------------
RUN chown -R shiny:shiny \
    /srv/shiny-server/data \
    /srv/shiny-server/input \
    /srv/shiny-server/output \
    /srv/shiny-server/app

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
