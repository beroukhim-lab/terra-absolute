FROM rocker/tidyverse:4.3.3

# Install missing CRAN packages (tidyverse already present)
RUN R -e 'install.packages(c("data.table","optparse"), repos="https://cloud.r-project.org")'

WORKDIR /opt/app
COPY capseg_conv.R /opt/app/capseg_conv.R

# Default entrypoint makes WDL command shorter (optional but convenient)
ENTRYPOINT ["Rscript", "/opt/app/capseg_conv.R"]