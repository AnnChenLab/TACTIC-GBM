FROM rocker/r-ver:4.5.3

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl cmake libuv1-dev \
    build-essential gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev \
    libicu-dev libcairo2-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# eulerr needs rustc >= 1.88.0; Ubuntu's apt rustc is too old, so install via rustup
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

RUN R -e "install.packages(c('shiny','plotly','ggplot2','data.table','dplyr','coin','openxlsx','tidyr','DT','htmltools','VennDiagram','readr','eulerr','bslib','survival','survminer','gridExtra','pheatmap','tibble','stringr'), repos='https://cloud.r-project.org')"

# OpenBLAS's runtime CPU-feature detection can pick an instruction set the host doesn't
# actually support, crashing with SIGILL (exit 132) inside BLAS/LAPACK-backed calls -
# hit in production (x86_64) via grDevices::convertColor, triggered by eulerr's
# Euler-diagram rendering. Pin to the reference BLAS/LAPACK implementation instead.
# (Multiarch triplet detected at build time so this also works for local arm64 testing.)
RUN ARCH=$(dpkg-architecture -qDEB_HOST_MULTIARCH) && \
    update-alternatives --set libblas.so.3-${ARCH} /usr/lib/${ARCH}/blas/libblas.so.3 && \
    update-alternatives --set liblapack.so.3-${ARCH} /usr/lib/${ARCH}/lapack/liblapack.so.3

WORKDIR /app

COPY app.R app.R
COPY data/ data/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app/app.R', host='0.0.0.0', port=3838)"]
