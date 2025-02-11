FROM ubuntu:22.04
MAINTAINER Chandana Tennakoon <chandana.tennakoon@pirbright.ac.uk>
ENV DEBIAN_FRONTEND=noninteractive 

# Update package list and install Python 3.10
RUN apt-get update && apt-get install -y \
    wget gpg \
    python3.10 \
    python3-pip && \
    ln -sf /usr/bin/python3.10 /usr/bin/python3
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/cran.list && \
    apt-get update && \
    apt-get install -y r-base r-base-dev && \

    R -e "install.packages('remotes', repos='http://cran.us.r-project.org')" && \
	R -e "remotes::install_version('adegenet', version='2.1.10', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('ape', version='5.8', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('ggExtra', version='0.10.1', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('ggnewscale', version='0.4.10', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('ggplot2', version='3.4.0', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('hierfstat', version='0.5-11', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('irlba', version='2.3.5.1', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('jsonlite', version='1.8.4', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('microseq', version='2.1.6', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('phangorn', version='2.11.1', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('phytools', version='2.3-0', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('Rcpp', version='1.0.9', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('RcppArmadillo', version='0.11.4.3.1', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('rentrez', version='1.2.3', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('restez', version='2.1.4', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('seqinr', version='4.2-36', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('stringr', version='1.5.0', repos='http://cran.us.r-project.org')"  && \
	R -e "remotes::install_version('tibble', version='3.1.8', repos='http://cran.us.r-project.org')"  

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('treeio')"
#	R -e "remotes::install_version('grid', version='4.2.2', repos='http://cran.us.r-project.org')"  && \

RUN apt-get install -y curl build-essential
WORKDIR /app

# Install Rust and Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
# Add Rust to PATH
ENV PATH="/root/.cargo/bin:${PATH}"
# Verify installation
RUN rustc --version && cargo --version
RUN cargo install gfautil

COPY requirements.txt requirements.txt 
#COPY script.py .

# Install dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
COPY bin/tree.R /app/
COPY bin/plot_variation_map.R /app/
COPY bin/plot_SNP_density.R /app/
COPY bin/visualisation.R /app/
COPY bin/plot_pavs.R /app/
COPY bin/annotate.R /app/
COPY bin/view_gml.R /app/
COPY bin/panaroo_viz.R /app/
