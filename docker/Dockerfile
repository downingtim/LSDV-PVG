FROM ubuntu:22.04
MAINTAINER Chandana Tennakoon <chandana.tennakoon@pirbright.ac.uk>
ENV DEBIAN_FRONTEND=noninteractive 

# Update package list and install Python 3.10
RUN apt-get update && apt-get install -y \
    wget gpg  libcurl4-openssl-dev libxml2-dev \
    python3.10 \
    python3-pip && \
    ln -sf /usr/bin/python3.10 /usr/bin/python3 && \
    ln -sf /usr/bin/python3.10 /usr/bin/python
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/cran.list && \
    apt-get update && \
    apt-get install libssl-dev && \
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
RUN R -e "BiocManager::install(update=TRUE,ask=FALSE)"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('treeio')"

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
COPY bin.docker/tree.R /bin/
COPY bin.docker/plot_variation_map.R /bin/
COPY bin.docker/plot_SNP_density.R /bin/
COPY bin.docker/visualisation.R /bin/
COPY bin.docker/annotate.R /bin/
COPY bin.docker/view_gml.R /bin/
COPY bin.docker/panaroo_viz.R /bin/

RUN cargo install fastix
RUN apt-get update && apt-get install git -y && git clone https://github.com/gi-bielefeld/pangrowth && cd pangrowth && make
RUN mv pangrowth/pangrowth /usr/local/bin/
COPY bin.docker/raxml-ng /usr/local/bin

RUN git clone https://github.com/vgl-hub/gfastats.git --recursive
RUN cd gfastats && make -j
RUN cd /app && cp gfastats/build/bin/gfastats /usr/local/bin
COPY bin.docker/lengths.py /bin
COPY bin.docker/plot_single_hist.py /bin
COPY bin.docker/plot_growth.py /bin
COPY bin.docker/plot_core.py /bin
COPY bin.docker/afs.py /bin
COPY bin.docker/parseGB.py /bin
COPY bin.docker/plot_hist.py /bin
RUN apt-get install ncbi-entrez-direct -y
RUN apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra -y
RUN wget https://github.com/vgteam/vg/releases/download/v1.64.1/vg && chmod +x vg && mv vg /bin/vg
RUN apt-get install graphviz -y
RUN apt-get install mafft -y
RUN apt-get install samtools -y
RUN apt-get install -y gffread
COPY bin.docker/plot_pavs.R /bin/
COPY bin.docker/plot_pavs.R /usr/bin/






