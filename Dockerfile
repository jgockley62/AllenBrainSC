FROM rocker/tidyverse:4.0.0

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libglpk40\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y git\
&& R -e "install.packages('BiocManager')"\
&& R -e "BiocManager::install('biomaRt')"\
&& R -e "BiocManager::install('ComplexHeatmap')"\
&& R -e "install.packages('data.table')"\
&& R -e "install.packages('doParallel')"\
&& R -e "BiocManager::install('edgeR')"\
&& R -e "install.packages('foreach')"\
&& R -e "BiocManager::install('GenomicRanges')"\
&& R -e "devtools::install_github('brian-bot/githubr')"\
&& R -e "install.packages('ggplot2')"\
&& R -e "BiocManager::install('IRanges')"\
&& R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
&& R -e "devtools::install_github('Sage-Bionetworks/knit2synapse')"\
&& R -e "install.packages('knitr')"\
&& R -e "BiocManager::install('limma')"\
&& R -e "install.packages('mclust')"\
&& R -e "install.packages('psych')"\
&& R -e "install.packages('RColorBrewer')"\
&& R -e "install.packages('rlang')"\
&& R -e "install.packages('R.utils')"\
&& R -e "install.packages('statmod')"\
&& R -e "install.packages('stringr')"\
&& R -e "BiocManager::install('sva')"\
&& R -e "install.packages('utils')"

run mkdir /home/nperumal/<AllenBrainSC \ 
&& git clone https://github.com/NitheshPerumal/AllenBrainSC.git /home/nperumal/AllenBrainSC/