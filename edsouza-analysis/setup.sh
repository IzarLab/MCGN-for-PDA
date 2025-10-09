# This script isn't meant to be run
# If an environment ever breaks, just uncomment the conda create blocks
# and reproduce the envs from first principles
# Also contains a few other package installations that can't be done in conda

cd /home/ubuntu/edsouza-summer2023/cxcr4-pdac

################################################################################
# conda create -n cellbender
mamba env create -f cellbender.env.yml  # Optimized for aws p2.xlarge w/CUDA
conda activate cellbender
git clone https://github.com/broadinstitute/CellBender
pip install -e CellBender
rm -rf ./CellBender
conda deactivate


################################################################################
# conda create -n numbat
# mamba install -n numbat -c bioconda -c conda-forge \
#     cellsnp-lite samtools  r-essentials r-base 'r-seurat<5.0.0' \
#     r-irkernel Jupyter r-devtools r-circlize\
#     bioconductor-ggtree bioconductor-genomicranges  r-gifski r-ggpubr \
#     bioconductor-milor r-beeswarm bioconductor-destiny bioconductor-mast \
#     bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-ucell \
#     bioconductor-deseq2 r-matrix.utils bioconductor-fgsea
# conda activate numbat && conda env export > numbat.env.yml && conda deactivate
mamba env create -f numbat.env.yml
sudo apt update && sudo apt install unzip moreutils libgomp1 libcairo2-dev # Important to run before installing numbat!

# Fix errors in irlba
R -e 'install.packages("Matrix", type = "source")'
R -e 'install.packages("irlba", type = "source")'

sudo apt-get install libfftw3-dev  # Apparently metap requires this?
R -e 'install_github("wjawaid/enrichR")'


NUMBAT="/home/ubuntu/InstallTemp/numbat"  # Local dir for numbat requirements
mkdir -p "${NUMBAT}"
wget -O "${NUMBAT}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz" \
    https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
wget -O "${NUMBAT}/1000G_hg38.zip" \
    http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
wget -O "${NUMBAT}/Eagle_v2.4.1.tar.gz" \
    https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xvzf "${NUMBAT}/Eagle_v2.4.1.tar.gz"
gunzip "${NUMBAT}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
unzip "${NUMBAT}/1000G_hg38.zip"
rm "${NUMBAT}/Eagle_v2.4.1.tar.gz"; rm "${NUMBAT}/1000G_hg38.zip"
ln -s "${NUMBAT}/Eagle_v2.4.1/eagle" "${CONDA_PREFIX}/bin/eagle"  # IMPORTANT!


conda activate numbat
Rscript -e 'install.packages("numbat", dependencies=TRUE, repos="http://cran.us.r-project.org")'
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("immunogenomics/presto")'  # Makes Seurat::FindAllMarkers() faster
Rscript -e 'devtools::install_github("carmonalab/STACAS")'  # Required for projecTILs
Rscript -e 'devtools::install_github("carmonalab/ProjecTILs")'
Rscript -e 'devtools::install_github("bnprks/BPCells/r")'
# Rscript "${CONDA_PREFIX}/lib/R/library/numbat/bin/pileup_and_phase.R" to run preprocessing script


conda deactivate

################################################################################
# conda create -n infercnv
# conda activate infercnv
# mamba install -n infercnv -c bioconda -c conda-forge r-essentials r-base r-seurat bioconductor-infercnv r-viridis
# conda activate infercnv && conda env export > infercnv.env.yml && conda deactivate
mamba env create -f infercnv.env.yml
INFERCNV="/home/ubuntu/InstallTemp/infercnv"
mkdir -p "${INFERCNV}"
aws s3 cp s3://cxcr4-pdac/infercnv/input/refdata-gex-GRCh38-2020-A_gen_pos.txt "${INFERCNV}/refdata-gex-GRCh38-2020-A_gen_pos.txt"

################################################################################
# conda create -n cxcr4
# mamba install -n cxcr4 -c anaconda jupyter beautifulsoup4 pytables
# mamba install -n cxcr4 -c conda-forge -c bioconda scanpy python-igraph leidenalg scrublet r-seurat rpy2 r-hdf5r r-essentials r-base r-ggrastr r-pheatmap bioconductor-dropletutils bioconductor-singler bioconductor-celldex bioconductor-singlecellexperiment bioconductor-scater r-devtools
# conda activate cxcr4 && conda env export > cxcr4.env.yml && conda deactivate
mamba env create -f cxcr4.env.yml
conda activate cxcr4
Rscript -e 'devtools::install_github("navinlabcode/copykat")'
conda deactivate
sudo apt-get install libblas-dev liblapack-dev


##############################################################################
# conda create -n integration
# mamba install -c conda-forge -c bioconda r-essentials r-seurat -n integration
# conda activate integration && conda env export > integration.env.yml && conda deactivate