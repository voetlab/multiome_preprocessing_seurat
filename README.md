# Multiome_preprocessing_seurat

Workflow to preprocess single nucleus multiome data (gene expression + ATAC) produced with 10X Genomics Cellranger ARC in Seurat. Using a sequential set of scripts, go from a folder with cellranger-arc outputs to a Seurat object with dimensionality reduction and clustering for gene expression data.


# Getting started

## Environment

A conda environment with both [Seurat](https://satijalab.org/seurat/) and [Tidyverse](https://www.tidyverse.org/) installed is recommended. Preprocessing uses dedicated packages [SoupX](https://github.com/constantAmateur/SoupX) and [scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html), while integration uses [harmony](https://github.com/immunogenomics/harmony).

## Data organization

The script assumes cellranger arc output has been reorganized to have all filtered_feature_bc_matrix.h5 files for all samples into one folder with the names changed to sample ids. Alongside this, raw_feature_bc_matrix.h5 files are also placed in a single folder with names changed to sample ids. The workflow expects pairs of matrices to be available and have the same name for every sample.

## Running the workflow

The workflow consists of 5 scripts that sequentially perform the following:

    1. Run SoupX and doublet finder
    2. Ensure cells are renamed to include sample id
    3. Apply quality control thresholds and exclude low-quality samples
    4. Perform dimensionality reduction after count normalization
    5. Integrate samples across the most relevant biological unit and find clusters

The output is a ready-to-annotate seurat object.

06_commands.sh is the configuration script and is the only one that should be adapted. Key adaptations would be the following:

* out-folder: where should outputs be saved/collected from?
* matrix/raw_matrices: the folder with the collected filtered/raw h5 files.
* QC: a comma separated vector of quality control metrices in the following order
    * min number of features
    * max number of features
    * min number of read, max number of reads
    * max mitochondrial percentage
    * max ribosomal percentage
* exclude: sample ids to exclude.

# Contact

Made by [Markus Boesch](markus.boesch@kuleuven.be), [Christos Vlachos](), [Emiel Geeraerts](emiel.geeraerts@kuleuven.be)

# Citation
10.5281/zenodo.13620951
