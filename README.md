# sliceanddice
Under development.

This is the repository for a tool and project currently called ```sliceanddice```, a tool with the goal of working with large genomic data and especially single-cell RNA-seq data in mind. The tool aims to simplify download and manipulation of sc-RNA-seq data from the Gene Expression Omnibus and perform light-weight gene coexpression analysis. Given a dataset and list of genes of interest, the tool will provide information about correlated genes across samples of the gene expression matrix. This tool is being written entirely in python and currently exists on jupyter-notebooks as well as a .py file which will ultimately become the basis for API and CLI usage.

### Installation

fork and clone into this repository to utilize ```sliceanddice```.


### Usage Examples

Example usage can be found in the notebooks directory. There are two notebooks. ```alex-MVN-simulation.ipynb``` includes an implementation of pymc3 on a simulated dataset. ```alex_final_notebook.ipynb``` shows the application of the code to a single-cell RNA-seq dataset derived from patients with Head and Neack Squamous Cell Carcinoma (HNSCC). This notebook contains the current state of the code and will be undergoing active development.

