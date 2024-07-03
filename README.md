This repo contains all computational analyses, in the form of [Jupyter Notebooks](https://jupyter.org/), used in the paper: *Resolving the source of branch length variation in the Y chromosome phylogeny*.

**Running the Notebooks**:
- Before running the notebooks, the input VCF data should first be downloaded (Zenodo link) and extracted to `data/input_data/`.

- The notebooks are intended to be run in numerical order. The full computational analysis should be completed in less than a day.

- The final size of the directory is ~100GB.

- Data generated from the analyses, which is used for plotting, will be stored in `data/output_data/`.

**Additional Resources**:
- The `human_chimp_divergence` directory includes Python scripts adapted from [Great Ape Y Evolution](https://github.com/makovalab-psu/great-ape-Y-evolution)

- The `BEAST` directory contains necessary .xml files (generated by Beauti) for BEAST analyses. It also contains .log and .trees files generated by BEAST. The .treeannotator files within are used for plotting summary trees as shown in the paper.

- The [y_chr_utils](https://github.com/Yaniv42/y_chr_utils) repository is included as a submodule.

## Installation
The required python packages are `pandas`, `numpy`, `statsmodels`, `matplotlib` and `pydantic`

The required R packages are `phangorn`, `lemon`, `ggplot2`, `ggtree` and `treeio`