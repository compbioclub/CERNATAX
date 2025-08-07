# CERNATAX

Welcome to the official documentation of **CERNATAX**!

## Overview

**CERNATAX** detects ceRNA network `mRNA-miRNA-lncRNA triplet axis` from RNA expression profiles.

We have manually curated a ceRNA network database by integrating `TargetSCAN 8.0`, `miRTarBase 9.0`, `miRDB 6.0`, `NPInter 4.0`, `ENCORI/starBase 2.0`, `miRWalk V3`, and `RNAInter in 2020`. 

To maximize both the accuracy and biological relevance of ceRNA regulatory inference in CERNATAX, we employed a comprehensive and multi-layered strategy for integrating RNA interaction data from diverse sources and disease contexts. All ceRNA interactions were curated exclusively from well-established, high-confidence databases with robust experimental or computational support, and then standardized to a unified gene annotation system to ensure compatibility and reduce nomenclature discrepancies. While the initial background network is constructed broadly to enhance coverage and reliability, we rigorously filter this network for each disease dataset, retaining only those interactions that involve differentially expressed genes identified in the specific disease context. This two-step process ensures that the final ceRNA network is highly disease- and dataset-specific, thereby minimizing the inclusion of irrelevant or non-specific interactions. Additionally, to further safeguard data integrity and comparability, we strongly recommend users normalize RNA expression data and remove batch effects prior to analysis. Collectively, these measures ensure that the ceRNA networks constructed by CERNATAX are both accurate and tailored to the biological context under investigation.

To run CERNATAX with full reference ceRNA network, please download the full db file [`ceRNA_database.csv`](https://doi.org/10.5281/zenodo.15357964), and place it into `cernatax/db`.


## Getting Started

Want to start using it immediately? Check out the [Installation Guide](installation.md).


## Tutorial Guide
The followings are tutorials of how to use CERNATAX on the SCZ cohort:

-   [Basic operations and stats for the reference ceRNA network](tutorial/reference_ceRNA_network.ipynb)
-   [Use DEG to get cohort-specific and disease relatated ceRNA axis](tutorial/ceRNA_axis_from_DEG.ipynb)
-   [ceRNA expression visualization for a cohort](tutorial/plot_ceRNA_exp_for_cohort.ipynb)
-   [ceRNA-axis correlation analysis for a cohort](tutorial/ceRNA_correlation_analysis.ipynb)
-   [Predict disease phenotype using ceRNA axis](tutorial/prediction.ipynb)

## Citation

If you use **CERNATAX** in your research, please cite the following paper:

APA format:

```
Liu, X., Jiang, A., Lyu, C., & Chen, L. (2025). Knowledge-driven annotation for gene interaction enrichment analysis. bioRxiv, 2025-04. https://doi.org/10.1101/2025.04.15.649030
```

BibTeX format:

```bibtex
@article{grea,
  title={Knowledge-driven annotation for gene interaction enrichment analysis},
  author={Liu, Xiaoyu and Jiang, Anna and Lyu, Chengshang and Chen, Lingxi},
  journal={bioRxiv},
  pages={2025--04},
  year={2025},
  doi={10.1101/2025.04.15.649030},
  publisher={Cold Spring Harbor Laboratory}
}
```


<div style="display:none;">
<script type='text/javascript' id='clustrmaps' src='//cdn.clustrmaps.com/map_v2.js?cl=ffffff&w=a&t=n&d=s_zp3a_kJX2eHlUNurnH4Jti8lf7sMFpJyhQQnn21MU'></script>
<script defer src='https://static.cloudflareinsights.com/beacon.min.js' data-cf-beacon='{"token": "aec36b862a47431a979dc263a1f98d74"}'></script>
</div>
