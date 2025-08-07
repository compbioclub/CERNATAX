# CERNATAX

CERNATAX detects ceRNA network mRNA-miRNA-lncRNA triplet axis from RNA expression profiles.

We have manually curated a ceRNA network database by integrating TargetSCAN 8.0, miRTarBase 9.0, miRDB 6.0, NPInter 4.0, ENCORI/starBase 2.0, miRWalk V3, and RNAInter in 2020. 

To maximize both the accuracy and biological relevance of ceRNA regulatory inference in CERNATAX, we employed a comprehensive and multi-layered strategy for integrating RNA interaction data from diverse sources and disease contexts. All ceRNA interactions were curated exclusively from well-established, high-confidence databases with robust experimental or computational support, and then standardized to a unified gene annotation system to ensure compatibility and reduce nomenclature discrepancies. While the initial background network is constructed broadly to enhance coverage and reliability, we rigorously filter this network for each disease dataset, retaining only those interactions that involve differentially expressed genes identified in the specific disease context. This two-step process ensures that the final ceRNA network is highly disease- and dataset-specific, thereby minimizing the inclusion of irrelevant or non-specific interactions. Additionally, to further safeguard data integrity and comparability, we strongly recommend users normalize RNA expression data and remove batch effects prior to analysis. Collectively, these measures ensure that the ceRNA networks constructed by CERNATAX are both accurate and tailored to the biological context under investigation.


To run CERNATAX with full reference ceRNA network, please download the full db file [`ceRNA_database.csv`](https://doi.org/10.5281/zenodo.15357964), and place it into `cernatax/db'.

The tutorial and case study of the SCZ cohort can be found at http://compbioclub.github.io/CERNATAX/.
