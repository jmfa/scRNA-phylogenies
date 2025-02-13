## **Estimating phylogenetic trees from (smart-seq2) scRNA-seq data.**

In this repository, we provide our recommended pipeline for performing single-cell phylogenetic analysis using scRNA-seq data. Our workflow includes:

- **Step 1:** [FASTQ Mapping and BAM Processing](https://github.com/jmfa/scRNA-phylogenies/tree/main/1.Mapping%2BProcessing)
- **Step 2&3:** [Variant calling + Genotyping](https://github.com/jmfa/scRNA-phylogenies/tree/main/2.VariantCalling%2BGenotyping)
- **Step 4:** [Phylogenetic Tree Reconstruction (with bootstrap support estimation)](https://github.com/jmfa/scRNA-phylogenies/tree/main/3.TreeReconstruction)
- **Step 5:** [Read Count Normalization](https://github.com/jmfa/scRNA-phylogenies/tree/main/4.ReadCountsNormalization)
- **Step 6:** [Phylogenetic Signal Estimation](https://github.com/jmfa/scRNA-phylogenies/tree/main/5.PhylogeneticSignal)

For further details on the methodology and its application, please refer to our publication:
**"Unraveling the phylogenetic signal of gene expression from single-cell RNA-seq data"**
https://www.biorxiv.org/content/10.1101/2024.04.17.589871v1
