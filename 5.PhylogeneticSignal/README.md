## Step 6: Phylogenetic signal estimation

Now that the gene expression data is processed and transformed, we can proceed to estimate the phylogenetic signal. In our study, we use the **sensiPhy**[https://www.rdocumentation.org/packages/sensiPhy/versions/0.8.5] R package, which is designed to account for phylogenetic uncertainty when estimating the signal. Specifically, we will estimate the phylogenetic signal for each gene by leveraging the bootstrap trees generated in the previous step (i.e., `BC01-CellPhy.raxml.bootstraps`), alongside the transformed gene expression data (`BC01_VSD.txt`).

#### Input files:
*Gene Expression Data:* `BC01_VSD.txt` – This file contains the transformed (VSD) gene expression profiles.
*Bootstrap Trees:* `BC01-CellPhy.raxml.bootstraps` – This file contains the 100 bootstrap phylogenetic trees generated from the CellPhy analysis.

#### Phylogenetic signal estimation process:
Phylogenetic signal will be measured by comparing the expression data across different cells in the context of their evolutionary relationships. We will use the *tree_physig* function to estimate the lambda score and its confidence intervals:

```r

# Load libraries
library(ape)
library(phytools)
library(sensiPhy)

# Read in bootstrap trees and gene expression data
multi=read.tree("BC01-CellPhy.raxml.bootstraps")
geneexp=read.table("BC01_VSD.tx", head=T, row.names=1)

# Create an empty list to store results
scAllele_phylosignal <- list()

# Iterate over each gene:
for (i in 1:dim(geneexp)[1]) {

# Extract expression values for the current gene
temp_matrix <- t(geneexp[i, ])
t=as.data.frame(temp_matrix[,1])
names(t)=colnames(temp_matrix)
t$Gene=t[,1]

# Estimate phylogenetic signal for the gene
res=tree_physig(trait.col="Gene", data=t, phy=multi, n.tree="all", method="lambda")
scAllele_phylosignal[[i]] <- cbind(colnames(temp_matrix), res$stats$mean[1], res$stats$CI_low[1], res$stats$CI_high[1],res$stats$mean[2], res$stats$CI_low[2], res$stats$CI_high[2])
}

# Convert to data frame, rename columns and export results:
phylosignal_df=as.data.frame(matrix(unlist(scAllele_phylosignal), ncol=7, byrow=T))
names(phylosignal_df)=c("GeneID","lambda_mean","lambda_CI_low","lambda_CI_high","pvalue_mean","pvalue_CI_low","pvalue_CI_high")
write.table(phylosignal_df, "BC01-Phylosignal", quote=F, row.names=F, sep="\t")
```

