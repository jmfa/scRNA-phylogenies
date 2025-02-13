## Step 4: Tree Reconstruction with CellPhy
Using the VCF file generated in the previous step, we can now run [CellPhy](https://github.com/amkozlov/cellphy) to infer the single-cell phylogeny.
Here, we use the default strategy to perform a full analysis, including tree search and bootstrapping:

```sh
./cellphy.sh FULL BC01-scAllele-Genotyped.Missingness.Indv.recode.vcf \
  --m GT16+FO+E \
  --bs-metric tbe \
  --prefix BC01-CellPhy \
  -bs-trees 100
```
> *Note*: The --m GT16+FO+E flag specifies the default genotype model (GT16), while --bs-metric tbe enables transfer bootstrap expectation (TBE) support values.

Once CellPhy completes, it will generate several output files. For this tutorial, we will focus on two key files:
- *BC01-CellPhy.raxml.support* → The single-cell phylogeny with bootstrap support (used for plotting).
- *BC01-CellPhy.raxml.bootstraps* → 100 bootstrap trees (used for phylogenetic signal estimation).
