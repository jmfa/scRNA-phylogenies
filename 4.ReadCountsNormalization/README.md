## Step 5: Read counts transformation

For all datasets analyzed in our study, the raw single-cell expression profiles were retrieved from the Gene Expression Omnibus (GEO) (e.g., [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)). The expression data is typically provided in a matrix format, where rows represent genes and columns represent individual samples. Below is an example of what the input matrices may look like:

```bash
head BC01-Samples
SRR2973279    BC01_02	Tumor
SRR2973280    BC01_03	Tumor
SRR2973281    BC01_04	Tumor
SRR2973282    BC01_05	Tumor
SRR2973283    BC01_06	Tumor

head BC1-RawCounts.tsv | cut -f 1-6 
BC01_02 BC01_03 BC01_04 BC01_05 BC01_06
DPM1  51 120  7  181  32
SCYL3 3  1 92 6  60
C1orf112  0 13 1  1  0
FUCA2 10  36 0 0 8
GCLC  3  5  108  2  53
```

The raw counts matrix (BC1-RawCounts.tsv) contains gene expression counts for each sample, and the sample metadata (BC01-Samples) provides information about each sample's ID and cell type (e.g., Tumor or Healthy). The raw counts can now be loaded into R for normalization and transformation using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) R package:

```r

# Load library
library(DESeq2)

# Load matrices
cts=read.table("BC1-RawCounts.tsv", head=T, row.names=1)
info=read.table("BC01-Samples", row.names=2, head=F)
names(info)=c("SRA","Type")

# Create DESeq2 dataset from raw counts
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = info, 
                              design = ~ Type)

# Run DESeq normalization
dds <- DESeq(dds)

# Apply variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Export results:
geneexp <- as.data.frame(assay(vsd))
write.table(geneexp, "BC01_VSD.txt", quote=F, row.names=T, sep="\t", col.names=T)
```

We can now proceed with phylogenetic signal estimation.
