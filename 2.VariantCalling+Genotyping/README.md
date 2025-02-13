## Step 2: Variant Calling & Genotyping

### Calling Variants with scAllele
Once the BAM files are ready, we can proceed with variant calling using [scAllele](https://github.com/gxiaolab/scAllele/).
First, we generate a list of all deduplicated BAM files:

```sh
ls *.dedup.bam > bamList
```
Now, we run scAllele using the following command:

```sh
./scAllele -b bamList -o BC01-scAllele -n 46 -g ${reference}/hg38.fa 
```

> **Note:** Adjust -n 46 to match the number of threads available for computation.


### Single-Cell Genotyping
Once variant calling is completed, we perform genotyping using a read count-based strategy (only SNVs are considered).
We iterate over each cell ID and extract single-cell genotypes:

```
while read cellID
do
    vcftools --vcf BC01-scAllele.vcf --remove-indels --remove-filtered-all --indv $cellID --recode --out temp

    # Create VCF header
    printf "##fileformat=VCFv4.2\n" > temp.head
    printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus genotype\">\n" >> temp.head
    printf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of reads overlapping the variant position\">\n" >> temp.head
    printf "##FORMAT=<ID=AC,Number=R,Type=Integer,Description=\"Allele counts\">\n" >> temp.head
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" "$cellID" >> temp.head

    # Genotype sites based on read counts
    grep -v "#" temp.recode.vcf | cut -f1,2,4,5,10 | sed 's/:/\t/g' | awk 'BEGIN{FS="\t"}{
        if ($6!=".") print $0
    }' | awk 'BEGIN{FS="\t"}{
        if ($7>=5 && $8>=2 && $9!=0) print $1"\t"$2"\t.\t"$3"\t"$4"\t10\tPASS\t.\tGT:DP:AC\t0/1:"$7":"$9","$8;
        else if ($7>=5 && $8==0) print $1"\t"$2"\t.\t"$3"\t"$4"\t10\tPASS\t.\tGT:DP:AC\t0/0:"$7":"$9","$8;
        else if ($7>=5 && $9==0) print $1"\t"$2"\t.\t"$3"\t"$4"\t10\tPASS\t.\tGT:DP:AC\t1/1:"$7":"$9","$8
    }' > temp.tail

    # Combine header and genotyped data
    cat temp.head temp.tail > ${cellID}.genotyped-scAllele.vcf
    rm temp.*

    # Compress and index VCF files
    bgzip ${cellID}.genotyped-scAllele.vcf
    tabix -p vcf ${cellID}.genotyped-scAllele.vcf.gz
done < BC01-CellIDs
```

### Merging Individual VCFs
Once all cells are genotyped, we merge the individual VCF files into a single dataset:

```sh
ls *.genotyped-scAllele.vcf.gz > mergeList
bcftools merge -l mergeList -O v -o BC01-scAllele-Genotyped.vcf
```

## Step 3: Low-Quality Variant & Cell Filtering

To ensure high-quality data, we apply filtering to remove sites and cells with excessive missing data. The thresholds used (as described in our methods section) are:
- Per-site missingness: Maximum 75% per cell type (examples of different cell types: Tumor and Healthy).
- Per-cell missingness: Maximum 75% genotypes missing.

### Filtering Low-Quality Variants (Per Site)
We filter sites separately for Tumor and Healthy cells:

```sh
vcftools --vcf BC01-scAllele-Genotyped.vcf --keep BC01-nonTumor_samples --recode --out temp.nonTumor
vcftools --vcf BC01-scAllele-Genotyped.vcf --keep BC01-Tumor_samples --recode --out temp.Tumor

vcftools --vcf temp.nonTumor.recode.vcf --extract-FORMAT-info GT --out temp.nonTumor
vcftools --vcf temp.Tumor.recode.vcf --extract-FORMAT-info GT --out temp.Tumor

exp=0.75
nb_tumor=$(wc -l < BC01-Tumor_samples)
nb_nontumor=$(wc -l < BC01-nonTumor_samples)
min_tumor=$(echo "$exp * $nb_tumor" | bc)
min_nontumor=$(echo "$exp * $nb_nontumor" | bc)

cat temp.nonTumor.GT.FORMAT | sed 's#0/0#R#g; s#0/1#M#g; s#1/1#M#g; s#\./\.#N#g' | \
awk 'BEGIN{OFS="\t"} {
    countM=0; countN=0;
    for(i=3; i<=NF; i++) {if($i=="M") countM++; if($i=="N") countN++;}
    print $1, $2, countM, countN
}' | awk -v min=$min_nontumor '{if ($4<=min) print $1, $2}' > pos_keep_nonTumor

cat temp.Tumor.GT.FORMAT | sed 's#0/0#R#g; s#0/1#M#g; s#1/1#M#g; s#\./\.#N#g' | \
awk 'BEGIN{OFS="\t"} {
    countM=0; countN=0;
    for(i=3; i<=NF; i++) {if($i=="M") countM++; if($i=="N") countN++;}
    print $1, $2, countM, countN
}' | awk -v min=$min_tumor '{if ($4<=min) print $1, $2}' > pos_keep_Tumor

cat pos_keep_nonTumor pos_keep_Tumor | sort | uniq -c | awk '{if ($1==2) print $2, $3}' > pos_PASS

vcftools --vcf BC01-scAllele-Genotyped.vcf --positions pos_PASS --recode --out BC01-scAllele-Genotyped.Missingness
rm temp.* pos_*
```

### Filtering Low-Quality Cells
Next, we remove cells with more than 75% missing genotypes:

```sh
vcftools --vcf BC01-scAllele-Genotyped.Missingness.recode.vcf --missing-indv --out BC01
awk '$5<=0.75 {print $1}' BC01.imiss > indv_keep
vcftools --vcf BC01-scAllele-Genotyped.Missingness.recode.vcf --keep indv_keep --recode --out BC01-scAllele-Genotyped.Missingness.Indv
```

## Final Thoughts
At this stage, we have a filtered VCF file containing high-quality variants and cells, ready for downstream phylogenetic analysis.
Next, we proceed with **Phylogenetic Tree Reconstruction**.

