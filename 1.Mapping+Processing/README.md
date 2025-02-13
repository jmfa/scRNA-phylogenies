## Step 1: Mapping Raw Reads
The first step in our pipeline is mapping raw reads to the **human genome reference** (e.g., hg38). Here, we assume that all data has already been retrieved from SRA and that quality control (QC) checks have been performed to trim low-quality reads and remove adapter sequences.

For mapping, we use STAR aligner (v2.7.9a), followed by MarkDuplicates from Picard Tools (v2.25.5) to remove PCR duplicates.

### Defining Input Samples
Let's assume we have a sample file (BC01-Samples) with the following structure:

```
$ head -n 5 BC01-Samples
SRR2973279    BC01_02   Tumor
SRR2973280    BC01_03   Tumor
SRR2973281    BC01_04   Tumor
SRR2973282    BC01_05   Tumor
SRR2973283    BC01_06   Tumor
```

Each row contains:
- SRA Accession ID (Column 1)
- Cell Identifier (Column 2)
- Cell type (Column 3)

We will use a simple while loop to extract these variables and use them as input for the Mapping_Processing.sh script:
```
while read file
do
    sraID=$(echo $file | cut -d " " -f 1)
    cellID=$(echo $file | cut -d " " -f 2)
    ./Mapping_Processing.sh 
done < BC01-Samples
```
This loop iterates through each line of the `BC01-Samples` file, extracts the relevant IDs, and passes them as arguments to the `Mapping_Processing.sh` script.
Once this step is complete, we can move on to the next section: `Variant Calling + Genotyping`.
