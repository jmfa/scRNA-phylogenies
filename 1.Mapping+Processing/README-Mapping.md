## Step 1: Mapping Raw Reads
The first step in our pipeline is mapping raw reads to the **human genome reference** (e.g., hg38). Here, we assume that all data has already been retrieved from SRA and that quality control (QC) checks have been performed to trim low-quality reads and remove adapter sequences.

For mapping, we use STAR aligner (v2.7.9a), followed by MarkDuplicates from Picard Tools (v2.25.5) to remove PCR duplicates.

### Defining Input Samples
Let's assume we have a sample file (BC01-Samples) with the following structure:
