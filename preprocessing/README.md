# Preprocessing

This step-by-step is based on [DIY Transcriptomics](https://diytranscriptomics.com/data). 
We thanks all efforts made by [Beiting lab](https://hostmicrobe.org/) to develop the 
transcriptomics field.

---

## Download your fastq files from Illumina basespace

Basespace cli download in the 
[Basespace CLI messy site](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview)
(The -O option doesn't work without mkdir before), this simplify the process: 
``` Linux
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
sudo mv bs /bin/bs
sudo chmod u+x /bin/bs
```
### Check your project id
> bs list projects 
\or 
> bs list run

> bs download projects -id XXXX --exclude=* --include=*fastq.gz

---

## Download the reference and build the index

Download the cDNA GRCh38 (release 112) from ensembl site:
https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/
>wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

Make the index based on transcripts file:
>kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa.gz

---

Run the scripts:

1. 1_run_qc.sh
2. 2_run_kallisto.sh
