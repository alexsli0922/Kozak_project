## Data sources & downloads

**Human genome**
* **gencode.v22.annotation.gtf**: Comprehensive gene annotation (CHR); Version 22
    * Downloaded from: http://www.gencodegenes.org/releases/current.html
* **hg38.fa**: Human genome sequence
    * Downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

**COSMIC mutation database** (release v73)
* To download:
    * ``` sftp "your_email_address"@sftp-cancer.sanger.ac.uk ```
    * ```sftp> get /files/grch38/cosmic/v73/VCF/name_of_file```
* **CosmicCodingMuts.vcf**: All coding mutations in COSMIC
* **CosmicNonCodingVariants.vcf**: All noncoding mutations in COSMIC
* **COSMIC Whole Genome** (coding & noncoding): Contains only mutations found with whole genome sequencing experiments. 
* 