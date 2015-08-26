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

## Analysis
**Getting the TIS coordinations in the human genome annotation**

Extract start codon features from gencode.v19.annotation.gtf: 
```shell
awk '$3 ~/start_codon/' gencode.v22.annotation.gtf > gencode.v22.annotation_startCodon.gtf
```
Extend start condon into TIS (-9 to +6) and move strand column to the end:
```python
import re
f = open("gencode.v22.annotation_startCodon.gtf","rU")
output = open("gencode.v22.annotation_TIS_raw_ss.gtf","w")
file_list = f.readlines()
for line_list in file_list:
    line_list = line_list[:-1]
    line = re.split("\t",line_list)
    strand = line[6]
    if strand == "+":
        start = int(line[3])
        start-=9
        line[3] = str(start)
        end = int(line[4])
        end+=3
        line[4] = str(end)
    else:
        start = int(line[3])
        start-=3
        line[3] = str(start)
        end = int(line[4])
        end+=9
        line[4] = str(end)
    line.remove(line[6])
    line.append(strand)
    newline =  '\t'.join(line)
    newline += "\n"
    output.write(newline)
print("done")
```
Modify the first column of the gtf file to match the vcf file (for the purpose of running bedtools):
```python
import re
file = open("gencode.v22.annotation_TIS_raw_ss.gtf")
output = open("gencode.v22.annotation_TIS_ss.gtf","w")
file_list = file.readlines()
for line in file_list:
    newline = line[3:]
    output.write(newline)
print ("done")
```
Get strand info in the vcf file:
```python
import re
f = open("CosmicCodingMuts.vcf","rU")
output = open("CosmicCodingMuts_s.vcf","w")
file_list = f.readlines()
for line_list in file_list:
    if line_list[0] != "#":
        line = re.split("\t",line_list)
        s = re.split(";",line[7])
        strand = s[1][-1]
        newline = line_list[:-1]
        newline += "\t"
        newline += strand
        newline += "\n"
        output.write(newline)
    else:
        output.write(line_list)

print("done")
```
Run **bedtools intersect** to find mutations in the COSMIC database that occur in the TIS region:
```shell
bedtools intersect -a ~/kozak/data/15_06_16_gencode_v22/gencode.v22.annotation_TIS_ss.gtf -b CosmicNonCodingVariants.vcf -wa -wb > TIS_mutation_strand_ss.txt
```


s

s
s
s
s
s

s
s

sss


s
s
s



ss

s

s
s

s
s