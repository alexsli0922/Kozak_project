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

**Kozak strength**
* Idea came from this paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299517/
* **Supplementary Table S2** contains TIS strength prediction for all permutations of TIS from -6 to +5
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
Get strand information in the vcf file to a seperate column (for the purpose of running bedtools):
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
Get mutation and TIS that are on the same strand of DNA:
```python
import re
f = open("TIS_mut_strand_ss.txt","rU")
f1 = open("CosmicMutantExport.tsv")
output = open("TIS_mut_strand_transcript.txt","w")
file_list = f.readlines()
file1_list = f1.readlines()
for line_list in file_list:
    line_list = line_list[:-1]
    line = re.split("\t",line_list)
    if line[8] == line[17]:
        line_list += "\n"
        output.write(line_list)
print "done"
```
Get the transcript ID numbers for all mutations:
```python
## get the transcipt id for all mutations
## use binary search
## get all mutation info as well
import re
f = open("TIS_mut_strand_ss.txt","rU")
f1 = open("CosmicMutantExport.tsv")
output = open("TIS_mut_strand_transcript.txt","w")
output1 = open("unfound_mutations.txt","w")
file_list = f.readlines()
file1_list = f1.readlines()
temp = [0] * len(file1_list)
tdict = {}
for tline_list in file1_list:
    tline = re.split("\t",tline_list)
    tdict[tline[12]] = [tline[1],tline[7],tline[8],tline[9],tline[10],tline[15]]
for line_list in file_list:
    line_list = line_list[:-1]
    line = re.split("\t",line_list)
    try:
        transcript = tdict[line[11]][0]
        psite = tdict[line[11]][1]
        ssite = tdict[line[11]][2]
        phist = tdict[line[11]][3]
        shist = tdict[line[11]][4]
        description = tdict[line[11]][5]
        line_list += "\t" + transcript + "\t" + psite + "\t" + ssite + "\t" + phist + "\t" + shist + "\t" + description + "\n"
        output.write(line_list)
    except KeyError:
        nl = line_list + "\n"
        output1.write(nl)
print "done"
```
Get only mutations with correctly matched transcription ID number:
```python
## check if gtf and vcf has the same gene name
import re
f = open("TIS_mut_strand_transcript.txt","rU")
output = open("TIS_mut_strand_tm.txt","w")
output1 = open("TIS_mut_strand_tum.txt","w")
file_list = f.readlines()
for line_list in file_list:
    line = re.split("\t",line_list)
    gtf= re.split(";",line[7])
    gname = gtf[1][16:-3]
    tname = line[18]
    if gname == tname:
        output.write(line_list)
    else:
        output1.write(line_list)
print "done"
```
Count number of mutation occurrence for each gene:
```python
## count number of occurrence for each gene
## using the name from gencode
import re
f = open("TIS_mut_strand_tm.txt","rU")
output = open("count_gene_tm.txt","w")
file_list = f.readlines()
names = {}
for line_list in file_list:
    line = re.split("\t",line_list)
    gtf= re.split(";",line[7])
    gname = gtf[4][12:-1]
    if gname in names.keys():
        names[gname] = names[gname] + 1
    else:
        names[gname] = 1
for n in names:
    output.write( n+"\t"+str(names[n])+'\n')
print len(names)
print "done"
```
Get counts of different characterizations of mutations: 
```python
## count number of occurrence for each gene
## using the name from gencode

import re
f = open("TIS_mut_strand_tm.txt","rU")
file_list = f.readlines()
muts = {}
for line_list in file_list:
    line = re.split("\t",line_list)
    m = line[23][:-1]
    if m in muts.keys():
        muts[m] += 1
    else:
        muts[m] = 1
for t in muts:
    print t + "\t" + str(muts[t]) + "\n"
print "done"
```
Get all substitution mutations and record number of mutations at each position from -9 to +6:
```python
## count number of mutation at each position
import re
f = open("TIS_mut_strand_tm.txt","rU")
output = open("count_position_tm_1.txt","w")
file_list = f.readlines()
position = [0]*15
for line_list in file_list:
    line = re.split("\t",line_list)
    types = line[23]
    if "Substitution" in types or "substitution" in types:
        if line[8] == "+":
            for i in range (0,len(line[12])):
                if line[12][i] != line[13][i]:
                    pos = int(line[10]) - int(line[3]) + i
                    try:
                        position[pos] +=1
                    except IndexError:
                        print pos
                        print line_list
        else:
            for i in range (0,len(line[12])):
                if line[12][i] != line[13][i]:
                    pos = int(line[4]) - int(line[10]) - i
                    try:
                        position[pos] +=1
                    except IndexError:
                        print pos
                        print line_list
for i in range (0,len(position)):
    output.write(str(i) + "\t" + str(position[i]) + "\n")
print "done"
```

## Noncoding mutations
For noncoding mutations, there is no strand information. To reduce the possibility of redundancy in the reported mutations, same mutations at the same position are counted as one mutation collectively.
```python
import re
f = open("TIS_mutation_strand_ss.txt","rU")
output = open("TIS_mut_strand_red.txt","w")
file_list = f.readlines()
temp_p = ""
temp_r = ""
temp_a = ""
for line_list in file_list:
    line = re.split("\t",line_list)
    pos = line[10]
    ref = line[12]
    alt = line[13]
    if temp_p != pos or temp_r !=ref or temp_a != alt:
        output.write(line_list)
        temp_p = pos
        temp_r = ref
        temp_a = alt
print "done"
```
## Kozak strength prediction
I used dinucleotide PWM from the kozak paper (link to the paper)
* Get the genomic sequences from -6 to +5 for both the reference sequences and the mutated sequences.
* Use the Kozak strength table to get the predicted kozak strengths

Get only mutations from -6 to +5 excluding the start codon
```python
import re
f = open("TIS_mut_strand_tm.txt","rU")
output = open("TIS_mut_kozak.txt","w")
file_list = f.readlines()
for line_list in file_list:
    line = re.split("\t",line_list)
    ref = line[12]
    alt = line[13]
    if len(ref) == len(alt):
        if line[8] == "+":
            for i in range (0,len(line[12])):
                if line[12][i] != line[13][i]:
                    pos = int(line[10]) - int(line[3]) + i
        else:
            for i in range (0,len(line[12])):
                if line[12][i] != line[13][i]:
                    pos = int(line[4]) - int(line[10]) - i
        if pos > 2 and pos < 14 and pos != 9 and pos != 10 and pos != 11:
            output.write(line_list)
print "done"
```
Create a bed file with coordinates of mutated TIS sites 
```python
## start codon to TIS
## move strand info to the last column
import re
f = open("/Users/xuanyi/kozak/analysis/gencode_v22/TIS_mut_kozak.txt","rU")
output = open("kozak_mut_cosmic.bed","w")
file_list = f.readlines()
for line_list in file_list:
    line = re.split("\t",line_list)
    strand = line[8]
    if strand == "+":
        start = int(line[3])
        start+=2
        line[3] = str(start)
        end = int(line[4])
        end-=1
        line[4] = str(end)
    else:
        start = int(line[3])
        line[3] = str(start)
        end = int(line[4])
        end-=3
        line[4] = str(end)
    output.write("chr" + line[0] + "\t" + line[3] + "\t" + line[4]  + "\n")
print("done")
```
Use **bedtools getfasta** to get the TIS sequences:
* ```~/kozak/bin/bedtools2/bin/bedtools getfasta -fi /Users/xuanyi/kozak/data/15_07_13_gencode_fasta/hg38.fa -bed kozak_mut_cosmic.bed -fo kozak_mut_cosmic.fa```

To capitalize all bases:
```python
import re
f = open("kozak_mut_cosmic.fa","rU")
output = open("kozak_ori_cosmic.txt","w")
file_list = f.readlines()
for line_list in file_list:
    if line_list[0] != ">":
        temp = line_list.upper()
        output.write(temp)
print("done")
```
Get a file in this format: chr	start	end	mut	ref	alt
```python
## start codon to TIS
## move strand info to the last column
import re
f = open("/Users/xuanyi/kozak/analysis/gencode_v22/TIS_mut_kozak.txt","rU")
output = open("kozak_mut_pos_cosmic.txt","w")
file_list = f.readlines()
for line_list in file_list:
    line = re.split("\t",line_list)
    strand = line[8]
    if strand == "+":
        start = int(line[3])
        start+=2
        line[3] = str(start)
        end = int(line[4])
        end-=1
        line[4] = str(end)
    else:
        start = int(line[3])
        line[3] = str(start)
        end = int(line[4])
        end-=3
        line[4] = str(end)
    output.write("chr" + line[0] + "\t" + line[3] + "\t" + line[4] + "\t" + line[10] + "\t" + line[12] + "\t" + line[13]  + "\n")
print("done")
```
Get the all the mutated TIS sequences:
```python
import re
f = open("/Users/xuanyi/kozak/analysis/gencode_v22/TIS_mut_kozak.txt","rU")
f1=open("kozak_ori_cosmic.txt","rU")
f2=open("kozak_mut_pos_cosmic.txt","rU")
output = open("kozak_mut_seq_cosmic.txt","w")
file_list = f.readlines()
file_list1 = f1.readlines()
file_list2 = f2.readlines()
names = {}
for i in range (0,len(file_list)):
    line = re.split("\t",file_list[i])
    if file_list2[i][-1] == "\n":
        file_list2[i] = file_list2[i][:-1]
    line2 = re.split("\t",file_list2[i])
    gtf= re.split(";",line[7])
    gname = gtf[4][12:-1]
    if gname in names.keys():
        newline = ""
        pos = int(line2[3]) - int(line2[1]) -1
        j = 0
        while j < len(file_list1[i]):
            if j == pos:
                for k in line2[5]:
                    newline += k
                    j += 1
            else:
                newline += file_list1[i][j]
                j += 1
        if newline[-1] == "\n":
            newline = newline[:-1]
        names[gname].append(newline)
    else:
        names[gname] = []
        if file_list1[i][-1] == "\n":
            nl = file_list1[i][:-1]
        names[gname].append(nl)
        newline = ""
        pos = int(line2[3]) - int(line2[1]) -1
        j = 0
        while j < len(file_list1[i]):
            if j == pos:
                for k in line2[5]:
                    newline += k
                    j += 1
            else:
                newline += file_list1[i][j]
                j += 1
        if newline[-1] == "\n":
            newline = newline[:-1]
        names[gname].append(newline)
for n in names:
    output.write( n+ "\t")
    for m in names[n]:
        output.write(m + "\t")
    output.write("\n")
print "done"
```
Get the kozak score for all TIS:
```python
import re
f = open("kozak_mut_seq.txt","rU")
f1=open("kozak_strength.txt","rU")
output = open("kozak_score_test.txt","w")
file_list = f.readlines()
file_list1 = f1.readlines()
sdict = {}
for line_list1 in file_list1:
    sline = re.split("\t", line_list1)
    sdict[sline[0]] = sline[1]
for i in range (0,len(file_list)):
    line = re.split("\t",file_list[i])
    output.write(line[0] + "\t" + str(len(line)-3) + "\t")
    j = 1
    while j < len(line) - 1:
        temp = line[j]
        if temp[6:9] == "ATG":
            newline = ""
            for k in range(0,len(temp)):
                if temp[k] == "T":
                    newline += "U"
                else:
                    newline += temp[k]
        else:
            newline = ""
            k = len(temp)
            while (k):
                t = k-1
                if temp[t] == "A":
                    newline+="U"
                if temp[t] == "T":
                    newline+="A"
                if temp[t] == "C":
                    newline+="G"
                if temp[t] == "G":
                    newline+="C"
                k -= 1
        try:
            score = sdict[newline]
            output.write(score + "\t")
        except KeyError:
            print temp
        j += 1
    output.write("\n")
print "done"
```
Get the difference between ref and mut for each genes
```python
## start codon to TIS
## move strand info to the last column
import re
f = open("kozak_score.txt","rU")
output = open("diff_score.txt","w")
file_list = f.readlines()
for line_list in file_list:
    line = re.split("\t",line_list)
    ref = line[2]
    diff = 0
    for i in range (3,len(line)-1):
        mut = line[i]
        try:
            temp = int(ref) - int(mut)
        except ValueError:
            continue
        if temp > diff:
            diff = temp
    output.write(line[0] + "\t" + str(diff)+ "\n")
print("done")
```
