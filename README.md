# Week 3
## Q2
scp /Users/daliaalsaihati/Desktop/ncbi_dataset 2.zip
alsaihdh@ilogin.ibex.kaust.edu.sa:/home/Alsaihdh/
### log into ibex
ssh Alsaihdh@ibex.kaust.sa

cd /home/Alsaihdh/

### unzip file
unzip ncbi_dataset 2.zip

## Q3 Answer 
To find the line with the smallest genome size and print the full line:

### command
```
sort -k2,2n data.tsv | head -n 1
```
### output
```
Aquifex aeolicus VF5    1.591   1,551,335       2       43.5
```

To output only the genome size for the smallest genome:
### command
```
sort -k2,2n data.tsv | head -n 1 | awk '{print $(NF-3)}'
```
### output
```
1.591
```

To find the line with the largest genome size and print the full line:
### command
```
sort -k2,2n data.tsv | tail -n 1
```
### output
```
Vibrio cholerae O1 biovar El Tor str. N16961    4.033   2,961,149       2       47.5
```

To output only the genome size for the largest genome:
### command
```
sort -k2,2n data.tsv | tail -n 1 | awk '{print $(NF-3)}'
```
### output
```
4.033
```
***
## Q4 Answer 

To find the number of genomes where the species name contains at least two "c"s
### command
```
tail -n +2 data.tsv | grep -i "c.*c" | wc -l
```
### output
```
7
```

To find the number of genomes where the species name contains at least two "c"s, but does not contain the word "coccus," we use the following shell command:
### command
```
 tail -n +2 ours.tsv | grep -i "c.*c" | grep -iv "coccus" | wc -l
 tail -n +2 data.tsv | grep -i "c.*c" | grep -iv "coccus" | wc -l
 ```

### output
```
5
```
***
## Q5 Answer 
To find all genome files (FASTA) larger than 3 megabytes
### command
```
 find data/ncbi_dataset/data -name "*.fna" -size +3M | wc -l
```

### output
```
3
```
***
# Week 4+5
## Q1 Answer 
### Bash Script
```
#!/bin/bash
peptide="KVRMFTSELDIMLSVNGPADQIKYFCRHWT"

num_amino_acids=${#peptide}
echo "Number of amino acids in the peptide: $num_amino_acids"

num_bases=$(( ($num_amino_acids) * 3 ))
echo "Number of bases in the ORF: $num_bases"
```
### Output
```
Number of amino acids in the peptide: 30
Number of bases in the ORF: 90
```
***
## Q2 Answer 
### Bash Script
```
#!/bin/bash

prodigal -i GCA_000006825.1/GCA_000006825.1_ASM682v1_genomic.fna -o genes_output.gbk -a proteins_output.faa

gene_count=$(grep -c "CDS" genes_output.gbk)

echo "Number of genes annotated: $gene_count"
```
### Output
```
Number of genes annotated: 2032
```
***
## Q3 Answer 
### Bash Script
I first used this script to run prodigal on all GCA files
```
#!/bin/bash

max_genes=0
best_genome=""

for genome_dir in GCA_*; do
  genome_file=$(find "$genome_dir" -name "*_genomic.fna")
  
  if [[ -f "$genome_file" ]]; then
    prodigal_output="${genome_dir}_genes.gbk"
    prodigal -i "$genome_file" -o "$prodigal_output" -a "${genome_dir}_proteins.faa"
    gene_count=$(grep -c "CDS" "$prodigal_output")
    echo "Genome: $genome_file - Number of genes: $gene_count"
    
    if [[ $gene_count -gt $max_genes ]]; then
      max_genes=$gene_count
      best_genome=$genome_file
    fi
  else
    echo "No genomic file found in $genome_dir"
  fi
done

echo "Genome with the highest number of genes: $best_genome with $max_genes genes."
```
The used this script to count the CDS's in the .gbk files generated by prodigal:
```
#!/bin/bash

for gbk_file in $(find . -name "*_genes.gbk"); do
  cds_count=$(grep -c "CDS" "$gbk_file")
  echo "File: $gbk_file - Number of CDS: $cds_count"
done
```
### Output
```
File: ./GCA_000006745.1_genes.gbk - Number of CDS: 3594
File: ./GCA_000006825.1_genes.gbk - Number of CDS: 2032
File: ./GCA_000006865.1_genes.gbk - Number of CDS: 2383
File: ./GCA_000007125.1_genes.gbk - Number of CDS: 3152
File: ./GCA_000008525.1_genes.gbk - Number of CDS: 1579
File: ./GCA_000008545.1_genes.gbk - Number of CDS: 1866
File: ./GCA_000008565.1_genes.gbk - Number of CDS: 3248
File: ./GCA_000008605.1_genes.gbk - Number of CDS: 1009
File: ./GCA_000008625.1_genes.gbk - Number of CDS: 1776
File: ./GCA_000008725.1_genes.gbk - Number of CDS: 897
File: ./GCA_000008745.1_genes.gbk - Number of CDS: 1063
File: ./GCA_000008785.1_genes.gbk - Number of CDS: 1505
File: ./GCA_000027305.1_genes.gbk - Number of CDS: 1748
File: ./GCA_000091085.2_genes.gbk - Number of CDS: 1063
```
***
## Q4 Answer 
### Bash Script

### Output
