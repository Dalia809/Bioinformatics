# Week 3

## Q2

scp /Users/daliaalsaihati/Desktop/ncbi_dataset 2.zip  
alsaihdh@ilogin.ibex.kaust.edu.sa:/home/Alsaihdh/

### Log into Ibex

ssh Alsaihdh@ibex.kaust.sa

cd /home/Alsaihdh/

### Unzip file

unzip ncbi_dataset 2.zip

## Q3 Answer

To find the line with the smallest genome size and print the full line:

### Command

```bash
sort -k2,2n data.tsv | head -n 1
```

### Output

```
Aquifex aeolicus VF5    1.591   1,551,335       2       43.5
```

To output only the genome size for the smallest genome:

### Command

```bash
sort -k2,2n data.tsv | head -n 1 | awk '{print $(NF-3)}'
```

### Output

```
1.591
```

To find the line with the largest genome size and print the full line:

### Command

```bash
sort -k2,2n data.tsv | tail -n 1
```

### Output

```
Vibrio cholerae O1 biovar El Tor str. N16961    4.033   2,961,149       2       47.5
```

To output only the genome size for the largest genome:

### Command

```bash
sort -k2,2n data.tsv | tail -n 1 | awk '{print $(NF-3)}'
```

### Output

```
4.033
```

***

## Q4 Answer

To find the number of genomes where the species name contains at least two "c"s:

### Command

```bash
tail -n +2 data.tsv | grep -i "c.*c" | wc -l
```

### Output

```
7
```

To find the number of genomes where the species name contains at least two "c"s, but does not contain the word "coccus":

### Command

```bash
tail -n +2 data.tsv | grep -i "c.*c" | grep -iv "coccus" | wc -l
```

### Output

```
5
```

***

## Q5 Answer

To find all genome files (FASTA) larger than 3 megabytes:

### Command

```bash
find data/ncbi_dataset/data -name "*.fna" -size +3M | wc -l
```

### Output

```
3
```

***

# Week 4+5

## Q1 Answer

### Bash Script

```bash
#!/bin/bash
peptide="KVRMFTSELDIMLSVNGPADQIKYFCRHWT*"

num_amino_acids=${#peptide} -1
echo "Number of amino acids in the peptide: $num_amino_acids"

num_bases=$(( ($num_amino_acids) * 3 )+3) 
echo "Number of bases in the ORF: $num_bases"
```

### Output

```
Number of amino acids in the peptide: 30
Number of bases in the ORF: 93
```

***

## Q2 Answer

### Bash Script

```bash
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

I first used this script to run Prodigal on all GCA files:

```bash
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

Then I used this script to count the CDS's in the .gbk files generated by Prodigal:

```bash
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

Similarly, I first used this script to run Prokka on all GCA files:

```bash
#!/bin/bash

for genome_dir in GCA_*; do
  genome_file=$(find "$genome_dir" -name "*_genomic.fna")
  
  if [[ -f "$genome_file" ]]; then
    prokka_output_dir="${genome_dir}_prokka_output"
    prokka --outdir "$prokka_output_dir" --prefix "${genome_dir}" "$genome_file"
    
    prokka_output_file="$prokka_output_dir/${genome_dir}.txt"
    
    if [[ -f "$prokka_output_file" ]]; then
      cds_count=$(grep -w "CDS" "$prokka_output_file" | awk '{print $2}')
      echo "Genome: $genome_file - Number of CDS annotated: $cds_count"
    else
      echo "Prokka output file not found for $genome_file"
    fi
  else
    echo "No genomic file found in $genome_dir"
  fi
done
```

Then I collected all CDS counts from within the `.txt` file in the `prokka_output` folder generated for each GCA file as follows:

```bash
#!/bin/bash

for prokka_dir in *_prokka_output; do
  txt_file=$(find "$prokka_dir" -name "*.txt")
  
  if [[ -f "$txt_file" ]]; then
    cds_count=$(grep "CDS:" "$txt_file" | awk '{print $2}')
    echo "Directory: $prokka_dir - CDS count: $cds_count"
  else
    echo "No .txt file found in $prokka_dir"
  fi
done
```

### Output

```
Directory: GCA_000006745.1_prokka_output - CDS count: 3589
Directory: GCA_000006825.1_prokka_output - CDS count: 2028
Directory: GCA_000006865.1_prokka_output - CDS count: 2383
Directory: GCA_000007125.1_prokka_output - CDS count: 3150
Directory: GCA_000008525.1_prokka_output - CDS count: 1577
Directory: GCA_000008545.1_prokka_output - CDS count: 1861
Directory: GCA_000008565.1_prokka_output - CDS count: 3245
Directory: GCA_000008605.1_prokka_output - CDS count: 1001
Directory: GCA_000008625.1_prokka_output - CDS count: 1771
Directory: GCA_000008725.1_prokka_output - CDS count: 892
Directory: GCA_000008745.1_prokka_output - CDS count: 1058
Directory: GCA_000008785.1_prokka_output - CDS count: 1504
Directory: GCA_000027305.1_prokka_output - CDS count: 1748
Directory: GCA_000091085.2_prokka_output - CDS count: 1056
```

## Q5 Answer


### Command to combine all `.gbk` files into one:
```
$ find _prokka_output -name ".gbk" -exec cat {} + >> combined_annotations.gbk
```
### Command to extract unique gene names from the combined `.gbk` file:
```
$ grep "/gene=" combined_annotations.gbk | cut -d'=' -f2 | tr -d '"' | sort | uniq > unique_gene_names.txt
```
### Command to display the first five unique gene names:
```
$ head -n 5 unique_gene_names.txt
```
### Output:
```
aaaT aaeA aaeA_1 aaeA_2 aaeB
```

