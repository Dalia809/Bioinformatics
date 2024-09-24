# Q2
scp /Users/daliaalsaihati/Desktop/ncbi_dataset 2.zip
alsaihdh@ilogin.ibex.kaust.edu.sa:/home/Alsaihdh/
### log into ibex
ssh Alsaihdh@ibex.kaust.sa

cd /home/Alsaihdh/

### unzip file
unzip ncbi_dataset 2.zip

# Q3 Answer 
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
# Q4 Answer 

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
 tail -n +2 data.tsv | grep -i "c.*c" | grep -iv "coccus" | wc -l
 ```

### output
```
5
```
***
# Q5 Answer 
To find all genome files (FASTA) larger than 3 megabytes
### command
```
 find data/ncbi_dataset/data -name "*.fna" -size +3M | wc -l
```

### output
```
3
```
