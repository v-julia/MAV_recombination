# *Mamastrovirus* species are shaped by recombination and can be reliably distinguished in ORF1b genome region
Data and code for preprint "*Mamastrovirus* species are shaped by recombination and can be reliably distinguished in ORF1b genome region" by Y. Aleshina and A. Lukashev

Alignments, phylogenetic trees, and metadata are located in `data`.

`scripts/Figures_code.R` - code to reproduce the figures.


## Pipeline for generation of alignment of concatenated ORFs

- parser_gb.py
- get_orfs_coord.py
- split_genome.py
- TranslatorX
- join_al.py
- trimal v1.2rev57
- remove_similar.py

### parser_gb.py

Converts file with nucleotide sequences in genbank format to fasta format with sequence names in the following format: GenbankAccessionNumber_annotation1_annotation2_..._annotationN, where annotation1, ..., annotationN are metadata from genbank qualifiers indicated in *-features* option. 

Sequences of length  beyond or above defined threshold will not be included in output file. 

Saves output file to the same directory as the input file.


#### Usage

```
parser_gb.py [-h] -input INPUT_FILE -min MIN_LENGTH -max MAX_LENGTH -f FEATURES
                    FEATURES


  -input <str>          Path to input file in GenBank format

                        
  -min   <int>          Minimal length of sequence. Sequences shorter than min
                        length will not be included in the final dataset
                        
  -max   <int>          Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset
  -features   <str>     string with qualifiers to retrieve from GenBank annotation,\
                        e.g. 'country,host,collection_date'
```

### get_orfs_coord.py

Retrieves the coordinates of ORFs from file with nucleotide sequences in genbank format.

#### Usage
```
get_orfs_coord.py [-h] -input INPUT_FILE -orf_map ORF_MAP_FILE [-r]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in genbank format
  -orf_map ORF_MAP_FILE, --orf_map_file ORF_MAP_FILE
                        Csv-file with short codes for ORFs
  -r, --remove_exceptions
                        Removes sequences listed in file. Only for noroviruses yet
```


### split_genome.py

Excises ORF sequences from sequences in *input_file* in genbank format according to coordinates in *coord_file*

#### Usage
```
split_genome.py [-h] -input INPUT_FILE -coord COORD_FILE

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in fasta-format
  -coord COORD_FILE, --coord_file COORD_FILE
                        csv-file with coordinates of ORFs retrived from CDS field of genbank file
```


### gap_in_row.py

Removes sequences with many gaps in row from file in fasta format. Output file in fasta format is saved in the save folder as input file.

#### Usage
```
gap_in_row.py [-h] -input INPUT_FILE -gap_count GAP_COUNT

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file
  -gap_count GAP_COUNT, --gap_count GAP_COUNT
                        Amount of gaps in a row
```

### remove_similar.py

Calculates p-distances between nucleotide sequences pairs in loop. If p-distance < min_distance and p-distance > max_distance 
removes the sequence with higher serial number. 

Saves output file with reduced sequence set to the directory of input file.


#### Usage

```
remove_similar.py [-h] -input INPUT_FILE -min MIN_DISTANCE -max MAX_DISTANCE


  -input INPUT_FILE, --input_file INPUT_FILE
                        Path to input file in fasta format
  -min MIN_DISTANCE, --min_distance MIN_DISTANCE
                        Minimal pairwise distance between sequences. If
                        p-distance is lower than min_distance sequence with
                        higher serial number will be removed from the dataset
  -max MAX_DISTANCE, --max_distance MAX_DISTANCE
                        Maximal pairwise distance between sequences. If
                        p-distance is higher than max_distance sequence with
                        higher serial number will be removed from the dataset
```


### join_al.py

Concatenates alignments in fasta format

#### Usage
```
join_al.py [-h] -input_list INPUT_LIST -out_name OUT_NAME

optional arguments:
  -h, --help            show this help message and exit
  -input_list INPUT_LIST, --input_list INPUT_LIST
                        List with names of input file
  -out_name OUT_NAME, --out_name OUT_NAME
                        Name of output file
```