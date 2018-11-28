# ProtMapPep
Map the peptides to its corresponding protein sequence and locate the modification sites

## Usage

    python Pepmap.py -c "config file" -s "species" -i "infile" -o "outfile" -p "pepindex" -m "modindex" -a "prot_acc_index"

#### Description 
    -c, --conf         config file used for searching
    
    -s,--sp            species name(as in database): case sensitive
  
    -i,--ifile         infile name with extension
  
    -o,--ofile         outfile name with extension
  
    -p,--pep           Column index of peptide in infile

    -m,--mod           column index of modified peptide in infile
  
    -a,--acc           column index of protein accession in infile

    -h, --help         show this help message and exit

#### example

    python Pepmap.py -c conf.txt -s HUMAN -i infile.txt -o outfile.txt -p 2 -m 5 -a 7
    

## Information written to outfile

Results will be wirtten as tab delimited to the specified outfile name. Last three columns of your outfile contains the mapped information as shown below

    Column[-3]          | Column[-2]                           | Column[-1]
    --------------------|--------------------------------------|-----------
    Peptide start in UP | Modification__modified position in UP| Indels(Yes/No)
 
