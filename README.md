# MSOAR 2.0 Ortholog Assignment System
## version 2.0, Oct 29, 2009

# What is MSOAR and MSOAR 2.0

MSOAR is a high-throughput system for assigning orthologs between
closely related species on a genome scale. It attempts to reconstruct
the evolutionary history of input genomes in terms of genome rearrangement
and gene duplication events. However, MSOAR is only able to deal with
randomly duplicated genes. 

MSOAR 2.0 incorporates tandem duplications into MSOAR, and combines gene 
phylogeny and genome rearrangement together to assign orthologs. For a pair 
of input genomes, MSOAR 2.0 first focuses on the tandemly duplicated genes 
of each genome and tries to identify among them those that were duplicated 
after the speciation (i.e., the so-called inparalogs), using a phylogenetic
tree reconciliation method. Then it removes most of the duplicated inparalogs 
in tandem position on each genome before MSOAR is invoked. Finally, it goes 
through a post-processing step to further improve the prediction accuracy of
the output ortholog pairs.


# Reference

If you use MSOAR 2.0, please cite the following paper:

```
Guanqun Shi, Liqing Zhang, Tao Jiang (2009)
  MSOAR 2.0: Incorporating tandem duplications into ortholog assignment
  based on genome rearrangement.
```


# Directories and Files

Directory:   Programs/   Software and programs used in the pipeline of MSOAR 2.0
Directory:   tools/      Tools used in MSOAR 2.0
Directory:   Examples/   Examples of input data files
File:        MSOAR2.0    Executable script file to run MSOAR 2.0
File:        README      This file


Programs/
----------------------------
Subdirectory:

 * `blast-2.2.13/` BLAST program package (Since the blast programs are too large to
 be included in the package, the user may need to download the programs by
 themselves and intall the programs under the current directory.)

 * `MAFFT/` MAFFT program package (Since the MAFFT programs are too large to be
 included in the package, the user may need to download the programs by
 themselves and intall the programs under the current directory.)

 * `mcl-06-058/` MCL program package

 * `MSOAR/` MSOAR program package

 * `PAL2NAL/` PAL2NAL program package phylip-3.67/ Phylip program package

Source code files:

```
    AddOutgroup.cpp        Add an outgroup to the distance matrix
    FamilySeparator.cpp    Extract gene information for each gene family
    GetTopHits.cpp         Extract the top hits from the output file of BLASTp
    MapToPhylip.cpp        Map gene names to Phylip input format
    Postprocessing.cpp     Post-processing step in MSOAR 2.0
    RemoveInparalogs.cpp   Remove inparalogs in TAGs
    NormalizeScores.cpp    Normalize the bit scores for the top hits
    TagGenerator.cpp       Find the inparalogs in TAGs for each gene family
    TreeTran.cpp           Use the specified outgroup to root the gene tree
```


tools/
-----------------------------
This directory contains some tools written in perl that are used in the
MSOAR2.0 script file.


Examples/
-----------------------------
```
G1.pep    Genome1   peptide     sequence
G2.pep    Genome2   peptide     sequence
G1.nuc    Genome1   nucleotide  sequence
G2.nuc    Genome2   nucleotide  sequence
G1.info   Genome1   gene        positional information
G2.info   Genome2   gene        positional information
```

## File Format:

### \*.pep files (G1.pep):
```
>gene_id
MVTEFIFLGLSDSQELQTFLFMLFFVFYGG
```

### \*.nuc files (G1.nuc):
```
>gene_id
ATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTCCTATTTATGTTGTTTTTTGTAT
```

### \*.info files:
```
gene_id        gene_symbol        chromosome        sign        start_position
```


# How to run MSOAR 2.0

MSOAR 2.0 is written on the platform of Linux and tested on Linux (version 
CentOS-5), we haven't tested it on other platforms. So it is better for the 
user to install and run MSOAR 2.0 on Linux. 

Since MSOAR 2.0 is an integrated pipeline for ortholog assignment, it uses
some other software and programs as a subroutine. So the user needs to install
the following software and programs before running MSOAR2.0.


Install Other Software
-----------------------
```
    blast-2.2.13/                    
    MAFFT/        
    mcl-06-058/                    
    MSOAR/                         
    PAL2NAL/                
    phylip-3.67/
```

In order to install these software, the user needs to download them under the 
directory /Programs/. Some of them need to be compiled, for example, to compile
the MCL program, simply type:

```
$ cd mcl-06-058/
$ make
```

# NOTICE:
# All software have already been included in the package except the blast-2.2.13/
# and MAFFT/ due to their large sizes. User needs to install these two software 
# under the corresponding directories before running MSOAR 2.0.


Compile Programs used in MSOAR 2.0
----------------------------------

Go to the directory /Programs/, simply type:

```
$ make
```

to compile all the programs used in MSOAR 2.0.



Run MSOAR 2.0
-------------

After installing all the software and compiling all the programs used in 
MSOAR 2.0, the user need to go to the directory containing the MSOAR2.0 file
and run the MSOAR2.0 executable script as follows:

```
$ MSOAR2.0 G1.pep G2.pep G1.nuc G2.nuc G1.info G2.info
```

To run MSOAR 2.0, the MSOAR2.0 script needs six input data files, which are
described below:

 - `G1.pep` and `G2.pep` are the input peptide sequence files for each genome.

 - `G1.nuc` and `G2.nuc` are the corresponding nucleotide sequence files for each genome.

 - `G1.info` and `G2.info` are files containing positional information of genes on each genome.

Formats of the input files are described in the Examples/ section, and examples are 
given in the directory Examples/.

The output of MSOAR 2.0 is stored in the file MSOAR2_result after running the script. 



# MSOAR 2.0 Website

http://msoar.cs.ucr.edu/MSOAR2.0/


# Contact Info

If you have any questions or comments about MSOAR 2.0, please contact me at:


  Guanqun Shi  
  Ph.D. student  
  CSE Department, Univ. of California  
  Riverside, CA 92507  
  Email: gshi@cs.ucr.edu  
