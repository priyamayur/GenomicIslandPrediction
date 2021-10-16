# TREASUREISLAND

TreasureIsland python package is a machine learning-based Genomic Island prediction software, that uses an unsupervised representation of DNA for its prediction.

TreasureIsland was contructed from [Benbow dataset](https://github.com/priyamayur/GenomicIslandPrediction/tree/master/Benbow). 

## Dependency :

Python >= 3.7

## Installation:

### Option1 - Use pip to install the package :
TreasureIsland can be installed by python package management system "pip" :

    pip install treasureisland --pre 

### Option2 - Locally install package:
    git clone https://github.com/priyamayur/GenomicIslandPrediction.git
    python -m pip install -e GenomicIslandPrediction
    
    
## Usage:

The treasureisland package is used to find genomic island predictions which can be downloaded as csv, xlsx, txt files demonstrated in [TreasureIsland package](#TreasureIsland-package)

Or, run script locally to get predicitons quickly:

Clone the github repository if not cloned before:   

    git clone https://github.com/priyamayur/GenomicIslandPrediction.git
    cd GenomicIslandPrediction
    python test_ti <DNA file>     
    
### Input file:

DNA sequence files in formats - fasta(recommended) or genbank with a sequenceID.

example: >NC_002620.2 Chlamydia muridarum str. Nigg, complete sequence
CACATAGCAAAACACTCAAAGTTTTTCAGCAAAAAAGCTTGTTGAAAAAATTGTTGACCGCCTGTTCACA....

### Performance:

TreasureIsland takes 2-5 mins to run depending on the size of the input.

### Output :

Can be downloaded in csv, xlsx, txt formats.
The results are shown in the following format for each genomic island:
<sequenceID> <start> <end> <confidence>

example : NC_002620.2 1.0 130000.0 0.95597
    
The sample outputs can be found in the repository - output_NC_002620.2.txt, output_NC_002620.2.csv, output_NC_002620.2.xlsx     
    
### Testing:
    
Repository contains some [sample DNA files](https://github.com/priyamayur/GenomicIslandPrediction/tree/master/genome) that can be used to test the TreasureIsland.
    
example :
    
    cd GenomicIslandPrediction
    python test_ti genome/ecoli.fasta    


## TreasureIsland package:

import the sequence class from treasureisland package:

    from treasureisland.dna_sequence import sequence 

Instantiate the sequence with the DNA sequence file path as the argument. 
The DNA file used can be a fasta or genbank file.

    seq = sequence("C:/Users/USER/GenomicIslandPrediction/genome/bsub.fasta") # enter local path for sequence file

Get prediction data frame from sequence by running the predict method.

    pred = seq.predict(seqfile)

The predictions can be downloaded in text, csv, excel formats.

    seq.predictions_to_csv(pred)
    seq.predictions_to_excel(pred)
    seq.predictions_to_text(pred)




## Contact:

Feel free to contact at banerjee.p1104@gmail.com


