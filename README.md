# TreasureIsland

TreasureIsland python package is a machine learning-based Genomic Island prediction software, that uses an unsupervised representation of DNA for its prediction.

TreasureIsland is constructed from the [Benbow dataset](https://github.com/priyamayur/GenomicIslandPrediction/tree/master/Benbow). 

## Dependencies :

Python >= 3.7

## Installation:

### Option1 - Use pip:


    python -m pip install treasureisland
    
 if treasureisland is already installed :
 
    python -m pip install treasureisland --upgrade

### Option2 - Locally install package:
    git clone https://github.com/priyamayur/GenomicIslandPrediction.git
    python -m pip install -e GenomicIslandPrediction
    
    
## Usage:

### Option1 - Run TreasureIsland directly from commandline  :
Run TreasureIsland from commandline on your DNA fasta file (example DNA files provided [here](https://github.com/priyamayur/GenomicIslandPrediction/tree/master/genome)), output is given in csv format:

    treasureisland mypath/<DNA file>.fasta [-o <output_file_path>] [-ut <upper threshold value>] 
    
### Option2 - Run TreasureIsland from python :
The TreasureIsland package is used to find genomic island predictions which can be downloaded in csv, xlsx, txt file formats demonstrated in [TreasureIsland package](#Running-TreasureIsland-package-from-python)

### Input file:

DNA sequence files in fasta format with a sequenceID.

example: >NC_002620.2 Chlamydia muridarum str. Nigg, complete sequence
CACATAGCAAAACACTCAAAGTTTTTCAGCAAAAAAGCTTGTTGAAAAAATTGTTGACCGCCTGTTCACA....

### Performance:

TreasureIsland takes 2-5 mins to run depending on the size of the input.

### Output :

The results are shown in the following format for each genomic island:
<sequenceID> <start> <end> <probability of GEI>

example : NC_002620.2 1.0 130000.0 0.95597

### Upper Threshold:

User also has the ability to change the upper threshold value to change the precision and recall tradeoff. 
upper threshold is set to 0.80 by default.

Example : 

treasureisland ecoli.fasta -o gei_output -ut 0.95
Setting the upper threshold to 0.95 would increase the precision and decrease the recall performance.
    
### Testing:
    
A repository containing some [sample DNA files](https://github.com/priyamayur/GenomicIslandPrediction/tree/master/genome) that can be downloaded to test TreasureIsland. 
Note : github downloads fasta file in txt format (filename.fasta.txt). 
    
example :
    
    treasureisland ecoli.fasta -o gei_output -ut 0.95 


### Running TreasureIsland package from python:

import the Predictor class from treasureisland package:

    from treasureisland.Predictor import Predictor

Instantiate the sequence with the DNA sequence file path as the argument. 
The DNA file used can be a fasta file.

    seq = Predictor("<Path to DNA fasta file>/ecoli.fasta", "<output_file_path>") 

Optionally, change the upper threshold value.

    seq.change_upper_threshold(0.9)

Get prediction data frame from sequence by running the predict method.

    pred = seq.predict()

The predictions can be downloaded in text, csv, excel formats.

    seq.predictions_to_csv(pred)
    seq.predictions_to_excel(pred)
    seq.predictions_to_text(pred)

## Contact:

Feel free to contact at banerjee.p1104@gmail.com


