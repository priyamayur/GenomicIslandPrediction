# TREASUREISLAND

TreasureIsland python package is a machine learning-based Genomic Island prediction software, that uses an unsupervised representation of DNA for its prediction.

## Dependency :

Python >= 3.7

## Installation and Usage:

#Use pip to install the package :

TreasureIsland can be installed by python package management system "pip" :
   pip install treasureisland --pre

The treasureisland package can be used to find predictions and downloaded as csv, xlsx, txt files demonstrated in [TreasureIsland package](#TreasureIsland-package)

#Locally install package and use it:
   git clone https://github.com/priyamayur/GenomicIslandPrediction.git
   python -m pip install -e GenomicIslandPrediction
   cd GenomicIslandPrediction
   python test_ti <DNA file> 


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

The sample outputs can be found in the repository - output_NC_002620.2.txt, output_NC_002620.2.csv, output_NC_002620.2.xlsx 





