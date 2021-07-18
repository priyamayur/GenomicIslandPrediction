# TREASUREISLAND

TreasureIsland is a machine learning-based Genomic Island prediction software, that uses an unsupervised representation of DNA for its prediction.

Installation :

Use pip to install the package :

    pip install treasureisland

## Sample code:

Sample code can be found in test.py 

import the gi_driver from treasure island package:

    from treasureisland.gi_driver import gi_driver 

Instantiate the gi_driver with the DNA sequence file path as the argument. 
The DNA file used can be a fasta or genbank file.

    driver = gi_driver("C:/Users/USER/GenomicIslandPrediction/genome/bsub.fasta") # enter local path for sequence file

Get prediction data frame from gi_driver by running the get predictions.

    pred = driver.get_predictions()

The predictions can be downloaded in text, csv, excel formats.

    driver.predictions_to_csv(pred)

    driver.predictions_to_excel(pred)

    driver.predictions_to_text(pred)

The sample outputs can be found in the repository - output.txt, output.csv, output.xlsx 




