#!/usr/bin/env python
import sys
from treasureisland.gi_driver import gi_driver


# driver = gi_driver("C:/Users/USER/GenomicIslandPrediction/genome/test.fasta") # enter local path for sequence file

def main(seqfile):
    driver = gi_driver(seqfile)
    pred = driver.get_predictions()
    driver.predictions_to_csv(pred)
    driver.predictions_to_excel(pred)
    driver.predictions_to_text(pred)

if __name__ == '__main__':
    main(sys.argv[1])

