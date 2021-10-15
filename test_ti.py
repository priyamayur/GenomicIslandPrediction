#!/usr/bin/env python
import sys
from treasureisland.Program import Program

def main(seqfile):
    driver = Program(seqfile)
    pred = driver.get_predictions()
    driver.predictions_to_csv(pred)
    driver.predictions_to_excel(pred)
    driver.predictions_to_text(pred)

if __name__ == '__main__':
    main(sys.argv[1])

