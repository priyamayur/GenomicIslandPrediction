#!/usr/bin/env python
import sys
import argparse
from treasureisland.Predictor import Predictor


parser = argparse.ArgumentParser()
parser.add_argument('infile',type=argparse.FileType('r'))


def main(seqfile):
    seq = Predictor(seqfile)
    pred = seq.predict()
    seq.predictions_to_csv(pred)


if __name__ == '__main__':
    parser.parse_args()
    main(sys.argv[1])