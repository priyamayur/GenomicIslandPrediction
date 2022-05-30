#!/usr/bin/env python
import sys
import argparse
from treasureisland.Predictor import Predictor


parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('r'))
parser.add_argument("-o", "--output", required=False, default="output", help="output directory name")
parser.add_argument("-ut", "--upperthreshold", required=False, default="0.80", help="set upper threshold value")


def main():
    args = parser.parse_args()
    seqfile = args.infile
    output = args.output
    upper_threshold = float(args.upperthreshold)

    seq = Predictor(seqfile, output)
    seq.change_upper_threshold(upper_threshold)
    pred = seq.predict()
    seq.predictions_to_csv(pred)
