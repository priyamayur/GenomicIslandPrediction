#!/usr/bin/env python
import sys
import argparse
from treasureisland.Predictor import Predictor


parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('r'))
parser.add_argument("-o", "--output", required=False, default="output", help="output directory name")


def main():
    args = parser.parse_args()
    seqfile = args.infile
    output = args.output

    seq = Predictor(seqfile, output)
    pred = seq.predict()
    seq.predictions_to_text(pred)
