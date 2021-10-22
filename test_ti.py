#!/usr/bin/env python
import sys
from treasureisland.dna_sequence import sequence


def main(seqfile):
    seq = sequence(seqfile)
    pred = seq.predict()
    seq.predictions_to_csv(pred)


if __name__ == '__main__':
    main(sys.argv[1])

