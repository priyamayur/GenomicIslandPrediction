import sys
import pandas as pd
from treasureisland.Predictor import Predictor
from Bio import SeqIO


def flatten_result(pred):
    flat_result = []
    for org in pred.keys():
        org_result = pred[org]
        flat_result.extend(org_result)
    return flat_result


def main(seqfile):
    seq = Predictor(seqfile)
    #sequences = list(SeqIO.parse(seqfile, "fasta"))
    #for s in sequences:
    #    print(s.id)
    upper = 0.90
    lowers = [0.60, 0.70, 0.80, 0.90]
    for lower in lowers:
        seq.change_thresholds(upper,lower)
        pred = seq.predict()
        flat_pred = flatten_result(pred)
        df = pd.DataFrame(flat_pred, columns=['accession', 'start', 'end', 'probability'])
        filename = 'result_20_10000_' + upper + '-' + lower + '_revise_paper_gc.xlsx'
        pd.DataFrame(df).to_excel(filename)


if __name__ == '__main__':
    main(sys.argv[1])


