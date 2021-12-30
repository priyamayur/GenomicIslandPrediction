import sys
import pandas as pd
from treasureisland.Predictor import Predictor


def flatten_result(pred):
    flat_result = []
    for org in pred.keys():
        org_result = pred[org]
        flat_result.extend(org_result)
    return flat_result


def main(seqfile):
    seq = Predictor(seqfile)
    pred = seq.predict()
    flat_pred = flatten_result(pred)
    df = pd.DataFrame(flat_pred, columns=['accession', 'start', 'end', 'probability'])
    filename = 'result_104_10000_85_50_gc.xlsx'
    pd.DataFrame(df).to_excel(filename)


if __name__ == '__main__':
    main(sys.argv[1])


