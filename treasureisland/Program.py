from treasureisland.IdentifyGI import IdentifyGI
import pickle
import pandas as pd
from Bio import SeqIO
import os
from . import models
from importlib import resources
import time


class Program:

    def __init__(self, input_file):
        self.input = input_file

    def format_input(self):
        file_extension = os.path.splitext(self.input)[1][1:]
        sequences = list(SeqIO.parse(self.input, file_extension))
        return sequences

    def process_output(self, output):
        all_gi_dict = {}
        for seq in output:
            for gi in seq.keys():
                gi_result = seq[gi]
                id = gi_result[0]
                start = gi_result[1] + 1
                end = gi_result[2]
                pred = gi_result[3]
                if id in all_gi_dict.keys():
                    all_gi_dict[id].append([id, start, end, pred])
                else :
                    all_gi_dict[id] = [[id, start, end, pred]]
        #print("output")
        #print(all_gi_dict)
        return all_gi_dict

    def get_predictions(self):
        start_time = time.time()
        print("--- %s seconds ---" % (time.time() - start_time))
        dna_sequence = self.format_input()
        dna_emb_model, classifier = self.get_models()
        genome = IdentifyGI(dna_sequence, dna_emb_model, classifier)
        fine_tuned_pred = genome.find_gi_predictions()

        output = self.process_output(fine_tuned_pred)
        print("--- %s seconds ---" % (time.time() - start_time))
        return output

    @staticmethod
    def get_models():
        read_classifier = resources.read_binary(models, "svm_upgrade_gensim_sklearn")
        classifier = pickle.loads(read_classifier)

        read_emb_model = resources.read_binary(models, "doc2vec_upgrade_gensim")
        dna_emb_model = pickle.loads(read_emb_model)

        return dna_emb_model, classifier

    def predictions_to_excel(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output_' + org + '.xlsx'
            pd.DataFrame(df).to_excel(filename)

    def predictions_to_csv(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output_' + org + '.csv'
            pd.DataFrame(df).to_csv(filename)


    def predictions_to_text(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output_' + org + '.txt'
            pd.DataFrame(df).to_csv(filename, header=None, index=None, sep=' ', mode='a')

