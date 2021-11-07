from treasureisland.IdentifyGI import IdentifyGI
import pickle
import pandas as pd
from Bio import SeqIO
import os
from . import models
from importlib import resources
import time


class Predictor:

    def __init__(self, file_path):
        self.input = file_path

    def __format_input(self, input):
        file_extension = os.path.splitext(input)[1][1:]
        sequences = list(SeqIO.parse(input, file_extension))
        return sequences

    def __process_output(self, output):
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, r'output')
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)
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
                else:
                    all_gi_dict[id] = [[id, start, end, pred]]
        return all_gi_dict


    def __get_models(self):
        read_classifier = resources.read_binary(models, "svm_upgrade_gensim_sklearn_1.0")
        classifier = pickle.loads(read_classifier)

        read_emb_model = resources.read_binary(models, "doc2vec_upgrade_gensim")
        dna_emb_model = pickle.loads(read_emb_model)

        return dna_emb_model, classifier


    def predict(self):
        '''
        :return: dictionary of genomic island predictions
        '''

        start_time = time.time()
        print("--- start predicting ---")
        dna_sequence = self.__format_input(self.input)
        dna_emb_model, classifier = self.__get_models()
        genome = IdentifyGI(dna_sequence, dna_emb_model, classifier)
        fine_tuned_pred = genome.find_gi_predictions()

        output = self.__process_output(fine_tuned_pred)
        print("--- finished predicting ---")
        print("--- %s seconds ---" % (time.time() - start_time))
        return output


    def predictions_to_excel(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output/output_' + org + '.xlsx'
            pd.DataFrame(df).to_excel(filename)


    def predictions_to_csv(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output/output_' + org + '.csv'
            pd.DataFrame(df).to_csv(filename)


    def predictions_to_text(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            filename = 'output/output_' + org + '.txt'
            pd.DataFrame(df).to_csv(filename, header=None, index=None, sep=' ', mode='a')
