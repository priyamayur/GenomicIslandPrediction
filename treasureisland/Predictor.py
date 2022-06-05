from treasureisland.IdentifyGI import IdentifyGI
import pickle
import pandas as pd
from Bio import SeqIO
import os
from . import models
from importlib import resources
import time
from treasureisland.Parameters import Parameters

class Predictor:

    def __init__(self, input_file_path, output_file_path="output"):
        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.parameters = Parameters()

    def __format_input(self, input):
        sequences = list(SeqIO.parse(input, "fasta"))
        return sequences

    def change_upper_threshold(self, upper_threshold = 0.8):
        if (upper_threshold < 0.5):
            raise Exception("upper threshold cannot be lesser than 0.5")
        self.parameters.set_upper_threshold(upper_threshold)

    def __process_output(self, output):
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, self.output_file_path)
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
        read_classifier = resources.read_binary(models, "svm_genome_centric")
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
        dna_sequence = self.__format_input(self.input_file_path)
        dna_emb_model, classifier = self.__get_models()

        genome = IdentifyGI(dna_sequence, dna_emb_model, classifier, self.parameters)
        fine_tuned_pred = genome.find_gi_predictions()

        output = self.__process_output(fine_tuned_pred)
        print("--- finished predicting ---")
        print("--- %s seconds ---" % (time.time() - start_time))
        return output


    def predictions_to_excel(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            filename = self.output_file_path + "/" + org_name + '.xlsx'
            pd.DataFrame(df).to_excel(filename)


    def predictions_to_csv(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            filename = self.output_file_path + "/" + str(org_name) + '.csv'
            pd.DataFrame(df).to_csv(filename)


    def predictions_to_text(self, predictions):
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            filename = self.output_file_path + "/" + org_name + '.txt'
            pd.DataFrame(df).to_csv(filename, header=None, index=None, sep=' ', mode='w')
