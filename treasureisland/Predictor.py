from gensim.models import Doc2Vec

from treasureisland.IdentifyGI import IdentifyGI
import pickle
import pandas as pd
from Bio import SeqIO
import os
from . import models
from importlib_resources import files, as_file
from importlib import resources
from gensim.test.utils import get_tmpfile
import time
from treasureisland.Parameters import Parameters
import warnings
warnings.filterwarnings('ignore')

class Predictor:

    def __init__(self, input_file_path, output_file_path="output"):
        '''
        Initialize variables
        :param input_file_path: path to input sequence/s
        :param output_file_path: path to output files for GEIs
        '''
        self.input_file_path = input_file_path
        self.output_file_path = output_file_path
        self.parameters = Parameters()
        self.out_of_distribution = []

    def __format_input(self, input):
        '''
        Parse the input sequences
        :param input: path to input sequence/s
        :return: sequences : list of parsed fasta files
        '''
        sequences = list(SeqIO.parse(input, "fasta"))
        return sequences

    def change_upper_threshold(self, upper_threshold = 0.80):
        '''
        Sets the upper threshold to a user-defined value
        '''
        if (upper_threshold < 0.5):
            raise Exception("upper threshold cannot be lesser than 0.5")
        self.parameters.set_upper_threshold(upper_threshold)

    def __process_output(self, output):
        '''
        Prcoesses GEI output as a dictionary
        :param output: list of GEI output for each sequence
        :return: GEI output as a dictionary
        '''
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
        '''
        Get embedding and SVM models
        :return: embedding and SVM models
        '''
        read_classifier = resources.read_binary(models, "svm_embedding_v4")
        classifier = pickle.loads(read_classifier)

        source = files(models).joinpath('embedding_v4')
        dna_emb_model = Doc2Vec.load(get_tmpfile(source))

        return dna_emb_model, classifier


    def predict(self):
        '''
        Main function that predicts all input GEI sequences
        :return: dictionary of genomic island predictions
        '''

        start_time = time.time()
        print("--- start predicting ---")
        dna_sequence = self.__format_input(self.input_file_path)
        dna_emb_model, classifier = self.__get_models()

        genome = IdentifyGI(dna_sequence, dna_emb_model, classifier, self.parameters)
        fine_tuned_pred, self.out_of_distribution = genome.find_gi_predictions()

        output = self.__process_output(fine_tuned_pred)
        print("--- finished predicting ---")
        print("--- %s seconds ---" % (time.time() - start_time))
        return output


    def predictions_to_excel(self, predictions):
        '''
        Outputs prediciton as an excel file
        :param predictions: GEI prediction dictionary
        '''
        count = 0
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            if self.out_of_distribution[count]:
                org_name = org_name + '_out_of_distribution_data'
            filename = self.output_file_path + "/" + org_name + '.xlsx'
            pd.DataFrame(df).to_excel(filename)
            count += 1


    def predictions_to_csv(self, predictions):
        '''
        Outputs prediciton as an csv file
        :param predictions: GEI prediction dictionary
        '''
        count = 0
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            if self.out_of_distribution[count]:
                org_name = org_name + '_out_of_distribution_data'
            filename = self.output_file_path + "/" + str(org_name) + '.csv'
            pd.DataFrame(df).to_csv(filename)
            count += 1


    def predictions_to_text(self, predictions):
        '''
        Outputs prediciton as an text file
        :param predictions: GEI prediction dictionary
        '''
        count = 0
        for org in predictions.keys():
            df = pd.DataFrame(predictions[org], columns=['accession', 'start', 'end', 'probability'])
            org_name = ''.join(e for e in str(org) if e.isalnum())
            if self.out_of_distribution[count]:
                org_name = org_name + '_out_of_distribution_data'
            filename = self.output_file_path + "/" + org_name + '.txt'
            pd.DataFrame(df).to_csv(filename, header=None, index=None, sep=' ', mode='w')
            count += 1
