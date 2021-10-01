from treasureisland.identifyGI import identifyGI
import pickle 
import pandas as pd
from Bio import SeqIO
from gensim.test.utils import get_tmpfile
from gensim.models.doc2vec import Doc2Vec
import os
import pkgutil
from . import models
from importlib import resources


class gi_driver:

  def __init__(self, input_file_path):
    self.input = input_file_path  


  def format_input(self):
    sequences = list(SeqIO.parse(self.input, "fasta"))
    return sequences

  def process_output(self, output):
    gi_list = []
    for seq in output:
      for gi in seq.keys():
        gi_result = seq[gi]
        id = gi_result[0]
        start = gi_result[1] + 1 
        end = gi_result[2]
        pred = gi_result[3]
        gi_list.append([id,start, end, pred])
        
    pred_df = pd.DataFrame(gi_list, columns=['id','start', 'end', 'probability'])
    return pred_df  

  def get_predictions(self):
    #Parameters 
    window_size = 10000
    kmer_size = 6
    upper_threshold = 0.75
    lower_threshold = 0.50
    tune_metric = 1000
    minimum_gi_size = 10000
    
    read_classifier = resources.read_binary(models, "svm_upgrade_gensim_sklearn")
    classifier = pickle.loads(read_classifier)
    
    read_emb_model = resources.read_binary(models, "doc2vec_upgrade_gensim")
    dna_emb_model = pickle.loads(read_emb_model)
    
    dna_sequence = self.format_input()

    genome = identifyGI(dna_sequence, window_size,kmer_size, dna_emb_model, classifier,upper_threshold,lower_threshold,tune_metric,minimum_gi_size)
    fine_tuned_pred = genome.find_gi_predictions()
    
    output_dataframe = self.process_output(fine_tuned_pred)    

    #print(fine_tuned_pred)

    return output_dataframe

  def predictions_to_excel(self,predictions):
    return pd.DataFrame(predictions).to_excel('output.xlsx')

  def predictions_to_csv(self,predictions):
    return pd.DataFrame(predictions).to_csv('output.csv')

  def predictions_to_text(self,predictions):
    return pd.DataFrame(predictions).to_csv('output.txt', header=None, index=None, sep=' ', mode='a') 


