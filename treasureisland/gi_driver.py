from .identifyGI import identifyGI
import pickle 
import pandas as pd
from Bio import SeqIO
from gensim.test.utils import get_tmpfile
from gensim.models.doc2vec import Doc2Vec

class gi_driver:

  def __init__(self, input_file_path):
    self.input = input_file_path  


  def format_input(self):
    input_format = self.input.split(".")[1]    
    sequences = list(SeqIO.parse(self.input, input_format))
    return sequences

  def process_output(self, output):
        gi_list = []
        for gi in output.keys():
            gi_result = output[gi]
            start = gi_result[0] +1 
            end = gi_result[1]
            pred = gi_result[2]
            gi_list.append([start, end, pred])
        pred_df = pd.DataFrame(gi_list, columns=['start', 'end', 'probability'])
        return pred_df  

  def get_predictions(self):
    #Parameters 
    window_size = 10000
    kmer_size = 6
    upper_threshold = 0.75
    lower_threshold = 0.50
    tune_metric = 1000
    minimum_gi_size = 10000

    fname = get_tmpfile("C:/Users/USER/GenomicIslandPrediction/treasureisland/models/doc2vec_dbow_50")   
    dna_emb_model = Doc2Vec.load(fname)
    with open('C:/Users/USER/GenomicIslandPrediction/treasureisland/models/' + 'svm_classifier' , 'rb') as handle:
          classifier = pickle.load(handle)
    dna_sequence = self.format_input()

    genome = identifyGI(dna_sequence, window_size,kmer_size, dna_emb_model, classifier,upper_threshold,lower_threshold,tune_metric,minimum_gi_size)
    fine_tuned_pred = genome.find_gi_predictions()
    
    output_dataframe = self.process_output(fine_tuned_pred)    

    print(fine_tuned_pred)

    return output_dataframe

  def predictions_to_excel(self,predictions):
    return pd.DataFrame(predictions).to_excel('output.xlsx')

  def predictions_to_csv(self,predictions):
    return pd.DataFrame(predictions).to_csv('output.csv')

  def predictions_to_text(self,predictions):
    return pd.DataFrame(predictions).to_csv('output.txt', header=None, index=None, sep=' ', mode='a')



    


