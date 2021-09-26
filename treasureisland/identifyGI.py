import pandas as pd
import numpy as np
from Bio.Seq import Seq

class identifyGI:

    def __init__(self,dna_sequence_list, window_size,kmer_size, dna_emb_model, classifier,upper_threshold,lower_threshold,tune_metric,minimum_gi_size):
        '''Initialize variables'''

        self.dna_sequence_list = dna_sequence_list
        self.window_size = window_size
        self.kmer_size = kmer_size
        self.dna_emb_model = dna_emb_model
        self.classifier = classifier
        self.upper_threshold = upper_threshold
        self.lower_threshold = lower_threshold
        self.tune_metric = tune_metric
        self.minimum_gi_size = minimum_gi_size
        

    def generate_kmers(self,segment):
        '''Generate overlapping kmers from a given DNA sequence
           segment : DNA sequence
        '''
        kmers = []
        for i in range(0, len(segment) - (self.kmer_size-1)):

            start = i
            end = i + self.kmer_size
            kmer = segment[start:end]
            kmers.append(kmer)

        return kmers       

    def split_dna_sequence(self, dna_sequence):
        '''
        Divides the DNA segment into equal small segments of sizes self.window_size 
        '''
        
        sequence = str(dna_sequence.seq).lower()
        processed_dna_seq = []
        segment_borders = []

        for i in range(0, len(sequence),self.window_size ):

            start = i
            end = i + self.window_size
            segment = sequence[start:end]
            kmers = self.generate_kmers(segment)
            processed_dna_seq.append(kmers)
            segment_borders.append([start,end])    

        return processed_dna_seq, segment_borders    

    def get_dna_vectors(self,processed_dna_seq):
        '''Uses the DNA embedding model to get the DNA vector of a DNA segment
           processed_dna_seq : DNA segment after preprocessing step
        '''
        
        dna_vectors = []
        for segments in processed_dna_seq:

            inferred_vector = self.dna_emb_model.infer_vector(segments,epochs=20)
            dna_vectors.append(inferred_vector)
            
        return dna_vectors    

    def get_dna_segment_probability(self, dna_vectors, segment_borders) :
        ''' 
        Gets the probability of each DNA segment to be a GI using the classifier
        dna_vectors : DNA vectors 
        segment_borders : start and end points of DNA segments  
        '''

        probability = self.classifier.predict_proba(dna_vectors)
        prob_df = pd.DataFrame(np.column_stack([segment_borders,probability]), columns = ['start', 'end','0','1'])
        prob_list = prob_df.values.tolist()

        return prob_list


    def get_GI_regions(self,dna_prob):
        ''' 
        Gets the DNA segments with probability above the upper threshold and 
        flanking sequences with probability between upper and lower thesholds 
        dna_prob : full list of DNA segments and their probabilities
        '''
        prev = -1
        gi_dict ={}
        gi_num = 1  
        prev_gi = False
        for row in dna_prob:
            found_gi = False
            if (row[3] >= self.upper_threshold):
                found_gi = True
                if 'gi_'+ str(gi_num) in gi_dict.keys():
                    gi_dict['gi_'+ str(gi_num)].append(row)
                else:
                    if (prev == -1) :      
                        gi_dict['gi_'+ str(gi_num)] = [row]
                    else:   
                        if prev[3] > self.lower_threshold :        
                            gi_dict['gi_'+ str(gi_num)] = [prev]
                            gi_dict['gi_'+ str(gi_num)].append(row)    
                        else :
                            gi_dict['gi_'+ str(gi_num)] = [row]     
                prev_gi = True         
            if found_gi == False and prev_gi == True:
                if row[3] > self.lower_threshold :
                    gi_dict['gi_'+ str(gi_num)].append(row)
                prev_gi = False
                gi_num += 1
            prev = row
          
        return gi_dict  

    def find_fragment_probability(self, GI_borders, dna_sequence):
        '''
        Get a DNA fragment probability for class GI
        GI_borders : set of start and end points
        '''
        gi_start = int(GI_borders[0])
        gi_end = int(GI_borders[1])  
        
        sequence = str(dna_sequence.seq).lower()        
        fragment = sequence[gi_start:gi_end]
        kmers = self.generate_kmers(fragment) 
        
        inferred_vector = [self.dna_emb_model.infer_vector(kmers, epochs=20)]
        prob = self.classifier.predict_proba(inferred_vector) 
        
        gi_prob = prob[0][1]  # gi_prob probability of region belonging to class 1(GI)

        return gi_prob

    def pre_fine_tune(self, gi_regions, gi, dna_sequence):
        '''
           Determines the start and end points and upper and lower limits for each GI fragment, 
           before the fine tuning step begins.
           Also, gets the merged GI fragments and their probabilities
           gi_regions : GI fragments along with its flanking sequence
           gi : a particular GI key
        '''
        first_frag_gi = 0

        # Determine the start, end, upper and lower limit for each GI

        if (gi_regions[gi][0][3] < self.upper_threshold ):     
            start = gi_regions[gi][1][0]
            lower_limit = (gi_regions[gi][0][0]) + (self.window_size/2)
        else :     
            start = gi_regions[gi][0][0]
            lower_limit = (gi_regions[gi][0][0]) 
            first_frag_gi = 1

        if (gi_regions[gi][-1][3] < self.upper_threshold):
            end = gi_regions[gi][-2][1]
            upper_limit = (gi_regions[gi][-1][1]) - (self.window_size/2)
        else:
            end = gi_regions[gi][-1][1]
            upper_limit = (gi_regions[gi][-1][1])   

        # merged_GI : GI segments after merging(if any)       

        merged_GI_borders = [start,end]

        # Determine probability of merged regions

        if end-start <= (self.window_size) :
            if first_frag_gi == 0:       
                merged_GI_prob = gi_regions[gi][1][3]
            if first_frag_gi == 1:       
                merged_GI_prob = gi_regions[gi][0][3]  
                first_frag_gi = 0 
        else:
            merged_GI_prob = self.find_fragment_probability(merged_GI_borders, dna_sequence)     

        tune_limit = [lower_limit,upper_limit]  # upper and lower limits for tuning borders for each GI

        return merged_GI_borders, tune_limit, merged_GI_prob
            

    def get_tuned_borders(self,left_border,right_border,tune_limit,left_tune_metric,right_tune_metric,tune_probs, dna_sequence):
        '''
        Fine tuning method to fine tune each GI fragment to return the list of possible new borders and their probability
        left_border : current GI left border
        right_border : current GI right border
        tune_limit : contains the left and right limits till which fine tune can be done
        left_tune_metric : length of fragment by which current left border shifts in the loop
        right_tune_metric : length of fragment by which current right border shifts in the loop
        tune_probs : list of all new GI fragments and their probabilities
        '''
        left_limit = tune_limit[0]
        right_limit =  tune_limit[1] 
        old_borders = [left_border,right_border] 
        new_borders = [left_border,right_border]  

        travelling_inward = True if left_tune_metric > 0 or right_tune_metric < 0 else False
        tune_metric = abs(left_tune_metric) if left_tune_metric != 0 else abs(right_tune_metric)
        border_left = True if right_tune_metric == 0  else False

        for i in range(int((self.minimum_gi_size/(2*tune_metric)))) :
            left_border += left_tune_metric
            right_border += right_tune_metric  
            if ((right_border - left_border + 1 ) < self.minimum_gi_size or left_border < (left_limit) or right_border > right_limit):      
                break
            frag_border = [left_border,right_border]
            frag_prob = self.find_fragment_probability(frag_border, dna_sequence)
            if (travelling_inward):
                previous_prob = tune_probs[-1][2]       
                if (frag_prob < self.lower_threshold):                    
                    if (border_left):
                        new_borders = [old_borders[0],tune_probs[-1][0]]
                    else: 
                        new_borders = [tune_probs[-1][1],old_borders[1]]
                    break  
                if (previous_prob > frag_prob):
                    break        
            else:
                if (frag_prob < self.upper_threshold):
                    break    
            tune_probs.append([left_border,right_border,frag_prob]) 
            new_borders = [left_border,right_border] 

        return tune_probs, new_borders         

    def fine_tune_borders(self, merged_GI_borders,tune_limit,merged_GI_prob, dna_sequence):
        '''
        Gets the fine-tuned GI borders and determines the best border for each GI
        merged_GI_borders : GI border after merging
        tune_limit : left and right limits till which it can be tuned
        merged_GI_prob : probability of the merged GI segment
        '''
        # each border and together -increase, decrease 
        tune_probs = []  
          
        left_border = merged_GI_borders[0]
        right_border = merged_GI_borders[1]  
        new_border_prob = (left_border,right_border,merged_GI_prob)

        #Final borders of each GI after the fine tuning        

        final_left = left_border     
        final_right = right_border

        tune_probs.append([left_border, right_border, merged_GI_prob])

        #Find left outer boundary
        
        tune_probs, new_limits = self.get_tuned_borders(left_border,right_border,tune_limit,-self.tune_metric,0,tune_probs, dna_sequence) 
        left_ob = new_limits[0]
        
        #Find left inner boundary
        
        tune_probs, new_limits = self.get_tuned_borders(left_border,right_border,tune_limit,self.tune_metric,0,tune_probs, dna_sequence)
        left_ib = new_limits[0]         
        
        #Find right outer boundary
        
        tune_probs, new_limits = self.get_tuned_borders(left_border,right_border,tune_limit,0,self.tune_metric,tune_probs, dna_sequence)
        right_ob =  new_limits[1]
        
        #Find right inner boundary
        
        tune_probs, new_limits = self.get_tuned_borders(left_border,right_border,tune_limit,0,-self.tune_metric,tune_probs, dna_sequence)
        right_ib = new_limits[1]
        
        #If left and/or right borders have moved in, then consider inner border, 
        #if left and/or right borders have moved out, then consider outer border,
        #if both, then give preference to outer border(i.e. increase recall over precision)
        
        if (left_ib != left_border):
            final_left = left_ib
        if (left_ob != left_border):
            final_left =  left_ob  

        if (right_ib != right_border):
            final_right = right_ib
        if (right_ob != right_border):
            final_right =  right_ob          

        
        frag_prob_final = self.find_fragment_probability([final_left,final_right], dna_sequence)        
       
        new_border_prob = (final_left,final_right,frag_prob_final)
        
        return new_border_prob

    def find_GI_borders(self,gi_regions,id, dna_sequence):
        ''' 
        Function to combine all the GI border finding steps
        gi_regions : GI fragments along with its flanking sequence
        '''
          
        gi_borders = {} 
          
        for gi in gi_regions.keys():

            merged_GI_borders, tune_limit, merged_GI_prob = self.pre_fine_tune(gi_regions, gi, dna_sequence)
            new_border_prob = self.fine_tune_borders(merged_GI_borders,tune_limit,merged_GI_prob, dna_sequence)
            gi_borders[gi] = [id, new_border_prob[0], new_border_prob[1], new_border_prob[2]]

        return gi_borders    


    def find_gi_predictions(self):
        ''' 
        Main function to call all other functions for identifying GI regions
        '''   
        all_gi_borders =[]
        for dna_sequence in self.dna_sequence_list:
            id = dna_sequence.id
            processed_dna_seq,segment_borders = self.split_dna_sequence(dna_sequence)
            dna_vectors = self.get_dna_vectors(processed_dna_seq)
            dna_prob = self.get_dna_segment_probability(dna_vectors, segment_borders)
            gi_regions = self.get_GI_regions(dna_prob) 
            gi_borders = self.find_GI_borders(gi_regions, id, dna_sequence)
            all_gi_borders.append(gi_borders)
            
        return all_gi_borders 