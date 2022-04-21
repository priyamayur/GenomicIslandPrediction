from evaluation.Evaluations import Evaluations
from importlib import resources
from . import reference_data
from . import predictions
import pandas as pd


class Program:

    def __init__(self, input_file):
        self.input = input_file

    def get_organism_info(self, data_set):
        organism = {}
        for i, data in data_set.iterrows():
            acc = data.Accession
            if acc in organism.keys():
                organism[acc].append([acc, data.Start, data.End])
            else:
                organism[acc] = [[acc, data.Start, data.End]]
        return organism

    def get_organism_info_predictors(self, data_set):
        organism = {}
        for i, data in data_set.iterrows():
            acc = data.accession
            if acc in organism.keys():
                organism[acc].append([acc, data.start, data.end])
            else:
                organism[acc] = [[acc, data.start, data.end]]
        return organism

    def get_dictionary(self):
        literature_pos_data_table, neg_data_table, pos_data_table, neg_test_data_table, pos_test_data_table = self.get_reference_data()
        alien_hunter, islander, islandpath_dimob, islandpick, islandviewer, sigi_hmm = self.get_prediciton_data()

        organism_pos_dict = self.get_organism_info(pos_data_table)
        organism_neg_dict = self.get_organism_info(neg_data_table)
        organism_pos_test_dict = self.get_organism_info(pos_test_data_table)
        organism_lit_dict = self.get_organism_info(literature_pos_data_table)
        organism_neg_test_dict = self.get_organism_info(neg_test_data_table)

        alien_hunter_dict = self.get_organism_info_predictors(alien_hunter)
        islandviewer_dict = self.get_organism_info_predictors(islandviewer)
        islandpath_dimob_dict = self.get_organism_info_predictors(islandpath_dimob)
        sigi_hmm_dict = self.get_organism_info_predictors(sigi_hmm)
        islandpick_dict = self.get_organism_info_predictors(islandpick)
        islander_dict = self.get_organism_info_predictors(islander)
        treasureisland = pd.read_excel(self.input)
        treasureisland_dict = self.get_organism_info_predictors(treasureisland)

        model_dict = {'islandviewer_dict': islandviewer_dict, 'islandpath_dimob_dict': islandpath_dimob_dict,
                      'sigi_hmm_dict': sigi_hmm_dict,
                      'islandpick_dict': islandpick_dict, 'islander_dict': islander_dict,
                      'alien_hunter_dict': alien_hunter_dict,
                      'treasureisland_dict': treasureisland_dict}

        return organism_pos_dict, organism_neg_dict, organism_pos_test_dict, organism_neg_test_dict, organism_lit_dict, model_dict

    def get_evaluation_result(self):
        eval = Evaluations()
        organism_pos_dict, organism_neg_dict, organism_pos_test_dict, organism_neg_test_dict, organism_lit_dict, model_dict = self.get_dictionary()

        print("evaluation of complete 104 organism")
        print(model_dict['treasureisland_dict'].keys())
        eval.evaluations_main_104(model_dict['treasureisland_dict'].keys(), model_dict, organism_pos_dict, organism_neg_dict, "complete")
        #eval.evaluations_main_104(organism_pos_dict.keys(), model_dict, organism_pos_dict, organism_neg_dict)
        # ---------------------------------------------------------------------------------------
        all_predictor = [model_dict['islandpath_dimob_dict'], model_dict['sigi_hmm_dict'],
                         model_dict['islandpick_dict'], model_dict['islander_dict'], model_dict['alien_hunter_dict'],
                         model_dict['treasureisland_dict']]
        names = ['islandpath_dimob_dict', 'sigi_hmm_dict', 'islandpick_dict', 'islander_dict', 'alien_hunter_dict',
                 'treasureisland_dict']
        # ---------------------------------------------------------------------------------------
        print("evaluation of novel GEIs for complete 104 organism")
        for all in range(len(all_predictor)):
            predictor = all_predictor[all]
            reference_list = all_predictor[0:all] + all_predictor[all + 1:]
            print(names[all])
            print(names[0:all] + names[all + 1:])
            organism_list = model_dict['treasureisland_dict'].keys()
            print(eval.calculate_novel_104(reference_list, predictor, organism_list))
        # ---------------------------------------------------------------------------------------
        print("evaluation of literature data")
        eval.evaluations_main_104(organism_lit_dict.keys(), model_dict, organism_lit_dict, organism_neg_dict, "literature")
        # ---------------------------------------------------------------------------------------
        print("evaluation of test data")
        pos_orgs = set(organism_pos_test_dict.keys())
        print(len(pos_orgs))
        neg_orgs = set(organism_neg_test_dict.keys())
        total_orgs = pos_orgs.union(neg_orgs)
        eval.evaluations_main_104(pos_orgs, model_dict, organism_pos_test_dict, organism_neg_test_dict, "test")
        #eval.evaluations_test(total_orgs, model_dict, organism_pos_test_dict, organism_neg_test_dict)


    @staticmethod
    def get_reference_data():
        read_literature = resources.read_binary(reference_data, "GI_literature_set_table.xlsx")
        literature_pos_data_table = pd.read_excel(read_literature)
        read_negative = resources.read_binary(reference_data, "GI_negative_set_table.xlsx")
        neg_data_table = pd.read_excel(read_negative)
        read_positive = resources.read_binary(reference_data, "GI_positive_set_table.xlsx")
        pos_data_table = pd.read_excel(read_positive)
        #read_negative_test = resources.read_binary(reference_data, "negative_test_table.xlsx")
        read_negative_test = resources.read_binary(reference_data, "negative_test_table_gc.xlsx")
        neg_test_data_table = pd.read_excel(read_negative_test)
        #read_positive_test = resources.read_binary(reference_data, "positive_test_table.xlsx")
        read_positive_test = resources.read_binary(reference_data, "positive_test_table_gc.xlsx")
        pos_test_data_table = pd.read_excel(read_positive_test)

        return literature_pos_data_table, neg_data_table, pos_data_table, neg_test_data_table, pos_test_data_table

    @staticmethod
    def get_prediciton_data():
        read_ah = resources.read_binary(predictions, "alien_hunter.xlsx")
        alien_hunter = pd.read_excel(read_ah)
        read_is = resources.read_binary(predictions, "islander.xlsx")
        islander = pd.read_excel(read_is)
        read_ipd = resources.read_binary(predictions, "islandpath_dimob.xlsx")
        islandpath_dimob = pd.read_excel(read_ipd)
        read_ip = resources.read_binary(predictions, "islandpick.xlsx")
        islandpick = pd.read_excel(read_ip)
        read_iv = resources.read_binary(predictions, "islandviewer.xlsx")
        islandviewer = pd.read_excel(read_iv)
        read_sh = resources.read_binary(predictions, "sigi_hmm.xlsx")
        sigi_hmm = pd.read_excel(read_sh)
        # read_ti = resources.read_binary(predictions, "my_model_predictions_2.xlsx")
        # treasureisland = pd.read_excel(read_ti)

        return alien_hunter, islander, islandpath_dimob, islandpick, islandviewer, sigi_hmm
