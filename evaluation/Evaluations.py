import math as math
import pickle

class Evaluations:

    def getOverlap(self, a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

    def get_overlap_helper(self, true_list, predicted_list):
        tot_true = 0
        tot_overlap = 0
        for true in true_list:
            tot_true += true[2] - true[1]
            true_data = [true[1], true[2]]
            for pred in predicted_list:
                prediction_data = [pred[1], pred[2]]
                overlap = self.getOverlap(true_data, prediction_data)
                tot_overlap += overlap
        return tot_true, tot_overlap

    def calculate_score(self, predictor_dict, organism_pos_test_dict, organism_neg_dict, org_list):

        org_score = {}
        TTP = 0
        TFP = 0
        TTN = 0
        TFN = 0
        acc = 0
        prec = 0
        rec = 0
        f_score = 0
        mcc = 0
        for organism in org_list:
            pos_list = []
            neg_list = []
            predicted_list = []
            if organism in organism_pos_test_dict.keys():
                pos_list = organism_pos_test_dict[organism]
            if organism in organism_neg_dict.keys():
                neg_list = organism_neg_dict[organism]
            if organism in predictor_dict.keys():
                predicted_list = predictor_dict[organism]

            tot_pos, TP = self.get_overlap_helper(pos_list, predicted_list)
            FN = tot_pos - TP

            tot_neg, FP = self.get_overlap_helper(neg_list, predicted_list)
            TN = tot_neg - FP

            prec_org = TP / (TP + FP) if TP + FP != 0 else 1
            rec_org = TP / (TP + FN) if TP + FN != 0 else 0
            f_score_org = 2 * prec_org * rec_org / (prec_org + rec_org) if prec_org + rec_org != 0 else 0
            acc_org = (TP + TN) / (TP + TN + FP + FN) if TP + TN + FP + FN != 0 else 0
            mcc_org = ((TP * TN) - (FP * FN)) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) if math.sqrt(
                (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) != 0 else 0
            org_score[organism] = [{'TP': TP}, {'FP': FP}, {'TN': TN}, {'FN': FN}, {'Precision': prec_org},
                                   {'Recall': rec_org}, {'F-Score': f_score_org}, {'Accuracy': acc_org}, {'MCC': mcc_org}]

            TTP += TP
            TFP += FP
            TFN += FN
            TTN += TN

            acc += acc_org
            prec += prec_org
            rec += rec_org
            f_score += f_score_org
            mcc += mcc_org
        print(org_score)
        return [mcc, f_score, acc, prec, rec], [TTP, TFP, TFN, TTN], org_score

    def evaluations_test(self, total_orgs, models, gi_eval, organism_neg_dict):
        for key in models.keys():
            model = models[key]
            org_num = len(total_orgs)
            print("---------------")
            print(key)
            result1, result2 = self.calculate_score(model, gi_eval, organism_neg_dict, total_orgs)
            print("---------------")
            TTP = result2[0]
            TFP = result2[1]
            TFN = result2[2]
            TTN = result2[3]
            prec = TTP / (TTP + TFP) if TTP + TFP != 0 else 1
            rec = TTP / (TTP + TFN) if TTP + TFN != 0 else 0
            f_score = 2 * prec * rec / (prec + rec) if prec + rec != 0 else 0
            acc = (TTP + TTN) / (TTP + TTN + TFP + TFN) if TTP + TTN + TFP + TFN != 0 else 0
            mcc = ((TTP * TTN) - (TFP * TFN)) / math.sqrt(
                (TTP + TFP) * (TTP + TFN) * (TTN + TFP) * (TTN + TFN)) if math.sqrt(
                (TTP + TFP) * (TTP + TFN) * (TTN + TFP) * (TTN + TFN)) != 0 else 0
            print("mcc = ", mcc)
            print("f-score = ", f_score)
            print("accuracy = ", acc)
            print("precision = ", prec)
            print("recall = ", rec)

    def evaluations_main_104(self, total_orgs, models, gi_eval, organism_neg_dict, result_type):
        evaluation_result = []
        for key in models.keys():
            model = models[key]
            org_num = len(total_orgs)
            print("---------------")
            print(key)
            result1, result2, org_score = self.calculate_score(model, gi_eval, organism_neg_dict, total_orgs)
            print("---------------")
            print("mcc = ", result1[0] / org_num)
            print("f-score = ", result1[1] / org_num)
            print("accuracy = ", result1[2] / org_num)
            print("precision = ", result1[3] / org_num)
            print("recall = ", result1[4] / org_num)
            evaluation_result.append(org_score)
        self.result_download(evaluation_result, result_type)

    def result_download(self, org_score, result_type):
        with open("evaluation_result" + str(result_type),"wb") as fp:  # Pickling
            pickle.dump(org_score, fp)


    def mergeIntervals(self, arr):

        # Sorting based on the increasing order
        # of the start intervals
        arr.sort(key=lambda x: x[0])

        # array to hold the merged intervals
        m = []
        s = -10000
        max = -100000
        for i in range(len(arr)):
            a = arr[i]
            if a[0] > max:
                if i != 0:
                    m.append([s, max])
                max = a[1]
                s = a[0]
            else:
                if a[1] >= max:
                    max = a[1]

        # 'max' value gives the last point of
        # that particular interval
        # 's' gives the starting point of that interval
        # 'm' array contains the list of all merged intervals

        if max != -100000 and [s, max] not in m:
            m.append([s, max])
        # print("The Merged Intervals are :", end = " ")
        for i in range(len(m)):
            # print(m[i], end = " ")
            pass

        return m

    def calculate_novel(self,reference_list, predicted_list):
        merged_pred = self.mergeIntervals(reference_list)
        # print(predicted_list)
        novel = 0
        for pred in predicted_list:
            tot_overlap = 0
            pred_int = [pred[0], pred[1]]
            tot_pred = pred[1] - pred[0]
            for pos in merged_pred:
                pos_int = [pos[0], pos[1]]
                overlap = self.getOverlap(pred_int, pos_int)
                tot_overlap += overlap
            coverage = tot_overlap / tot_pred
            if (coverage == 0):
                novel += 1
        return novel

    def calculate_novel_104(self, reference_list, predictor, organism_list):
        total_novel = 0
        for org in organism_list:
            all_ref_pred = []
            all_pred = []

            if (org in predictor.keys()):
                pred_temp = predictor[org]
            else:
                pred_temp = []
            for pt in pred_temp:
                all_pred.append([pt[1], pt[2]])

            for ref in reference_list:
                if (org not in ref.keys()):
                    pos_list = []
                else:
                    pos_list = ref[org]
                for pos in pos_list:
                    all_ref_pred.append([pos[1], pos[2]])
            novel = self.calculate_novel(all_ref_pred, all_pred)
            total_novel += novel
        return total_novel

