from Bio import SeqIO

test_data = list(SeqIO.parse("C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/get_evaluation/main_organisms_104.fasta","fasta"))

test = test_data[0:30]

SeqIO.write(test, "C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/get_evaluation/test_104.fasta", "fasta")