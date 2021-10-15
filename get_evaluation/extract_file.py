from Bio import SeqIO

test_data = list(SeqIO.parse("C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/genome/main_organisms_104.fasta","fasta"))

test = test_data[0:2]

SeqIO.write(test, "C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/genome/test_104.fasta", "fasta")