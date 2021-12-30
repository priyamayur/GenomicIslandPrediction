from Bio import SeqIO

test_data = list(SeqIO.parse("C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/get_evaluation/main_organisms_104.fasta","fasta"))
test_org = ['NC_004116.1',
 'NC_009504.1',
 'NC_004603.1',
 'NC_010322.1',
 'NC_008024.1',
 'NC_010501.1',
 'NC_009708.1',
 'NC_008253.1',
 'NC_008563.1',
 'NC_008321.1',
 'NC_010515.1',
 'NC_007606.1',
 'NC_010334.1',
 'NC_010473.1',
 'NC_007958.1',
 'NC_009342.1',
 'NC_009512.1',
 'NC_009783.1',
 'NC_009665.1',
 'NC_005071.1']

test = []
for t in test_data:
    t_id = t.id
    t_org = t_id.split("_")[0] + "_" + t_id.split("_")[1]
    if (t_org in test_org):
        test.append(t)
print(len(test))
SeqIO.write(test, "C:/Users/kkp834/PycharmProjects/GenomicIslandPrediction/get_evaluation/test_104.fasta", "fasta")