from Bio import SeqIO
import os
res = []
for x in next(os.walk('.'))[1]:

    with open(x + '/contigs.fasta') as f:
        for r in SeqIO.parse(f, "fasta"):
            res.append((x[3:], len(r.seq)))
res.sort(key=lambda x:(int(x[0])))
for i in res: print(i[1])