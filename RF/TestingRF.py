import PreparingMatrix
import CombineLigandsProteins
import ReadingFasta

PreparingMatrix.export()
protmat = PreparingMatrix.proteins

CombineLigandsProteins.importdict()
logFC = CombineLigandsProteins.logFCdict
p = CombineLigandsProteins.pdict

ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logmat = []
pmat = []
for protein in proteins:
    logmat.append(logFC[protein.name])
    pmat.append(p[protein.name])
print(logmat)
print(pmat)