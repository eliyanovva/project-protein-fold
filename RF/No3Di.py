import PreparingMatrix
import CombineLigandsProteins
import ReadingFasta
import RandomForest

ReadingFasta.import_variables()
protmat = ReadingFasta.sequence_matrix

CombineLigandsProteins.exportdicts()
logFC = CombineLigandsProteins.citlog
p = CombineLigandsProteins.citp

ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logmat = []
pmat = []
for protein in proteins:
    logmat.append(logFC[protein.name])
    pmat.append(p[protein.name])

RandomForest.train(protmat, logmat)
RandomForest.train(protmat, pmat)