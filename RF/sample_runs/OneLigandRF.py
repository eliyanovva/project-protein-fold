import PreparingMatrix
import CombineLigandsProteins
import ReadingFasta
import RandomForest

PreparingMatrix.export()
protmat = PreparingMatrix.proteins

CombineLigandsProteins.exportdicts()
logFC = CombineLigandsProteins.citlog
p = CombineLigandsProteins.citp
cor = CombineLigandsProteins.citcor

ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logmat = []
pmat = []
for protein in proteins:
    logmat.append(cor[protein.name])
    pmat.append(p[protein.name])

RandomForest.train(protmat, logmat)
RandomForest.train(protmat, pmat)
