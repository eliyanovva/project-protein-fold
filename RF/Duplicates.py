def remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat):
    unique_seqs = set()
    unique_proteins = set()

    for i in range(4):
        for id in AA_seqvar[i]:
            for kmer in AA_feat[i]:
                if kmer not in AA_seqvar[i][id]:
                    AA_seqvar[i][id][kmer] = 0
    for i in range(4):
        for id in Di_seqvar[i]:
            for kmer in Di_feat[i]:
                if kmer not in Di_seqvar[i][id]:
                    Di_seqvar[i][id][kmer] = 0

    all_ids = []
    for id in AA_seqvar[0]:
        all_ids.append(id)

    for k in range(len(all_ids)):
        id = all_ids[k]
        freq_str = ""
        for i in range(4):
            for kmer in AA_feat[i]:
                freq_str += str(AA_seqvar[i][id][kmer])
        for i in range(4):
            for kmer in Di_feat[i]:
                freq_str += str(Di_seqvar[i][id][kmer])
        if freq_str not in unique_seqs:
            unique_seqs.add(freq_str)
            unique_proteins.add(id)
    return unique_proteins

def remove_ligands(ligand_counts):
    unique_seqs = set()
    unique_ligands = []

    for lig in ligand_counts:
        freq_str = ""
        for kmer in ligand_counts[lig]:
            freq_str += str(ligand_counts[lig][kmer])
        if freq_str not in unique_seqs:
            unique_seqs.add(freq_str)
            unique_ligands.append(lig)

    return unique_ligands