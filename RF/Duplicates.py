def remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat, pairs_by_prot):
    unique_seqs = {}

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
            unique_seqs[freq_str] = id

        else:
            old_id = unique_seqs[freq_str]
            if pairs_by_prot[id] > pairs_by_prot[old_id]:
                unique_seqs[freq_str] = id

    unique_proteins = list(unique_seqs.values())

    return unique_proteins

def remove_ligands(ligand_counts, total_by_lig):
    unique_seqs = {}

    for lig in ligand_counts:
        freq_str = ""
        for kmer in ligand_counts[lig]:
            freq_str += str(ligand_counts[lig][kmer])
        if freq_str not in unique_seqs:
            unique_seqs[freq_str] = lig
        else:
            old_lig = unique_seqs[freq_str]
            if total_by_lig[lig] > total_by_lig[old_lig]:
                unique_seqs[freq_str] = lig

    unique_ligands = list(unique_seqs.values())

    return unique_ligands