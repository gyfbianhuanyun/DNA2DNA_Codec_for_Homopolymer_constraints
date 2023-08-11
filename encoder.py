import numpy as np


# DNA to DNA sequence encoder
def homopolymer_coding(dna_seq, homopolymer):
    # dna_seq: list, homopolymer: int
    homo_count = 1
    dna_base = ''
    dna_homo_encoded_seq = []
    i = 0
    while i < len(dna_seq):
        # check homopolymer value
        if homo_count != homopolymer:
            if dna_seq[i] == dna_base:
                homo_count += 1
            else:
                dna_base = dna_seq[i]
                homo_count = 1
            dna_homo_encoded_seq.append(dna_seq[i])
            # upgrade i
            i += 1

        # homopolymer coding
        else:
            # if the number of base is 1 after base with homopolymer
            if i + 1 == len(dna_seq):
                if dna_base == dna_seq[i]:
                    if dna_base == 'A':
                        dna_homo_encoded_seq.append('C')
                    elif dna_base == 'C':
                        dna_homo_encoded_seq.append('A')
                    elif dna_base == 'G':
                        dna_homo_encoded_seq.append('T')
                    else:
                        dna_homo_encoded_seq.append('G')
                dna_homo_encoded_seq.append(dna_seq[i])
                i += 1
            else:
                if dna_base == 'A':
                    if dna_seq[i] == 'A':
                        dna_homo_encoded_seq.append('C')
                    elif dna_seq[i] == 'C' and dna_seq[i+1] == 'A':
                        dna_homo_encoded_seq.append('G')
                    dna_homo_encoded_seq.append(dna_seq[i])
                    dna_homo_encoded_seq.append(dna_seq[i + 1])
                    if dna_seq[i] == 'G' and dna_seq[i+1] == 'C':
                        dna_homo_encoded_seq.append('G')
                elif dna_base == 'C':
                    if dna_seq[i] == 'C':
                        dna_homo_encoded_seq.append('A')
                    elif dna_seq[i] == 'A' and dna_seq[i+1] == 'C':
                        dna_homo_encoded_seq.append('T')
                    dna_homo_encoded_seq.append(dna_seq[i])
                    dna_homo_encoded_seq.append(dna_seq[i + 1])
                    if dna_seq[i] == 'T' and dna_seq[i+1] == 'A':
                        dna_homo_encoded_seq.append('T')
                elif dna_base == 'G':
                    if dna_seq[i] == 'G':
                        dna_homo_encoded_seq.append('T')
                    elif dna_seq[i] == 'T' and dna_seq[i+1] == 'G':
                        dna_homo_encoded_seq.append('A')
                    dna_homo_encoded_seq.append(dna_seq[i])
                    dna_homo_encoded_seq.append(dna_seq[i + 1])
                    if dna_seq[i] == 'A' and dna_seq[i+1] == 'T':
                        dna_homo_encoded_seq.append('A')
                else:
                    if dna_seq[i] == 'T':
                        dna_homo_encoded_seq.append('G')
                    elif dna_seq[i] == 'G' and dna_seq[i+1] == 'T':
                        dna_homo_encoded_seq.append('C')
                    dna_homo_encoded_seq.append(dna_seq[i])
                    dna_homo_encoded_seq.append(dna_seq[i + 1])
                    if dna_seq[i] == 'C' and dna_seq[i+1] == 'G':
                        dna_homo_encoded_seq.append('C')

                # upgrade homo_count and dna_base
                if dna_seq[i] != dna_seq[i+1]:
                    homo_count = 1
                else:
                    homo_count = 2
                dna_base = dna_homo_encoded_seq[-1]

                # upgrade i
                i += 2

    return dna_homo_encoded_seq
