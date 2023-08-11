import numpy as np


# homopolymer coding decoder
def homo_decoder(dna_seq, homopolymer):
    # dna_seq: list, homopolymer: int
    homo_count = 1
    dna_base = ''
    dna_homo_decoded_seq = []
    i = 0
    while i < len(dna_seq):
        # check homopolymer value
        if homo_count != homopolymer:
            if dna_seq[i] == dna_base:
                homo_count += 1
            else:
                dna_base = dna_seq[i]
                homo_count = 1
            dna_homo_decoded_seq.append(dna_seq[i])
            # upgrade i
            i += 1

        # homopolymer coding
        else:
            if i + 1 == len(dna_seq):
                dna_homo_decoded_seq.append(dna_seq[i])
                i += 1
            elif i + 2 == len(dna_seq):
                if dna_base == 'A':
                    if dna_seq[i] == 'C' and dna_seq[i + 1] == 'A':
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                elif dna_base == 'C':
                    if dna_seq[i] == 'A' and dna_seq[i + 1] == 'C':
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                elif dna_base == 'G':
                    if dna_seq[i] == 'T' and dna_seq[i + 1] == 'G':
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                else:
                    if dna_seq[i] == 'G' and dna_seq[i + 1] == 'T':
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                i += 2
            else:
                if dna_base == 'A':
                    if dna_seq[i] == 'C':
                        if dna_seq[i+1] == 'A':
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            dna_homo_decoded_seq.append(dna_seq[i + 2])
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    elif dna_seq[i] == 'G':
                        if dna_seq[i+1] == 'C':
                            if dna_seq[i+2] == 'A':
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                                dna_homo_decoded_seq.append(dna_seq[i + 2])
                            elif dna_seq[i+2] == 'G':
                                dna_homo_decoded_seq.append(dna_seq[i])
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                            else:
                                raise ValueError('Decoding error, undecodable sequence, please reconfirm.')
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                        i += 2
                elif dna_base == 'C':
                    if dna_seq[i] == 'A':
                        if dna_seq[i + 1] == 'C':
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            dna_homo_decoded_seq.append(dna_seq[i + 2])
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    elif dna_seq[i] == 'T':
                        if dna_seq[i + 1] == 'A':
                            if dna_seq[i + 2] == 'C':
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                                dna_homo_decoded_seq.append(dna_seq[i + 2])
                            elif dna_seq[i + 2] == 'T':
                                dna_homo_decoded_seq.append(dna_seq[i])
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                            else:
                                raise ValueError('Decoding error, undecodable sequence, please reconfirm.')
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                        i += 2
                elif dna_base == 'G':
                    if dna_seq[i] == 'T':
                        if dna_seq[i + 1] == 'G':
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            dna_homo_decoded_seq.append(dna_seq[i + 2])
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    elif dna_seq[i] == 'A':
                        if dna_seq[i + 1] == 'T':
                            if dna_seq[i + 2] == 'G':
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                                dna_homo_decoded_seq.append(dna_seq[i + 2])
                            elif dna_seq[i + 2] == 'A':
                                dna_homo_decoded_seq.append(dna_seq[i])
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                            else:
                                raise ValueError('Decoding error, undecodable sequence, please reconfirm.')
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                        i += 2
                else:
                    if dna_seq[i] == 'G':
                        if dna_seq[i + 1] == 'T':
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            dna_homo_decoded_seq.append(dna_seq[i + 2])
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    elif dna_seq[i] == 'C':
                        if dna_seq[i + 1] == 'G':
                            if dna_seq[i + 2] == 'T':
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                                dna_homo_decoded_seq.append(dna_seq[i + 2])
                            elif dna_seq[i + 2] == 'C':
                                dna_homo_decoded_seq.append(dna_seq[i])
                                dna_homo_decoded_seq.append(dna_seq[i + 1])
                            else:
                                raise ValueError('Decoding error, undecodable sequence, please reconfirm.')
                            i += 3
                        else:
                            dna_homo_decoded_seq.append(dna_seq[i])
                            dna_homo_decoded_seq.append(dna_seq[i + 1])
                            i += 2
                    else:
                        dna_homo_decoded_seq.append(dna_seq[i])
                        dna_homo_decoded_seq.append(dna_seq[i + 1])
                        i += 2

                # upgrade homo_count and dna_base
                if dna_seq[i - 1] != dna_seq[i - 2]:
                    homo_count = 1
                else:
                    homo_count = 2
                dna_base = dna_seq[i - 1]

    return dna_homo_decoded_seq
