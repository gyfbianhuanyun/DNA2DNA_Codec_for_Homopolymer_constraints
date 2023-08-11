import random
import editdistance as ed
import numpy as np


def edit_dist(s1, s2):
    # Get edit distance between s1 and s2
    # Import edit distance function
    v_dist = np.vectorize(ed.eval)
    # Calculate the distance
    dist = v_dist(s1, s2)
    # Get mean of distance
    mean_dist = np.mean(dist)
    return mean_dist


# check homopolymer
def cal_homo(dna_seq):
    dna_base = ''
    homo_count = 1
    homo_max = 1
    for i in range(len(dna_seq)):
        dna_1_seq = dna_seq[i]
        for j in range(len(dna_1_seq)):
            if dna_1_seq[j] == dna_base:
                homo_count += 1
            else:
                dna_base = dna_1_seq[j]
                homo_count = 1
            if homo_count > homo_max:
                homo_max = homo_count
    return homo_max


# check GC content
def cal_gc_content(dna_data):
    dna_data_array = np.array(dna_data)
    gc_list = []
    for i in range(len(dna_data)):
        # Count the number of bases 'C' and 'G'
        bases_, count_ = np.unique(dna_data_array[i], return_counts=True)

        if count_[np.where(bases_ == 'C')]:
            count_c = count_[np.where(bases_ == 'C')]
        else:
            count_c = 0
        if count_[np.where(bases_ == 'G')]:
            count_g = count_[np.where(bases_ == 'G')]
        else:
            count_g = 0

        _gc_count = count_c + count_g
        _gc_percent = _gc_count / sum(count_)
        gc_list.append(_gc_percent)
    gc_percent = sum(gc_list) / len(gc_list)
    return gc_percent


# replace_function
def replace_func(dna_seq, homopolymer, target_seqs, replaced_seqs):
    # dna_seq: list, homopolymer: int
    homo_count = 1
    dna_base = ''
    dna_check = ''
    dna_replaced_seq = []
    i = 0
    while i < len(dna_seq):
        # check homopolymer value
        if homo_count != homopolymer:
            if dna_seq[i] == dna_base:
                homo_count += 1
                dna_check += dna_base
            else:
                dna_base = dna_seq[i]
                homo_count = 1
                dna_check = dna_base
            dna_replaced_seq.append(dna_seq[i])
            # upgrade i
            i += 1
        else:
            if i + 1 == len(dna_seq):
                dna_replaced_seq.append(dna_seq[i])
                i += 1
            elif i + 2 == len(dna_seq):
                if dna_base == 'A':
                    if dna_seq[i] == 'C' and dna_seq[i + 1] == 'A':
                        last_bases = 'G' + dna_seq[i] + dna_seq[i + 1]
                    elif dna_seq[i] == 'G' and dna_seq[i + 1] == 'C':
                        last_bases = dna_seq[i] + dna_seq[i + 1] + 'G'
                    else:
                        last_bases = dna_seq[i] + dna_seq[i + 1]
                elif dna_base == 'C':
                    if dna_seq[i] == 'A' and dna_seq[i + 1] == 'C':
                        last_bases = 'T' + dna_seq[i] + dna_seq[i + 1]
                    elif dna_seq[i] == 'T' and dna_seq[i + 1] == 'A':
                        last_bases = dna_seq[i] + dna_seq[i + 1] + 'T'
                    else:
                        last_bases = dna_seq[i] + dna_seq[i + 1]
                elif dna_base == 'G':
                    if dna_seq[i] == 'T' and dna_seq[i + 1] == 'G':
                        last_bases = 'A' + dna_seq[i] + dna_seq[i + 1]
                    elif dna_seq[i] == 'A' and dna_seq[i + 1] == 'T':
                        last_bases = dna_seq[i] + dna_seq[i + 1] + 'A'
                    else:
                        last_bases = dna_seq[i] + dna_seq[i + 1]
                else:
                    if dna_seq[i] == 'G' and dna_seq[i + 1] == 'T':
                        last_bases = 'C' + dna_seq[i] + dna_seq[i + 1]
                    elif dna_seq[i] == 'C' and dna_seq[i + 1] == 'G':
                        last_bases = dna_seq[i] + dna_seq[i + 1] + 'C'
                    else:
                        last_bases = dna_seq[i] + dna_seq[i + 1]
                dna_replaced_seq.append(last_bases)
                i += 2
            else:
                dna_check += dna_seq[i] + dna_seq[i + 1] + dna_seq[i + 2]
                # print(dna_check)
                if dna_check in target_seqs:
                    index = target_seqs.index(dna_check)
                    replaced_seq = replaced_seqs[index]
                    replaced_seq = replaced_seq[homopolymer:]
                    dna_replaced_seq.append(replaced_seq)
                    # back to initial state
                    homo_count = 1
                    dna_base = dna_seq[i + 2]
                    dna_check = dna_base
                    # upgrade i
                    i += 3
                else:
                    dna_replaced_seq.append(dna_seq[i])
                    # upgrade state
                    homo_count = 1
                    dna_base = dna_seq[i]
                    dna_check = dna_base
                    # upgrade i
                    i += 1
        # print(i, '/', len(dna_seq), dna_replaced_seq)
        # print('\t', homo_count, dna_base, dna_check)

    return dna_replaced_seq


# open file and change the undecodable string
def replace_base(filename, homo):
    # Get replace dna sequence 1
    replace_dna_seq_1 = []
    target_dna_dict_1 = {'A': ['GCG'], 'C': ['TAT'], 'G': ['ATA'], 'T': ['CGC']}
    target_dna_seq_1 = make_target_string_list(homo, target_dna_dict_1)
    for target_seq in target_dna_seq_1:
        if target_seq[0] == 'A':
            add_base = 'A'
        elif target_seq[0] == 'C':
            add_base = 'C'
        elif target_seq[0] == 'G':
            add_base = 'G'
        else:
            add_base = 'T'

        re_base = target_seq[:homo + 2] + add_base + target_seq[homo + 2:]
        replace_dna_seq_1.append(re_base)

    # Get replace dna sequence 2
    replace_dna_seq_2 = []
    target_dna_dict_2 = {'A': ['GCC', 'GCT', 'GCA'], 'C': ['TAA', 'TAG', 'TAC'],
                         'G': ['ATC', 'ATT', 'ATG'], 'T': ['CGA', 'CGG', 'CGT']}
    target_dna_seq_2 = make_target_string_list(homo, target_dna_dict_2)
    for target_seq in target_dna_seq_2:
        if target_seq[0] == 'A':
            add_base = 'G'
        elif target_seq[0] == 'C':
            add_base = 'T'
        elif target_seq[0] == 'G':
            add_base = 'A'
        else:
            add_base = 'C'
        re_base = target_seq[:homo+2] + add_base + target_seq[homo+2:]
        replace_dna_seq_2.append(re_base)

    target_dna_seq = target_dna_seq_2 + target_dna_seq_1
    replace_dna_seq = replace_dna_seq_2 + replace_dna_seq_1

    # 打开文件并读取所有行
    with open(filename, "r") as f:
        seqs = []
        for line in f:
            seq = list(line.rstrip('\n'))
            seqs.append(seq)

    replaced_seqs = []
    for i in range(len(seqs)):
        seq = seqs[i]
        new_seq = replace_func(seq, homo, target_dna_seq, replace_dna_seq)
        replaced_seqs.append(new_seq)

    # 将新的字符串写入到新文件中
    new_filename = filename.split('.')[0] + '_new.txt'
    with open(new_filename, "w") as new_file:
        for item in replaced_seqs:
            new_file.write(''.join(item) + '\n')

    return new_filename


# open file, read and count the number of lines by line
def check_ratio(file, string_to_search_list):
    right_list = []
    wrong_list = []
    with open(file, 'r') as file:
        _count = 0
        for line in file:
            flag = _count
            for string_to_search in string_to_search_list:
                if string_to_search in line:
                    _count += 1
                    right_list.append(line.strip('\n'))
                    break
            if flag == _count:
                wrong_list.append(line.strip('\n'))

    return _count, right_list, wrong_list


def make_target_string_list(homo, string_map):
    target_string_list = []
    dna_base = ['A', 'C', 'G', 'T']
    for base in dna_base:
        if base in string_map:
            _val_ = string_map[base]
            for _seq_ in _val_:
                target_seq = base * homo
                target_seq += _seq_
                target_string_list.append(target_seq)

    return target_string_list


# inverse replace
def inverse_replace_base(dna_seq, homo):
    seq_str = ''.join(dna_seq)

    # Get replace dna sequence
    replace_dna_seq_1 = []
    target_dna_dict_1 = {'A': ['GCGC', 'GCGT', 'GCGA'], 'C': ['TATA', 'TATG', 'TATC'],
                         'G': ['ATAC', 'ATAT', 'ATAG'], 'T': ['CGCA', 'CGCG', 'CGCT']}
    target_dna_seq_1 = make_target_string_list(homo, target_dna_dict_1)

    for target_seq in target_dna_seq_1:
        re_base = target_seq[:homo+2] + target_seq[homo+3:]
        replace_dna_seq_1.append(re_base)

    replace_dna_seq_2 = []
    target_dna_dict_2 = {'A': ['GCAG'], 'C': ['TACT'],
                         'G': ['ATGA'], 'T': ['CGTC']}
    target_dna_seq_2 = make_target_string_list(homo, target_dna_dict_2)

    for target_seq in target_dna_seq_2:
        re_base = target_seq[:homo + 2] + target_seq[homo + 3:]
        replace_dna_seq_2.append(re_base)

    replace_dna_seq = replace_dna_seq_1 + replace_dna_seq_2
    target_dna_seq = target_dna_seq_1 + target_dna_seq_2

    seq_inverse = inverse_func(list(seq_str), homo, target_dna_seq, replace_dna_seq)
    seq_str = ''.join(seq_inverse)

    return seq_str


# inverse function
def inverse_func(dna_seq, homopolymer, target_seqs, replaced_seqs):
    # dna_seq: list, homopolymer: int
    homo_count = 1
    dna_base = ''
    dna_check = ''
    dna_replaced_seq = []
    i = 0
    while i < len(dna_seq):
        # check homopolymer value
        if homo_count != homopolymer:
            if dna_seq[i] == dna_base:
                homo_count += 1
                dna_check += dna_base
            else:
                dna_base = dna_seq[i]
                homo_count = 1
                dna_check = dna_base
            dna_replaced_seq.append(dna_seq[i])
            # upgrade i
            i += 1
        else:
            if i + 1 == len(dna_seq):
                dna_replaced_seq.append(dna_seq[i])
                i += 1
            elif i + 2 == len(dna_seq):
                last_bases = dna_seq[i] + dna_seq[i + 1]
                dna_replaced_seq.append(last_bases)
                i += 2
            elif i + 3 == len(dna_seq):
                check_seq = dna_seq[i] + dna_seq[i + 1] + dna_seq[i + 2]
                if dna_base == 'A':
                    if check_seq == 'GCA':
                        last_bases = 'CA'
                    elif check_seq == 'GCG':
                        last_bases = 'GC'
                    else:
                        last_bases = check_seq
                elif dna_base == 'C':
                    if check_seq == 'TAC':
                        last_bases = 'AC'
                    elif check_seq == 'TAT':
                        last_bases = 'TA'
                    else:
                        last_bases = check_seq
                elif dna_base == 'G':
                    if check_seq == 'ATA':
                        last_bases = 'AT'
                    elif check_seq == 'ATG':
                        last_bases = 'TG'
                    else:
                        last_bases = check_seq
                else:
                    if check_seq == 'CGC':
                        last_bases = 'CG'
                    elif check_seq == 'CGT':
                        last_bases = 'GT'
                    else:
                        last_bases = check_seq
                i += 3
                dna_replaced_seq.append(last_bases)

            else:
                dna_check += dna_seq[i] + dna_seq[i + 1] + dna_seq[i + 2] + dna_seq[i + 3]
                # print(dna_check)
                if dna_check in target_seqs:
                    index = target_seqs.index(dna_check)
                    replaced_seq = replaced_seqs[index]
                    replaced_seq = replaced_seq[homopolymer:]
                    dna_replaced_seq.append(replaced_seq)

                    # back to initial state
                    homo_count = 1
                    dna_base = dna_seq[i + 3]
                    dna_check = dna_base

                    # upgrade i
                    i += 4

                else:
                    dna_replaced_seq.append(dna_seq[i])
                    # upgrade state
                    homo_count = 1
                    dna_base = dna_seq[i]
                    dna_check = dna_base
                    # upgrade i
                    i += 1

    return dna_replaced_seq


def gen_homo_seq(option, dna_length, homo_cons):
    seq = []
    base_dict = ['A', 'C', 'G', 'T']
    base_count = 1
    last_base = ''
    for i in range(dna_length):
        if base_count == homo_cons:
            base_dict.remove(seq[-1])
            base = random.choice(base_dict)
            base_dict.append(seq[-1])
        else:
            if option == 'iid':
                base = random.choice(base_dict)
            elif option == 'markov':
                if i == 0:
                    base = random.choice(base_dict)
                else:
                    weights = [0.4 if char == last_base else 0.2 for char in base_dict]
                    base = random.choices(base_dict, weights=weights)[0]
            else:
                raise ValueError("Check the option of Generating DNA data! (Option: iid or markov)")
        seq.append(base)
        if base == last_base:
            base_count += 1
        else:
            last_base = base
            base_count = 1

    return seq


def generator_seq(option, num_seq, seq_length, homo_cons):
    dna_seqs = []
    for i in range(num_seq):
        dna_seqs.append(gen_homo_seq(option, seq_length, homo_cons))
    return dna_seqs