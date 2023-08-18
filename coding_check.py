import os
import random
import datetime
import pandas as pd
import time

from decoder import homo_decoder
from encoder import homopolymer_coding
from utils import (generator_seq, check_ratio, replace_base, make_target_string_list,
                   edit_dist, inverse_replace_base, cal_gc_content, cal_homo)


current_date = datetime.datetime.now().strftime('%Y%m%d%H%M')


def transfer_coding_check(opt):
    target_dict = {'A': ['GCC', 'GCT'], 'C': ['TAA', 'TAG'], 'G': ['ATC', 'ATT'], 'T': ['CGA', 'CGG']}

    homo_list = opt.homo_list
    homo_target_list = opt.homo_t_list
    dna_len_list = opt.dna_l_list
    dna_num = opt.dna_num

    for dna_len in dna_len_list:
        for homo in homo_list:
            for homo_target in homo_target_list:
                homo_num = []
                homo_target_num = []
                len_list_dec = []
                len_list_enc = []
                len_list_rep = []
                len_list_ori = []
                coding_results = []
                pure_h_len_list_enc = []
                pure_h_results = []
                edit_dist_list = []
                edit_dist_list_pure = []

                if homo == homo_target:
                    continue

                random.seed(opt.random_seed)
                # Generate the random data
                ori_dna = generator_seq(opt.codec_type, opt.generated_type, dna_num, dna_len, homo)
                file_name = f"generated/generator_h{homo}2h{homo_target}_l{dna_len}_{dna_num}" \
                            f"_{opt.generated_type}.txt"
                with open(file_name, "a") as f:
                    for _enc_data_ in ori_dna:
                        f.write(''.join(_enc_data_) + "\n")

                # start time check point
                start_time = time.time()

                # Replace bases that cannot be decoded
                # Calculate the number of base sequences that cannot be decoded
                target_string = make_target_string_list(homo, target_dict)

                count_ori, _, _ = check_ratio(file_name, target_string)
                print(f'Homopolymer: {homo} to {homo_target}, Length: {dna_len}')
                print(f'\tThe number of lines containing {target_string} is {count_ori}')

                # Calculate the number of base sequences that cannot be decoded after performing the replacing
                new_file_name = replace_base(file_name, homo)
                count_replace, wrong_, _ = check_ratio(new_file_name, target_string)
                print(f'Homopolymer: {homo} to {homo_target}, Length: {dna_len}')
                print(f'\tThe number of lines containing {target_string} is {count_replace}')
                print(f'\tThe wrong sequence is {wrong_}')

                with open(new_file_name, "r") as f:
                    seqs = []
                    for line in f:
                        seq = list(line.rstrip('\n'))
                        seqs.append(seq)

                for i in range(len(seqs)):
                    # get replaced sequence length
                    len_list_rep.append(len(seqs[i]))

                    # get sequence by homopolymer decoder: original h1
                    dna_seq_dec = homo_decoder(seqs[i], homo)
                    len_list_dec.append(len(dna_seq_dec))

                    # get sequence by homopolymer target encoder: target h2
                    dna_seq_enc = homopolymer_coding(dna_seq_dec, homo_target)
                    len_list_enc.append(len(dna_seq_enc))

                    # get sequence by homopolymer target decoder: target h2
                    dna_seq_dec_h_t = homo_decoder(dna_seq_enc, homo_target)

                    # get sequence by homopolymer target encoder: original h1
                    dna_seq_enc_h = homopolymer_coding(dna_seq_dec_h_t, homo)

                    # get sequence by inverse replace function
                    dna_seq_reverse = inverse_replace_base(dna_seq_enc_h, homo)

                    # get original DNA sequence
                    dna_seq_ori = ''.join(ori_dna[i])

                    # get edit distance between original seq and Transfer encoding seq
                    v_ed_dist = edit_dist(dna_seq_ori, ''.join(dna_seq_enc))
                    edit_dist_list.append(v_ed_dist)

                    # Write encoded DNA sequence to file
                    if opt.write_encoded_data:
                        with open(opt.encoded_data_filename, "a") as f:
                            f.write(''.join(dna_seq_enc) + "\n")

                    # check the result of Transfer coding
                    if dna_seq_ori == dna_seq_reverse:
                        coding_results.append(True)
                    else:
                        coding_results.append(False)

                    homo_num.append(homo)
                    homo_target_num.append(homo_target)
                    len_list_ori.append(dna_len)

                    # pure h Homopolymer encoding
                    pure_h_dna_seq_enc = homopolymer_coding(ori_dna[i], homo_target)
                    pure_h_len_list_enc.append(len(pure_h_dna_seq_enc))

                    # get edit distance between original seq and Homopolymer encoding seq
                    v_ed_dist_p = edit_dist(dna_seq_ori, ''.join(pure_h_dna_seq_enc))
                    edit_dist_list_pure.append(v_ed_dist_p)

                    pure_h_dna_seq_dec = homo_decoder(pure_h_dna_seq_enc, homo_target)
                    _h3_decoded_ = ''.join(pure_h_dna_seq_dec)
                    if dna_seq_ori == _h3_decoded_:
                        pure_h_results.append(True)
                    else:
                        pure_h_results.append(False)

                # codec finished time check point
                end_time = time.time()
                print(f'Times: {end_time - start_time}')

                # remove the generated files
                if not opt.write_generated_data:
                    os.remove(file_name)
                os.remove(new_file_name)

                # collecting information of Transfer coding results
                dna_data = {'homo': homo_num, 'homo_t': homo_target_num,
                            'length_ori': len_list_ori,
                            'replaced_length': len_list_rep, 'decoding_length': len_list_dec,
                            'encoding_length': len_list_enc, 'coding_result': coding_results,
                            'edit_distance': edit_dist_list,
                            'pure_h3_encoding': pure_h_len_list_enc,
                            'pure_h3_results': pure_h_results, 'edit_distance_pure': edit_dist_list_pure}
                df = pd.DataFrame(data=dna_data)

                # check the coding results
                check_point = sum(df['coding_result'])
                if check_point == opt.dna_num:
                    print(f'Transfer encoding succeed: Homopolymer {homo} to {homo_target},'
                          f' Length {dna_len}, DNA num {dna_num}')
                else:
                    print(f'Transfer encoding failed: Homopolymer {homo} to {homo_target},'
                          f' Length {dna_len}, DNA num {dna_num}')

                # recoding the results
                df.to_csv(f'results/{current_date}_H{homo}toh{homo_target}_dna_length_{dna_len}'
                          f'_num_{dna_num}_{opt.generated_type}.csv')


def homopolymer_coding_check(opt):
    dna_homo_list = opt.homo_list
    dna_len_list = opt.dna_l_list
    dna_num = opt.dna_num
    for dna_homo in dna_homo_list:
        for dna_len in dna_len_list:
            if opt.dna_data_file:
                print(f'The experiment of homo: {dna_homo} data: {dna_len}, {dna_num}, '
                      f'{opt.dna_data_file}')
            else:
                print(f'The experiment of homo: {dna_homo} data: {dna_len}, {dna_num}, '
                      f'{opt.generated_type}')
            encoding_change = []
            enc_dec_success = []
            enc_length = []
            enc_homo_max = []
            ori_gc_content_percent = []
            enc_gc_content_percent = []
            edit_dist_list = []

            # read DNA data if it already exists
            if opt.dna_data_file:
                dna_ori = []
                with open(opt.dna_data_file, "r") as f:
                    for line in f:
                        seq = list(line.rstrip('\n'))
                        dna_ori.append(seq)

            # else generate the DNA data
            else:
                dna_ori = generator_seq(opt.codec_type, opt.generated_type, dna_num, dna_len, dna_homo)
                # write data to txt file
                if opt.write_generated_data:
                    file_name = f"generated/generator_h{dna_homo}_l{dna_len}_{dna_num}" \
                                f"_{opt.generated_type}.txt"
                    with open(file_name, "a") as f:
                        for _data_ in dna_ori:
                            f.write(''.join(_data_) + "\n")

            # start time check point
            start_time = time.time()

            for j in range(len(dna_ori)):
                # Encoding
                dna_seq_enc = homopolymer_coding(dna_ori[j], dna_homo)
                if dna_ori[j] == dna_seq_enc:
                    # print(f'\tIs the sequence the same before and after encoding? Yes')
                    encoding_change.append('No')
                else:
                    encoding_change.append('Yes')
                    # print(f'\tIs the sequence the same before and after encoding? No')

                # Calculate edit distance
                # get original DNA sequence
                dna_seq_ori = ''.join(dna_ori[j])
                v_ed_dist = edit_dist(dna_seq_ori, ''.join(dna_seq_enc))
                edit_dist_list.append(v_ed_dist)

                # Write encoded DNA sequence to file
                if opt.write_encoded_data:
                    with open(opt.encoded_data_filename, "a") as f:
                        f.write(''.join(dna_seq_enc) + "\n")

                # Decoding
                dna_seq_dec = homo_decoder(dna_seq_enc, dna_homo)
                if dna_ori[j] == dna_seq_dec:
                    enc_dec_success.append(True)
                    # print(f'\tWas decoding successful? Yes')
                else:
                    enc_dec_success.append(False)
                    # print(f'\tWas decoding successful? No')

                current_dna_len = len(dna_seq_enc)
                homo_max = cal_homo(dna_seq_enc)
                enc_length.append(current_dna_len)
                enc_homo_max.append(homo_max)

                if opt.cal_gc:
                    enc_gc_content = cal_gc_content(dna_seq_enc)
                    enc_gc_content_percent.append(enc_gc_content[0]*100)
                    ori_gc_content = cal_gc_content(dna_ori[j])
                    ori_gc_content_percent.append(ori_gc_content[0]*100)

            # start time check point
            end_time = time.time()
            times = end_time - start_time
            print(f'Time: {times}s')

            dna_data = {'changed': encoding_change, 'decoding_success': enc_dec_success,
                        'encoded_DNA_length': enc_length, 'Homopolymer': enc_homo_max,
                        'edit_distance': edit_dist_list}
            if opt.cal_gc:
                dna_data.setdefault('original_GC_content', ori_gc_content_percent)
                dna_data.setdefault('encoded_GC_content', enc_gc_content_percent)
            df = pd.DataFrame(data=dna_data)

            _mean_ = df['encoded_DNA_length'].mean()
            print(f'The mean of encoded DNA length is {_mean_}')
            df.to_csv(f'results/{current_date}_H_{dna_homo}_No_{dna_num}_dna_length_{dna_len}'
                      f'_{opt.generated_type}.csv')
