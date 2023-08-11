import argparse
import os
import random
import datetime
import pandas as pd

from decoder import homo_decoder
from encoder import homopolymer_coding
from utils import (generator_seq, check_ratio, replace_base, make_target_string_list,
                   edit_dist, inverse_replace_base)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(opt):
    target_dict = {'A': ['GCC', 'GCT'], 'C': ['TAA', 'TAG'], 'G': ['ATC', 'ATT'], 'T': ['CGA', 'CGG']}

    homo_list = opt.homo_list
    homo_target_list = opt.homo_t_list
    dna_len_list = opt.dna_l_list
    dna_num = opt.dna_num

    current_date = datetime.datetime.now().strftime('%Y%m%d%H%M')

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
                ori_dna = generator_seq(opt.generated_type, dna_num, dna_len, homo)
                file_name = f"generated/generator_h{homo}2h{homo_target}_l{dna_len}_{dna_num}.txt"
                with open(file_name, "a") as f:
                    for _enc_data_ in ori_dna:
                        f.write(''.join(_enc_data_) + "\n")

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

                    dna_seq_dec = homo_decoder(seqs[i], homo)
                    len_list_dec.append(len(dna_seq_dec))

                    dna_seq_enc = homopolymer_coding(dna_seq_dec, homo_target)
                    len_list_enc.append(len(dna_seq_enc))

                    dna_seq_dec_h_t = homo_decoder(dna_seq_enc, homo_target)
                    dna_seq_enc_h = homopolymer_coding(dna_seq_dec_h_t, homo)

                    dna_seq_reverse = inverse_replace_base(dna_seq_enc_h, homo)
                    dna_seq_ori = ''.join(ori_dna[i])

                    v_ed_dist = edit_dist(dna_seq_ori, ''.join(dna_seq_enc))
                    edit_dist_list.append(v_ed_dist)

                    if dna_seq_ori == dna_seq_reverse:
                        coding_results.append(True)
                    else:
                        coding_results.append(False)

                    homo_num.append(homo)
                    homo_target_num.append(homo_target)
                    len_list_ori.append(dna_len)

                    # pure h coding
                    pure_h_dna_seq_enc = homopolymer_coding(ori_dna[i], homo_target)
                    pure_h_len_list_enc.append(len(pure_h_dna_seq_enc))

                    v_ed_dist_p = edit_dist(dna_seq_ori, ''.join(pure_h_dna_seq_enc))
                    edit_dist_list_pure.append(v_ed_dist_p)

                    pure_h_dna_seq_dec = homo_decoder(pure_h_dna_seq_enc, homo_target)
                    _h3_decoded_ = ''.join(pure_h_dna_seq_dec)
                    if dna_seq_ori == _h3_decoded_:
                        pure_h_results.append(True)
                    else:
                        pure_h_results.append(False)

                os.remove(file_name)
                os.remove(new_file_name)

                dna_data = {'homo': homo_num, 'homo_t': homo_target_num,
                            'length_ori': len_list_ori,
                            'replaced_length': len_list_rep, 'decoding_length': len_list_dec,
                            'encoding_length': len_list_enc, 'coding_result': coding_results,
                            'edit_distance': edit_dist_list,
                            'pure_h3_encoding': pure_h_len_list_enc,
                            'pure_h3_results': pure_h_results, 'edit_distance_pure': edit_dist_list_pure}
                df = pd.DataFrame(data=dna_data)
                check_point = sum(df['coding_result'])
                if check_point == opt.dna_num:
                    print(f'Transfer encoding succeed: Homopolymer {homo} to {homo_target}, Length {dna_len}')

                df.to_csv(f'results/{current_date}_H{homo}toh{homo_target}_dna_length_{dna_len}'
                          f'_num_{dna_num}_{opt.generated_type}.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--homo_list", type=list, default=[3],
                        help="The list of original homopolymer constraints")
    parser.add_argument("--homo_t_list", type=list, default=[4],
                        help="The list of target homopolymer constraints")
    parser.add_argument("--dna_l_list", type=int, default=[100],
                        help="The list of DNA sequence lengths")
    parser.add_argument("--dna_num", type=int, default=10000,
                        help="The number of DNA sequences")
    parser.add_argument("--random_seed", type=int, default=111,
                        help="The seed of random generator")
    parser.add_argument("--generated_type", type=str, default='iid',
                        help="The type of generated DNA data (iid of Markov)")

    dna_2_dna_opt = parser.parse_args()

    main(dna_2_dna_opt)
