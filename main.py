import argparse
from coding_check import transfer_coding_check, homopolymer_coding_check


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
    if opt.codec_type == 'Transfer':
        transfer_coding_check(opt)
    elif opt.codec_type == 'Homopolymer':
        homopolymer_coding_check(opt)
    else:
        raise ValueError("Check the option of codec_type! (Option: Transfer or Homopolymer)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--codec_type", type=str, default='Transfer',
                        help="The type of codec will be checked (Transfer or Homopolymer)")
    parser.add_argument("--dna_l_list", type=int, default=[100],
                        help="The list of DNA sequence lengths")
    parser.add_argument("--dna_num", type=int, default=100000,
                        help="The number of DNA sequences")
    parser.add_argument("--dna_data_file", type=str, default='',
                        help="The file of DNA data")
    parser.add_argument("--generated_type", type=str, default='iid',
                        help="The type of generated DNA data (iid of markov)")
    parser.add_argument("--homo_list", type=list, default=[3],
                        help="The list of original homopolymer constraints")
    parser.add_argument("--homo_t_list", type=list, default=[4],
                        help="The list of target homopolymer constraints")
    parser.add_argument("--random_seed", type=int, default=111,
                        help="The seed of random generator")
    parser.add_argument("--write_encoded_data", type=str2bool, default=False,
                        help="Write encoded data to file")
    parser.add_argument("--write_generated_data", type=str2bool, default=False,
                        help="Write generated data to file")
    parser.add_argument("--cal_gc", type=str2bool, default=False,
                        help="Calculate the ratio of GC content in DNA data")
    parser.add_argument("--encoded_data_filename", type=str, default='data/dna_encoded_data.txt',
                        help="The filename of the encoded data")

    transfer_codec_opt = parser.parse_args()

    main(transfer_codec_opt)
