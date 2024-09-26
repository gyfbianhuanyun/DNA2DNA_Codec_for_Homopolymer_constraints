# DNA2DNA_Codec_for_Homopolymer_constraints
Adaptable DNA Storage Coding: An Efficient Framework for Homopolymer Constraint Transitions

## Download
Clone code from GitHub through ssh
```
git@github.com:gyfbianhuanyun/DNA2DNA_Codec_for_Homopolymer_constraints.git
```

## Encoding Processing
DNA to DNA Transfer encoding

1. Random generate DNA data using i.i.d or Markov Chain
2. Encoding

    A. Homopolymer constraint encoding

    B. Transfer encoding

    C. Check results

## Code
### Structure
DNA to DNA Coding
```
Main python file
1. main.py
    Codec runs Python file containing arguments 
2. coding_check.py
    The encoding processing with options (Transfer of Homopolymer)
3. encoder.py
    Encoder Function
4. decoder.py
    Decoder Function

```

Run codec
```
Python main.py --options information
```

### Options settings
```
# Coding method
--codec_type: The type of codec will be used (Transfer or Homopolymer)

# Data information
--dna_l_list: The list of DNA sequence lengths
--dna_num: The number of DNA sequences
--dna_data_file: The filename of DNA data
--generated_type: The generated DNA data (iid of Markov)

# DNA storage channel constraints
--homo_list: The list of original homopolymer constraints
--homo_t_list: The list of target homopolymer constraints

# Generated file in codec processing
--write_encoded_data: Write encoded data to file or not
--write_generated_data: Write generated data to file or not
--encoded_data_filename: The filename of the encoded data

# Other codec processing options
--random_seed: The seed of random generator
--cal_gc: Calculate the ratio of GC content in DNA data
```

### Example
#### 1. Using Homopolymer encoding to encode DNA data

DNA sequence length is 100, DNA data filename is 'DNA_data.txt', target homopolymer constraint is 3, and Encoded DNA data is 'Enc_DNA.txt'.
```
python main.py --codec_type Homopolymer --dna_l_list [100] --dna_data_file 'DNA_data.txt' --write_encoded_data True --encoded_data_filename 'Enc_DNA.txt' --homo_list [3]
```

#### 2. Using Transfer encoding to encode DNA data

DNA sequence length is 100, DNA data filename is 'DNA_data.txt', original homopolymer constraint is 4, target homopolymer constraint is 3, and Encoded DNA data is 'Trans_DNA.txt'.
```
python main.py --codec_type Transfer --dna_l_list [100] --dna_data_file 'DNA_data.txt' --write_encoded_data True --encoded_data_filename 'Trans_DNA.txt' --homo_list [4] --homo_t_list [3]
```

#### 3. If you need to generate random DNA sequences

The number of DNA bases is 10,000, and the random generation seed is 111.

##### A. Independent and identically distributed sequence generation (i.i.d)

Using Homopolymer encoding to encode DNA data, DNA sequence length is 100, target homopolymer constraint is 3, and Encoded DNA data is 'Enc_DNA.txt'.
```
python main.py --codec_type Homopolymer --dna_l_list [100] --dna_num 10000 --generated_type iid --homo_list [3] --write_encoded_data True --write_generated_data True --encoded_data_filename ''Enc_DNA.txt'' --random_seed 111
```

##### B. Markov sequence generation (Markov)

Using Transfer encoding to encode DNA data, DNA sequence length is 100, original homopolymer constraint is 4, target homopolymer constraint is 3, and Encoded DNA data is 'Trans_DNA.txt'.
```
python main.py --codec_type Transfer --dna_l_list [100] --dna_num 10000 --generated_type markov --homo_list [4] --homo_t_list [3] --write_encoded_data True --write_generated_data True --encoded_data_filename ''Trans_DNA.txt'' --random_seed 111
```
