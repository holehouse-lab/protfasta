#!/usr/bin/env python

###############################################

import random
import sys
from os import path
import protfasta
import argparse 
from argparse import RawTextHelpFormatter
from protfasta import __version__ as VERSION_MAJ


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================


def stdout(msg, silent):
    if not silent:
        print(msg)
        

####################################################################################################
#
#
def exit_error(msg):
    print('[FATAL ERROR]: %s'%(msg))
    exit(1)

####################################################################################################
#
#
def validate(instring, option):
    if instring.lower() not in option:
        exit_error('Could not find [%s] in list of valid options [%s]'%(instring, str(option)))

    return instring.lower()

####################################################################################################
#
#
def validate_int(val, min_val, param_name):
    try:
        val_i = int(val)

        if val_i < min_val:
            raise Exception
    except Exception:
        exit_error('%s must be a numerical value > %i'%(param_name, min_val))

    return val_i

def print_statistical_summary(data):
    
    length_list=[]
    for d in data:
        length_list.append(len(d[1]))

    length_list.sort()
    q25 = length_list[int(len(length_list)*0.25)]
    q50 = length_list[int(len(length_list)*0.50)]
    q75 = length_list[int(len(length_list)*0.75)]
    print('[STATS]: Total number of sequences : %i ' % (len(length_list)))
    print('[STATS]: 25th percentile lengt     : %i ' % (q25))
    print('[STATS]: Median length             : %i ' % (q50))
    print('[STATS]: 75th percentile lengt     : %i ' % (q75))
    print('[STATS]: Longest sequence          : %i ' % (length_list[-1]))
    print('[STATS]: Shortes sequence          : %i ' % (length_list[0]))

    

        
####################################################################################################
#
#
def main():

    dsc='pfasta is a simple command-line tool for parsing, sanitizing, and manipulating\nprotein-based FASTA files. It is the command-line utility from the package protfasta'

    parser = argparse.ArgumentParser(description=dsc,formatter_class=RawTextHelpFormatter)

    # note nargs means EITHER 0 or 1 arguments are accepted
    parser.add_argument("filename", nargs='?', help='Input FASTA file')

    parser.add_argument("-o", help="Output fasta file (is created)") 
    parser.add_argument("--non-unique-header", help="", action='store_true') 
    parser.add_argument("--duplicate-record", help="How to deal with duplicate records in the file.\nOptions are ['ignore', 'fail', 'remove'] (default = fail)") 
    parser.add_argument("--duplicate-sequence", help="How to deal with duplicate sequences in the file.\nOptions are ['ignore', 'fail', 'remove'] (default = ignore)") 
    parser.add_argument("--invalid-sequence", help="How to deal with duplicate sequences in the file. Available options are shown\nbelow and described (default = fail)\n\nignore : skip invalid residues \nfail   : throw exception on invalid sequences \nremove : remove sequences with invalid characters \nconvert-all : Convert B->N, U->C, X->G, Z->Q, '*'->'',\n             '-'->'' (and throw exception if remaining invalid characters exist)\nconvert-res : same as convert-all except ignore alignment character '-'\nconvert-all-ignore - same as convert-all except invalid characters left over are ignored.\nconvert-res-ignore - same as convert-res except invalid characters left over are ignored.   \nconvert-all-remove - same as convert-all except sequences with invalid characters are removed.\nconvert-res-remove - same as convert-res except sequences with invalid characters are removed.   ")  
    parser.add_argument("--number-lines", help="Number of lines for FASTA file") 
    parser.add_argument("--shortest-seq", help="Shortest sequence included ") 
    parser.add_argument("--longest-seq", help="Longest sequence included") 
    parser.add_argument("--random-subsample", help="Randomly sub-sample from")
    parser.add_argument("--print-statistics", help="Print information on the sequences",action='store_true') 
    parser.add_argument("--no-outputfile", help="Prevents pfasta from writing an output file ",action='store_true') 
    parser.add_argument("--silent", help="Generate no output at all to STDOUT", action='store_true') 
    parser.add_argument("--remove-comma-from-header", help="Flag that, if set, commas in FASTA headers will be replaced with ';' characters. Useful if you have code that parses FASTA headers as part of a CSV file", action='store_true') 
    parser.add_argument("--version", help="Flag that, if set, means we just print the version and exit.", action='store_true') 

    args = parser.parse_args()
    silent = args.silent

    if args.version:
        print(VERSION_MAJ)
        sys.exit(0)

    # this behavior phenocopies the default behavior if we did not allow --version to overide
    if not args.filename:
        parser.error("pfasta: error: the following arguments are required: filename")

    
    if not silent:

        print('')
        print('........................')
        print(f'pfasta version {VERSION_MAJ}')  
        print('Please report bugs to:\nhttps://github.com/holehouse-lab/protfasta')
        print('........................')


    if not path.exists(args.filename):
        exit_error('File %s does not exist'%(args.filename))
        

    # sanitize and set  outfile
    if args.o:
        outfile = args.o
    else:
        outfile = 'output.fasta'


    # sanitize and set uniuqe header
    if args.non_unique_header:
        expect_unique_header = False
    else:
        expect_unique_header = True

    # sanitize duplicate re
    if args.duplicate_record:
        duplicate_record = validate(args.duplicate_record, ['ignore','fail','remove'])
    else:
        duplicate_record = 'fail'

    if args.duplicate_sequence:
        duplicate_sequence = validate(args.duplicate_sequence, ['ignore','fail','remove'])
    else:
        duplicate_sequence = 'ignore'

    if args.invalid_sequence:
        invalid_sequence = validate(args.invalid_sequence, ['ignore',
                                                            'fail',
                                                            'remove',
                                                            'convert-all',
                                                            'convert-res',
                                                            'convert-all-ignore',
                                                            'convert-res-ignore',
                                                            'convert-all-remove',
                                                            'convert-res-remove'])


        if invalid_sequence == 'convert-all':
            invalid_sequence = 'convert'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':'',
                               '-':''}

        elif invalid_sequence == 'convert-res':
            invalid_sequence = 'convert'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':''}

        elif invalid_sequence == 'convert-all-ignore':
            invalid_sequence = 'convert-ignore'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':'',
                               '-':''}

        elif invalid_sequence == 'convert-res-ignore':
            invalid_sequence = 'convert-ignore'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':''}

        elif invalid_sequence == 'convert-all-remove':
            invalid_sequence = 'convert-remove'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':'',
                               '-':''}

        elif invalid_sequence == 'convert-res-remove':
            invalid_sequence = 'convert-remove'
            correction_dict = {'B':'N',
                               'U':'C',
                               'X':'G',
                               'Z':'Q',
                               '*':''}

        else:
            correction_dict = None
            

            
    else:
        invalid_sequence = 'fail'
        correction_dict = None

    if args.number_lines:
        number_of_lines = validate_int(args.number_lines, 5, '--number-lines')
    else:
        number_of_lines = 60

    if args.shortest_seq:
        shortest = validate_int(args.shortest_seq, 0, '--shortest-seq')
    else:
        shortest = False

    if args.longest_seq:
        longest = validate_int(args.longest_seq, 0, '--longest-seq')
    else:
        longest = False

    # sanitize and set uniuqe header
    if args.print_statistics:
        print_stats = True
    else:
        print_stats = False

    if args.random_subsample:
        random_subsample = validate_int(args.random_subsample, 0, '--random-subsample')
    else:
        random_subsample = False

    if longest and shortest:
        if longest < shortest:
            exit_error('--longest-seq must be longer than --shortest-seq')

    if args.remove_comma_from_header:
        def hp(s):
            return s.replace(',',';')
    else:
        hp = None
        


    # read in
    stdout('', silent)
    stdout('[INFO]: Reading in the file %s'%(args.filename), silent)
    if silent:
        verb=False
    else:
        verb=True
    data = protfasta.read_fasta(args.filename, 
                                expect_unique_header=expect_unique_header,
                                header_parser=hp,
                                duplicate_sequence_action=duplicate_sequence,
                                duplicate_record_action=duplicate_record,
                                invalid_sequence_action=invalid_sequence,
                                correction_dictionary = correction_dict,
                                return_list=True,
                                verbose=verb)
       
    
    # if length filters are 
    if longest:
        stdout('[INFO]: Filtering out sequences longer than %s'%(longest),silent)
        tmp = []
        for i in data:
            if len(i[1]) < longest:
                tmp.append(i)
        data = tmp
    
    if shortest:
        stdout('[INFO]: Filtering out sequences shorter than %s'%(shortest), silent)
        tmp = []
        for i in data:
            if len(i[1]) > shortest:
                tmp.append(i)
        data = tmp

    if len(data) < 1:
        if longest or shortest:
            stdout('[INFO]: 0 sequences found that match the longes/shortest filtering criterion',silent)
            sys.exit(0)
        else:
            stdout('[INFO]: 0 sequences found that match the longes/shortest filtering criterion',silent)
            sys.exit(0)


    if random_subsample:
        if len(data) < random_subsample:
            stdout('[INFO]: Cannot subsample as the requested number to subsample (%i) is more than\n        the total number of sequences (%i). Using all sequences'%(random_subsample, len(data)), silent)
        else:
            stdout('[INFO]: Subsampling %i sequences from the complete dataset (%i)'%(random_subsample, len(data)), silent)
        tmp = []
        x = list(range(0,len(data)))
        random.shuffle(x)
        idx = x[0:random_subsample]
        for i in idx:
            tmp.append(data[i])
        data = tmp



    if print_stats and not silent:
        print_statistical_summary(data)

    if args.no_outputfile:
        stdout('[INFO]: No outputfile requested ',silent)
    else:
        stdout('[INFO]: Writing new sequence file [%s]...'%(outfile),silent)
        protfasta.write_fasta(data, outfile, linelength=number_of_lines)

        
                   
            

