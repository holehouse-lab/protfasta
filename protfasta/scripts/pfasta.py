#!/usr/bin/env python
"""Command-line interface for protfasta.

``pfasta`` is a small CLI tool for parsing, sanitizing, and
manipulating protein FASTA files.  It is installed as an entry-point
when the *protfasta* package is installed.
"""

from __future__ import annotations

import random
import sys
from os import path
from typing import Callable, Optional, cast

import protfasta
import argparse
from argparse import RawTextHelpFormatter
from protfasta import __version__ as VERSION_MAJ


## ===================================================================================================
##                              Main Script - hold onto your hat!
## ===================================================================================================


def stdout(msg: str, silent: bool) -> None:
    """Print *msg* to stdout unless *silent* is ``True``.

    Parameters
    ----------
    msg : str
        Message to print.
    silent : bool
        If ``True``, suppress output.
    """
    if not silent:
        print(msg)
        

####################################################################################################
#
#
def exit_error(msg: str) -> None:
    """Print a fatal-error message and terminate the process.

    Parameters
    ----------
    msg : str
        Error description shown to the user.
    """
    print('[FATAL ERROR]: %s' % (msg))
    exit(1)

####################################################################################################
#
#
def validate(instring: str, option: list[str]) -> str:
    """Validate and normalise a CLI option against allowed values.

    Parameters
    ----------
    instring : str
        The value supplied by the user.
    option : list[str]
        Acceptable lower-case values.

    Returns
    -------
    str
        The lower-cased input if it is valid; otherwise the process
        exits via :func:`exit_error`.
    """
    if instring.lower() not in option:
        exit_error('Could not find [%s] in list of valid options [%s]' % (instring, str(option)))

    return instring.lower()

####################################################################################################
#
#
def validate_int(val: str, min_val: int, param_name: str) -> int:
    """Parse a string as an integer and enforce a minimum value.

    Parameters
    ----------
    val : str
        The raw string from the command line.
    min_val : int
        The minimum acceptable integer value (inclusive).
    param_name : str
        Name of the parameter (used in error messages).

    Returns
    -------
    int
        The validated integer.
    """
    try:
        val_i = int(val)

        if val_i < min_val:
            raise Exception
    except Exception:
        exit_error('%s must be a numerical value > %i'%(param_name, min_val))

    return val_i

def print_statistical_summary(data: list[list[str]]) -> None:
    """Print basic length statistics for a set of sequences.

    Parameters
    ----------
    data : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.
    """
    length_list: list[int] = []
    for d in data:
        length_list.append(len(d[1]))

    if not length_list:
        print('[STATS]: Total number of sequences : 0 ')
        return

    length_list.sort()
    q25 = length_list[int(len(length_list)*0.25)]
    q50 = length_list[int(len(length_list)*0.50)]
    q75 = length_list[int(len(length_list)*0.75)]
    print('[STATS]: Total number of sequences : %i ' % (len(length_list)))
    print('[STATS]: 25th percentile length    : %i ' % (q25))
    print('[STATS]: Median length             : %i ' % (q50))
    print('[STATS]: 75th percentile length    : %i ' % (q75))
    print('[STATS]: Longest sequence          : %i ' % (length_list[-1]))
    print('[STATS]: Shortest sequence         : %i ' % (length_list[0]))

    

        
####################################################################################################
#
#
def main() -> None:
    """Entry-point for the ``pfasta`` command-line tool.

    Parses arguments, reads the input FASTA file via
    :func:`protfasta.read_fasta`, applies any requested filters
    (length, sub-sampling, duplicate/invalid handling), and writes
    the result with :func:`protfasta.write_fasta`.
    """

    dsc = ('pfasta is a simple command-line tool for parsing, sanitizing, and manipulating\n'
           'protein-based FASTA files. It is the command-line utility from the package protfasta')

    parser = argparse.ArgumentParser(description=dsc, formatter_class=RawTextHelpFormatter)

    # note nargs means EITHER 0 or 1 arguments are accepted
    parser.add_argument("filename", nargs='?', help='Input FASTA file')

    parser.add_argument("-o", help="Output fasta file (is created)") 
    parser.add_argument("--non-unique-header", help="", action='store_true') 
    parser.add_argument("--duplicate-record", help="How to deal with duplicate records in the file.\nOptions are ['ignore', 'fail', 'remove'] (default = fail)") 
    parser.add_argument("--duplicate-sequence", help="How to deal with duplicate sequences in the file.\nOptions are ['ignore', 'fail', 'remove'] (default = ignore)") 
    parser.add_argument("--invalid-sequence", help="How to deal with invalid (non-standard) residues in the file. Available options\nare shown below and described (default = fail)\n\nignore : skip invalid residues \nfail   : throw exception on invalid sequences \nremove : remove sequences with invalid characters \nconvert-all : Convert B->N, U->C, X->G, Z->Q, '*'->'',\n             '-'->'' (and throw exception if remaining invalid characters exist)\nconvert-res : same as convert-all except ignore alignment character '-'\nconvert-all-ignore - same as convert-all except invalid characters left over are ignored.\nconvert-res-ignore - same as convert-res except invalid characters left over are ignored.   \nconvert-all-remove - same as convert-all except sequences with invalid characters are removed.\nconvert-res-remove - same as convert-res except sequences with invalid characters are removed.   ")  
    parser.add_argument("--number-lines", help="Number of residues per line in the output FASTA file (default = 60)")
    parser.add_argument("--shortest-seq", help="Minimum length filter; sequences shorter than or equal to this length are discarded")
    parser.add_argument("--longest-seq", help="Maximum length filter; sequences longer than or equal to this length are discarded")
    parser.add_argument("--random-subsample", help="Randomly sub-sample this many sequences from the final set")
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

    # note we compare against None (rather than relying on truthiness) so that
    # an explicitly-passed value of 0 is honoured rather than silently ignored
    if args.shortest_seq is not None:
        shortest = validate_int(args.shortest_seq, 0, '--shortest-seq')
    else:
        shortest = None

    if args.longest_seq is not None:
        longest = validate_int(args.longest_seq, 0, '--longest-seq')
    else:
        longest = None

    # sanitize and set uniuqe header
    if args.print_statistics:
        print_stats = True
    else:
        print_stats = False

    if args.random_subsample is not None:
        random_subsample = validate_int(args.random_subsample, 0, '--random-subsample')
    else:
        random_subsample = None

    if longest is not None and shortest is not None:
        if longest < shortest:
            exit_error('--longest-seq must be longer than --shortest-seq')

    hp: Optional[Callable[[str], str]]
    if args.remove_comma_from_header:
        def hp(s: str) -> str:
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
    # return_list=True is passed below, so the return value is always a list of
    # [header, sequence] pairs -- cast so the length/subsample filters that
    # follow are type-checkable.
    data: list[list[str]] = cast(list, protfasta.read_fasta(args.filename,
                                expect_unique_header=expect_unique_header,
                                header_parser=hp,
                                duplicate_sequence_action=duplicate_sequence,
                                duplicate_record_action=duplicate_record,
                                invalid_sequence_action=invalid_sequence,
                                correction_dictionary = correction_dict,
                                return_list=True,
                                verbose=verb))



    # if length filters are requested
    if longest is not None:
        stdout('[INFO]: Filtering out sequences longer than %s'%(longest),silent)
        tmp = []
        for i in data:
            if len(i[1]) < longest:
                tmp.append(i)
        data = tmp

    if shortest is not None:
        stdout('[INFO]: Filtering out sequences shorter than %s'%(shortest), silent)
        tmp = []
        for i in data:
            if len(i[1]) > shortest:
                tmp.append(i)
        data = tmp

    if len(data) < 1:
        stdout('[INFO]: 0 sequences remain after filtering',silent)
        sys.exit(0)


    if random_subsample is not None:
        if len(data) < random_subsample:
            stdout('[INFO]: Cannot subsample as the requested number to subsample (%i) is more than\n        the total number of sequences (%i). Using all sequences'%(random_subsample, len(data)), silent)
        else:
            stdout('[INFO]: Subsampling %i sequences from the complete dataset (%i)'%(random_subsample, len(data)), silent)
        tmp = []
        x = list(range(0,len(data)))
        random.shuffle(x)
        idx = x[0:random_subsample]
        for position in idx:
            tmp.append(data[position])
        data = tmp



    if print_stats and not silent:
        print_statistical_summary(data)

    if args.no_outputfile:
        stdout('[INFO]: No outputfile requested ',silent)
    else:
        stdout('[INFO]: Writing new sequence file [%s]...'%(outfile),silent)
        protfasta.write_fasta(data, outfile, linelength=number_of_lines)

        
                   
            

