"""
protfasta - A simple but robuts FASTA parser explicitly for protein sequences.

This file handles the input and output, including validation of input arguments

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other. 

"""

from .protfasta_exceptions import ProtfastaException


def check_inputs(expect_unique_header, 
                 header_parser, 
                 duplicate_record_action, 
                 duplicate_sequence_action,
                 invalid_sequence_action, 
                 return_list, 
                 output_filename, 
                 verbose):
    """
    Function that performs sanity validation for all input arguments. If arguments do not match the expected
    behaviour then this function will throw an exception.

    If new functionality is included, input must be santitized in this function!

    Parameters
    ------------
    expect unique : ?
        Checks it's a bool

    header_parser : ?
        Checks it's a callable function that returns a single string (and takes a single
        string as the input argument)

    duplicate_record_action : ?
        Checks its a string that matches one of a specific set of keywords

    duplicate_sequence_action : ?
        Checks its a string that matches one of a specific set of keywords

    invalid_sequence : ?
        Checks it's a string that matches a specific keyword

    return_list : ?
        Checks it's a bool

    output_filename : ?
        Checks it's a string

    verbose : ?
        Checks it's a bool

    Returns
    ---------

        No return - stateless function that checks things are OK
    
    """

    # check the expect_unique_header keyword
    if type(expect_unique_header) != bool:
        raise ProtfastaException("keyword 'expect_unique_header' must be a boolean")

    # validate the header_parser 
    if header_parser is not None:
        if not callable(header_parser):
            raise ProtfastaException("keyword 'header_parser' must be a function [tested with callable()]")
        
        try:
            a = header_parser('this test string should work')
            if type(a) != str:
                raise Exception
        except Exception:
            raise ProtfastaException('Something went wrong when testing the header_parser function.\nEnsure that the test example works and that the function returns a string [type str]')

    # check the duplicates_record_action 
    if duplicate_record_action not in ['ignore','fail','remove']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail','remove'")

    # check the duplicates_record_action 
    if duplicate_sequence_action not in ['ignore','fail','remove']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail', 'remove'")


    # check the invalid_sequence 
    if invalid_sequence_action not in ['ignore','fail','remove','convert','convert-ignore']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail','remove','convert','convert-ignore'")

    # check the return_list
    if type(return_list) != bool:
        raise ProtfastaException("keyword 'verbose' must be a boolean")

    # check the return_list
    if output_filename is not None:
        if type(output_filename) != str:
            raise ProtfastaException("keyword 'output_filename' must be a string")

    if type(verbose) != bool:
        raise ProtfastaException("keyword 'verbose' must be a boolean")

    if duplicate_record_action is 'ignore':
        if expect_unique_header is True:
            raise ProtfastaException('Cannot expect unique headers and ignore duplicate records')





####################################################################################################
#
#    
def internal_parse_fasta_file(filename, expect_unique_header=True, header_parser=None, verbose=False):
    """
    Base level FASTA file parser. Header lines must begin with a ">" and be a single line. 
    No other requirements are necessary.

    This is not the function that we expect users to use, but if they wanted to they could.

    Parameters
    ------------

    filename : string
        String representing the absolute or relative path of a FASTA file.

    expect_unique_header : boolean {True}
        Should the function expect each header to be unique? In general this is true for FASTA files, 
        but this is strictly not guarenteed. 

    header_parser : function {None}
        header_parser is a user-defined function that will be fed the FASTA header and whatever it returns
        will be used as the actual header as the files are parsed. This can be useful if you know your FASTA
        header has a consistent format that you want to take advantage of.

        A function provided here MUST:
            1. Take a single input argument (the header string)
            2. Return a single string

    verbose : boolean {False}
        If set to True, protfasta will print out information as it works its way through reading and
        parsing FASTA files. This can be useful for diagnosis.

    Returns
    ----------
    list of lists
        Retrns a list of lists, where each sublist contains a FASTA header and a sequence

    """
        
    # read in the file...
    try:
        with open(filename,'r') as fh:
            content = fh.readlines()
    except FileNotFoundError:
        raise ProtfastaException('Unable to find file: %s'%(filename))
    
    if verbose:
        print('[INFO]: Read in file with %i lines'%(len(content)))

    # note, we'll keep the ability to directly parse dictionaries
    return _parse_fasta_all(content, 'list', expect_unique_header=expect_unique_header, header_parser=header_parser, verbose=verbose)
    


####################################################################################################
#
#    
def _parse_fasta_all(content, mode, expect_unique_header=True, header_parser=None, verbose=False):
    """
    Internal function that actually performs parsing. Note that at present, the only code within
    protfasta calls ``_parse_fasta_all`` with mode = 'list'.
    

    Parameters
    -----------
    content : list
        A list generated by reading in a file (i.e. each line corresponds to a consecutive line 
        in a file).
   
    mode: {'list', 'dict'}
        Keyword selector that defines the return type

    expect_unique_header : boolean {True}
        Should the function expect each header to be unique? In general this is true for FASTA files, 
        but this is strictly not guarenteed. 

    header_parser : function {None}
        header_parser is a user-defined function that will be fed the FASTA header and whatever it returns
        will be used as the actual header as the files are parsed. This can be useful if you know your FASTA
        header has a consistent format that you want to take advantage of.

        A function provided here MUST:
            1. Take a single input argument (the header string)
            2. Return a single string

    verbose : boolean {False}
        If set to True, protfasta will print out information as it works its way through reading and
        parsing FASTA files. This can be useful for diagnosis.


    Returns 
    ----------
    list or dict 
        Depending on the ``mode`` selector either a list or a dictionary of sequences is returned. In the case
        of a dictionary, keys are FASTA headers and values are sequences. In the case of a lists, each element
        is a two-position sub-list, where element 0 is the header and element 1 is the sequence.    

    """

    
    # --------------------------------------------------------------
    ## Define local functions if we want to return a dictionary
    if mode == 'dict':
        def check_duplicate():
            if header in return_data:
                raise ProtfastaException('Found non-unique FASTA header [%s]'%(header))

        def update():
            return_data[header] = seq.upper()
        return_data={}

    # --------------------------------------------------------------
    # define local functions if we want to return a list
    elif mode == 'list':
        ## Define local functions if we want to return a dictionary        
        def check_duplicate():
            if expect_unique_header:
                if header in all_headers:
                    raise ProtfastaException('Found non-unique FASTA header [%s]'%(header))

        def update():
            return_data.append([header,seq.upper()])
            all_headers.append(header)

        return_data=[]
        all_headers=[]
    
    ## START OF PARSING FUNCTION    
    seq=''
    header=''

    for line in content:
        
        sline=line.strip()

        # if empty line just skip...
        if len(sline) == 0:
            continue

        # if  first non-whitespace character is a '>'
        if sline[0] == '>':

            # get the current header line
            h = sline[1:]

            
            # if we'd previously had a sequence assigned, means we have just started a 
            # 'new' sequence
            if len(seq) > 0:
                
                # see if duplicate header can be found - will raise exception if duplicate found 
                check_duplicate()
                update()
                                    
            # reset the header and sequence
            if header_parser:
                header = header_parser(h)
            else:
                header = h
            seq=''
            
        else:            
            # we're on a line that is neither empty nor started with
            # a '>' so it's treated as a sequence line
            seq = seq + sline

    # if we exit with a sequence in toe, there's one final sequence to
    # add...        
    if len(seq) > 0:
        check_duplicate()
        update()

    if verbose:
        print('[INFO]: Parsed file to recover %i sequences' %(len(return_data)))


    return return_data



