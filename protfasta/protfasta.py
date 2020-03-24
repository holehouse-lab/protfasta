"""
protfasta - A simple but robuts FASTA parser explicitly for protein sequences.

This file handles the main event!


.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other. 

"""

from . import io, utilities
from .protfasta_exceptions import ProtfastaException
    
def read_fasta(filename, 
               expect_unique_header=True,
               header_parser=None,
               duplicate_sequence_action='ignore',
               duplicate_record_action='fail',
               invalid_sequence_action='fail',
               return_list=False,               
               output_filename=None,
               correction_dictionary=None,
               verbose=False):

    """
    read_fasta is the main one of of only two user-facing functions associated with protfasta.

    It is designed as a catch-all function for reading in a FASTA file, performing santization,
    and returning a list or dictionary of sequences and their associated headers.

    There are a number of parameters which can be included, but as one might expect the simplest
    usage is just
    
    >>> x = read_fasta(filename)
    
    This will read in the file associated with filename and return a dictionary, where the keys
    are the FASTA file headers and the values are the amino acid sequences associated with each.

    Note that as of python 3.7 the order in which one adds items to a dictionary is guarenteed
    to be the order in which they're retrieved, so cycling through the resulting dictionary will
    in fact allow you to cycle through in order. 

    In addition to this simple usage, there are a number of keywords which are described in depth
    below and allow additional processing to be complete. 

    There is an order of options in which sanitization occurs:

        1. File is read in (and unique headers are tested for HERE)

        2. Check for duplicate records and respond appropriately (optional)

        3. Check for duplicate sequences and respond appropriately (optional)
    
        4. Invalid sequences dealt with (optional)

        5. Final set of sequences/headers written to a new FASTA file (optional)
    
        6. Dictionary/list returned to user.

    Understanding there is a specific order is important when considering what options to
    pass. If a set of options are incompatible, this will be caught before the file is read.
    
    Sequence corrections:

    
    Parameters
    ---------------

    expect_unique_header : boolean {True}
        Should the function expect each header to be unique? In general this is true for FASTA files, 
        but this is strictly not guarenteed. If this is set to True and a duplicate header is found
        then this means an error will be thrown. If it's set to false duplicate headers are dealt with,
        although for this to work return_list must also be set to True. Note that this won't happen
        automatically to avoid the scenario where you expect a dictionary to return and actually get
        a list


    header_parser : function {None}
        header_parser is a user-defined function that will be fed the FASTA header and whatever it returns
        will be used as the actual header as the files are parsed. This can be useful if you know your FASTA
        header has a consistent format that you want to take advantage of.

        A function provided here MUST:
            1. Take a single input argument (the header string)
            2. Return a single string

        When parsing this function the following test is applied
            a = header_parser('this test string should work')

        So ensure 

    duplicate_record_action : string {'ignore'}
        Selector that determines how to deal with duplicate entries. Note that duplicate records refers to
        entries in the fasta file where both the sequence and the header are identical. duplicate_record_action
        is only relevant keyword when expect_unique_header is False.

        Options are as follows:        
            'ignore'   : duplicate entries are allowed and ignored

            'fail'     : duplicate entries cause parsing to fail

            'remove'   : duplicate entries are removed, so there's only one copy of any duplicates
    

    duplicate_sequence_action : string {'ignore'}
        Selector that determines how to deal with duplicate sequences. This completely ignores the header
        and simply asks is two sequences are duplicated (or not). 


        Options are as follows:        
            'ignore'   : duplicate sequences are allowed and ignored

            'fail'     : duplicate sequences cause parsing to fail

            'remove'   : duplicate sequences are removed, guarenteeing that 
    

    invalid_sequence_action :  string {'fail'}
        Selector that determines how to deal with invalid sequences. Options are as follows:

            'ignore'   : invalid sequences are completely ignored

            'fail'     : invalid sequence cause parsing to fail

            'remove'  : invalid sequences are removed

            'convert'  : invalid sequences are converted to valid sequences (assuming the invalid

        By default it fails (which is the safest option).


    return_list : boolean {False}
        Flag that tells the function to return a list of 2-mer lists (where position 0 is the header
        and position 1 the sequence). If you have duplicate identical headers which you want to deal
        with, this is required.


    output_filename : stringe {None}
        If you are performing sanitization of the input file it is often useful to write out the 
        actual set of sequences you'll be analyzing, so you have a persistent copy of this data 
        for further analysis later on. If you provide a string to output filename it will cause
        a new FASTA file to be written with the final set of sequences returned.


    correction_dictionary : dict {None}
        protfasta can automatically correct non-standard amino acids to standard amino acids.
        This is useful if downstream analysis assumes fully standard amino acids. This is also 
        useful for removing '-'  from aligned sequences.
    
        The standard conversions used are:

        B -> N
        U -> C
        X -> G
        Z -> Q
        * -> <empty string>
        - -> <empty string>

        If alternative definitions are needed they can be passed via the correction_dictionary
        keyword. This should be a dictionary that maps characters to valid  amino acids (upper case).


    verbose : boolean {False}
        If set to True, protfasta will print out information as it works its way through reading and
        parsing FASTA files. This can be useful for diagnosis.

    
      
    """
    
    # first we sanity check all of the inputs provided. NOTE. If additional functionality is added, new
    # keywords MUST be sanity checked in this function
    io.check_inputs(expect_unique_header,
                    header_parser, 
                    duplicate_record_action,
                    duplicate_sequence_action,
                    invalid_sequence_action, 
                    return_list, 
                    output_filename,
                    verbose)

    # the actual file i/o happens here
    raw = io.internal_parse_fasta_file(filename, expect_unique_header=expect_unique_header, header_parser=header_parser, verbose=verbose)

    # first deal with duplicate records
    updated = _deal_with_duplicate_records(raw, duplicate_record_action, correction_dictionary)

    # deal with duplicate sequences
    updated = _deal_with_duplicate_sequences(updated, duplicate_sequence_action, verbose)

    # next decide how we deal with invalid amino acid sequences
    updated = _deal_with_invalid_sequences(updated, invalid_sequence_action, verbose)

    # finally, if we wanted to write the final set of sequences we're going to use...:
    if return_list:
        return updated
    else:
        return utilities.convert_list_to_dictionary(updated, verbose)
        return updated
        


# ------------------------------------------------------------------
#
def write_fasta(fasta_data, filename, linelength=60, verbose=False):
    """
    Simple function that takes a dictionary of key to sequence values
    and writes out a valid FASTA file. 

    No return type, but writes a file to disk according to the location
    defined by filename .


    Parameters
    -----------

    fasta_data : dictionary or list
        If a dictionary is passed then keys must be identifiers and the values are 
        amino acid sequences. If a list is passed it must be a 


    filename [string]
        Filename to write to. Should end with .fasta but this is not 
        enforced.

    linelength [int] {60}
        Length of line to be written for sequence (note this does
        not effect the header line. 60 is default used by Uniprot.
        If set to 0, None or False no line-length limit is used.


    """

    # This part of code means we can pass either a dictionary or a list of lists
    # in for write_fasta to deal with

    if type(fasta_data) == list:
        def get_sequence():
            return (entry[0], entry[1])

        # quick validate
        for i in fasta_data:
            if len(i)  != 2:
                raise ProtfastaException('While processing a list for write_fasta_file at least one of the elements was not a 2-position sublist:\n%s'%(str(i)))

    if type(fasta_data) == dict:
        def get_sequence():
            return (entry, fasta_data[entry])
        

    # override line length for sane input. N
    if linelength == False or linelength == None or linelength < 1:
        linelength = False
        
    else:        
        # cast linelength to an integer here as a soft type checking
        linelength = int(linelength)


    # open the file handler - all
    with open(filename,'w') as fh:

        # for each entry
        for entry in fasta_data:


            (header, seq) = get_sequence()

            # write the header line
            fh.write('>%s\n'%(header))
            
            # the $wrotenewline boolean is ONLY here to
            # avoid the unlikely scenario in which the last character
            # is also an integer number of linelength, such that you'd
            # get TWO spaces between sequences. This ensures there is ALWAYS
            # only one blank line between sequences
            for i in range(0,len(seq)):
                fh.write('%s'%seq[i])
                wrotenewline=False

                # if linelength is valid
                if linelength:

                    # if we reach an integer number of line-length write
                    # a newline character
                    if (i+1) % linelength == 0:
                        fh.write('\n')
                        wrotenewline=True
            
            if wrotenewline:
                fh.write('\n')
            else:
                fh.write('\n\n')




####################################################################################################
#
#
def _deal_with_invalid_sequences(raw, invalid_sequence_action='fail', verbose=False, correction_dictionary=None):

    # fail on an invalid sequence
    if invalid_sequence_action == 'fail':
        utilities.fail_on_invalid_sequences(raw)
        return raw

    # simply remove invalid sequences
    if invalid_sequence_action == 'remove':
        updated = utilities.remove_invalid_sequences(raw)
        
        if verbose:
            print('Removed %i of %i due to sequences with invalid characters' % (len(raw) - len(updated), len(raw)))

        return updated

    # convert invalid sequences
    if invalid_sequence_action == 'convert':
        (updated, count) = utilities.convert_invalid_sequences(raw, correction_dictionary)
        if verbose:
            print('Converted %i sequences to valid sequences'%(count))

        # note we then rescan in case there were still characters we couldn't deal with
        try:
            utilities.fail_on_invalid_sequences(updated)
        except ProtfastaException as e:
            raise ProtfastaException("******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s"%(str(e)))

        return updated

    # ignore invalid sequences
    if invalid_sequence_action == 'ignore':
        return raw
        
    raise ProtfastaException("Invalid option passed to the selector 'invalid_sequence_action': %s" %(invalid_sequence_action))



####################################################################################################
#
#
def _deal_with_duplicate_records(raw, duplicate_record_action='ignore', verbose=False):

    if duplicate_record_action == 'ignore':
        pass

    if duplicate_record_action == 'fail':
        utilities.fail_on_duplicates(raw)

    if duplicate_record_action == 'remove':
        updated = utilities.remove_duplicates(raw)
        if verbose:
            print('Removed %i of %i due to duplicate records ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw



####################################################################################################
#
#
def _deal_with_duplicate_sequences(raw, duplicate_sequence_action='ignore', verbose=False):


    if duplicate_sequence_action == 'ignore':
        pass

    if duplicate_sequence_action == 'fail':
        utilities.fail_on_duplicate_sequences(raw)
        
    if duplicate_sequence_action == 'remove':
        updated = utilities.remove_duplicate_sequences(raw)
        if verbose:
            print('Removed %i of %i due to duplicate sequences ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw



if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pass
