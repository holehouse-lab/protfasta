"""
"""
# Add imports here
from protfasta import utilities as _utilities
from protfasta import io as _io
from protfasta._configs import STANDARD_AAS, STANDARD_CONVERSION
from protfasta import protfasta as _protfasta
import os
from protfasta.protfasta_exceptions import ProtfastaException

# Handle versioning with versioningit
from ._version import __version__


_ROOT = os.path.abspath(os.path.dirname(__file__))
def _get_data(path):    
    """
    This function is used for getting the absolute path for loading
    test data. Pay no attention to the code behind the curtain.
    
    """
    return os.path.join(_ROOT, 'data', path)


def read_fasta(filename, 
               expect_unique_header=True,
               header_parser=None,
               check_header_parser=True,
               duplicate_sequence_action='ignore',
               duplicate_record_action='fail',
               invalid_sequence_action='fail',
               alignment=False,
               return_list=False,               
               output_filename=None,
               correction_dictionary=None,
               verbose=False):    

    """
    ``read_fasta`` is the main one of of only two user-facing functions associated with **protfasta**. 
    It is designed as a catch-all function for reading in a FASTA file, performing sanitization, 
    and returning a list or dictionary of sequences and their associated headers.

    There are a number of parameters which can be included, but as one might expect the simplest
    usage is just
    
    >>> x = read_fasta(filename)
    
    This will read in the file associated with filename and return a dictionary, where the keys
    are the FASTA file headers and the values are the amino acid sequences associated with each.
    
    Note that as of python 3.7 the order in which one adds items to a dictionary is guaranteed
    to be the order in which they're retrieved, so cycling through the resulting dictionary will
    in fact allow you to cycle through in order. 

    In addition to this simple usage, there are a number of keywords which are described in depth
    below and allow additional processing to be complete. 

    There is an order of options in which sanitization occurs:
    
    1. File is read in, custom headers are parsed, and unique headers are tested (if ``expect_unique = True``)
    

    2. Check for duplicate records and respond appropriately (**optional**)
    
    3. Check for duplicate sequences and respond appropriately (**optional**)
    
    4. Invalid sequences dealt with (**optional**)
    
    5. Final set of sequences/headers written to a new FASTA file (**optional**)
    
    6. Dictionary/list returned to user.
    
    Understanding there is a specific order is important when considering what options to
    pass. If a set of options are incompatible, this will be caught before the file is read.

    Parameters
    ----------

    expect_unique_header : bool
        [**Default = True**] Should the function expect each header to be unique? In general this is true for FASTA files, 
        but this is strictly not guarenteed. If this is set to True and a duplicate header is found
        then this means an error will be thrown. If it's set to false duplicate headers are dealt with,
        although for this to work ``return_list`` must also be set to True. Note that this won't happen
        automatically to avoid the scenario where you expect a dictionary to return and actually get
        a list. 

    header_parser : function
        [**Default = None**] ``header_parser`` allows a user-defined function that will be fed the FASTA header and 
        whatever it returns will be used as the actual header as the files are parsed. This can be useful if you 
        know your FASTA header has a consistent format that you want to take advantage of. A function provided here MUST        
        **(1)** Take a single input argument (the header string) and **(2)** Return a single string.
        When parsing this function the following test is applied, unless ``check_header_parser`` is set to false.
            >>> return_string = header_parser('this test string should work')
        Where ``return_string`` is tested to be a string. The function will show an exception if this test fails and ``check_header_parser`` is set to true.

    check_header_parser : bool
        [**Default = True**] Flag which - if set to false - will not test if the header_parser function returns a valid string. 
        This may lead to unexpected header values if the passed header_parser function is not well defined.
        
    duplicate_record_action : ``'ignore'``, ``'fail'``, ``'remove'``
        [**Default = 'fail'**] Selector that determines how to deal with duplicate entries. Note that duplicate records refers to
        entries in the fasta file where both the sequence and the header are identical. duplicate_record_action
        is only relevant keyword when expect_unique_header is False.
        Options are as follows:        
            * ``ignore``  - duplicate entries are allowed and ignored

            * ``fail``    - duplicate entries cause parsing to fail and throw an exception
  
            * ``remove`` - duplicate entries are removed, so there's only one copy of any duplicates     
    
    duplicate_sequence_action : ``'ignore'``, ``'fail'``, ``'remove'``
        [**Default = 'ignore'**] Selector that determines how to deal with duplicate sequences. This completely ignores the header
        and simply asks is two sequences are duplicated (or not). 
            * ``ignore``  - duplicate sequences are allowed and ignored

            * ``fail``    - duplicate sequences cause parsing to fail and throw an exception
  
            * ``remove`` - duplicate sequences are removed, so there's only one copy of any duplicates (1st instance kept)     
    
    invalid_sequence_action : ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``, ``'convert-ignore', ``'convert-remove'``
        [**Default = 'fail'**] Selector that determines how to deal with invalid sequences. If ``convert`` or ``convert-ignore`` are chosen, then conversion is completed with either the standard conversion table (shown under the ``correction_dictionary`` documentation) or with a custom conversion dictionary passed to ``correction_dictionary``. 
        Options are as follows: 
            * ``ignore``  - invalid sequences are completely ignored

            * ``fail``    - invalid sequence cause parsing to fail and throw an exception
  
            * ``remove`` - invalid sequences are removed

            * ``convert`` - invalid sequences are convert

            * ``convert-ignore`` - invalid sequences are converted to valid sequences and any remaining invalid residues are ignored

            * ``convert-remove`` - invalid sequences are converted to valid sequences where possible, and any remaining sequences with invalid residues are removed

    alignment : bool
        [**Default = False**] Flag which - if set to true - the Fasta file is treated as containing alignments (with dashes) such that '-' characters are not
        treated as invalid or converted. Works in concert with other flags.  
    
    return_list : bool 
        [**Default = False**] Flag that tells the function to return a list of 2-mer lists (where position 0 is the header
        and position 1 the sequence). If you have duplicate identical headers which you want to deal with, this is required.

    output_filename : string 
        [**Default = None**] If you are performing sanitization of the input file it is often useful to write out the 
        actual set of sequences you'll be analyzing, so you have a persistent copy of this data 
        for further analysis later on. If you provide a string to output filename it will cause
        a new FASTA file to be written with the final set of sequences returned.

    correction_dictionary : dict
        [**Default = None**] **protfasta** can automatically correct non-standard amino acids to standard amino acids using the
        ``invalid_sequence`` keyword. This is useful if downstream analysis assumes/requires fully standard amino acids. 
        This is also useful for removing '-'  from aligned sequences. The standard conversions used are:
        
            * ``B   -> N``
            * ``U   -> C``
            * ``X   -> G``
            * ``Z   -> Q``
            * ``" " -> <empty string>`` (i.e. a whitespace character)
            * ``*   -> <empty string>``
            * ``-   -> <empty string>``
        However, if alternative definitions are needed they can be passed via the ``correction_dictionary`` keyword.
        The ``correction_dictionary`` should be a dictionary that maps sequences characters to some other character (ideally
        valid amino acid characters). In principle this could be used to perform arbitrary coarse-graining if a sequence...
        
    verbose : bool 
        [**Default = False**] If set to True, **protfasta** will print out information as it works its way through reading and
        parsing FASTA files. This can be useful for diagnosis.

    Returns
    ----------        
        Return type is *list* or *dict*
        If ``return_list`` is set to ``True`` then the function returns a list of lists. In each sublist contains two elements, where the first is the FASTA record header and the second the sequence. The order of FASTA records will match the order they were read in from the FASTA file. If ``return_list`` is ``False`` then the function returns a dictionary where the keys are the FASTA record heades and the values are the sequences. NOTE the order of keys will match the order that the FASTA file was read in IF the Python version is 3.7 or higher.

    """

    # first we sanity check all of the inputs provided. NOTE. If additional functionality is added, new
    # keywords MUST be sanity checked in this function
    _io.check_inputs(expect_unique_header,
                     header_parser, 
                     check_header_parser,
                     duplicate_record_action,
                     duplicate_sequence_action,
                     invalid_sequence_action, 
                     alignment,
                     return_list, 
                     output_filename,
                     verbose,
                     correction_dictionary)
    

    # the actual file i/o happens here
    raw = _io.internal_parse_fasta_file(filename, expect_unique_header=expect_unique_header, header_parser=header_parser, verbose=verbose)

    # first deal with duplicate records
    updated = _protfasta._deal_with_duplicate_records(raw, duplicate_record_action, verbose)

    # deal with duplicate sequences
    updated = _protfasta._deal_with_duplicate_sequences(updated, duplicate_sequence_action, verbose)

    # next decide how we deal with invalid amino acid sequences

    ##
    ## If we're using the convert-remove action...
    if invalid_sequence_action == 'convert-remove':

        # first run a convert ignore
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          'convert-ignore',
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=correction_dictionary)

        # then a remove on those that are left
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          'remove',
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=False)

    else:
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          invalid_sequence_action,
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=correction_dictionary)
        

    # If we wanted to write the final set of sequences we're going to use...:
    if output_filename:
        write_fasta(updated, output_filename)

    # if we asked for a list...
    if return_list is True:
        pass    
    else:
        updated = _utilities.convert_list_to_dictionary(updated, verbose)

        
    return updated
        


# ------------------------------------------------------------------
#
def write_fasta(fasta_data, filename, linelength=60, verbose=False, append_to_fasta=False):
    """
    Simple function that takes a dictionary of key to sequence values
    and writes out a valid FASTA file. No return type, but writes a file 
    to disk according to the location defined by the variable ``filename``.
    
    Parameters
    -----------
    fasta_data : dict or list
        If a dictionary is passed then keys must be identifiers and the values are 
        amino acid sequences. If a list is passed it must be a list where each element
        contains two sub-elements, a header, and a sequence.

    filename: string
        Filename to write to. Should end with `.fasta` or `.fa` but this is not 
        enforced.

    linelength : int
        [**Default = 60**] Length of line to be written for sequence (note this does
        not effect the header line. 60 is default used by UniProt. If set to 0, None or 
        False no line-length limit is used. Note ``linelength`` must be > 5.

    append_to_fasta : bool
        Whether to append to a fasta file that already exists. If this is set to True,
        if the file does not exist, protfasta will create a new file. However, if the
        file does exist, protfasta will simply append additional fasta entries to the
        existing file. 
        Default=False

    Returns 
    ----------
    None 
        No return value is provided but a new FASTA file is written to disk

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
        # cast linelength to an integer here as a soft type checking, and
        # if it's shorter than 5 then reset to 5 
        linelength = int(linelength)
        if linelength < 5:
            linelength = 5

    # set the 'mode' for open. If append_to_file==False, use 'w' and overwrite
    # existing .fasta file. Otherwise use 'a' and add to existing file if it exists.
    if append_to_fasta==False:
        open_mode='w'
    else:
        open_mode='a'
    
    # open the file handle.
    with open(filename, open_mode) as fh:

        # for each entry
        for entry in fasta_data:

            (header, seq) = get_sequence()
            if len(seq) < 1:
                raise ProtfastaException('Seqence associated with [%s] is empty'%(header))
                

            # write the header line
            fh.write('>%s\n'%(header))
            
            # the $wrotenewline bool is ONLY here to
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
    
