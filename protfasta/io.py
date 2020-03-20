"""
io.py



"""


from .protfasta_exceptions import ProtfastaException


####################################################################################################
#
#    
def parse_fasta_file(filename, expect_unique=True, header_parser=None, return_list=False, verbose=False):
    """
    Base level FASTA file parser. Header lines must begin with a ">" and be a single line. 
    No other requirements are necessary.

    Parameters
    ------------

    filename : string
        String representing the absolute or relative path of a 


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
        raise ProtFastaException('Unable to find file: %s'%(filename))
    
    if verbose:
        print('Read in file with %i lines'%(len(content)))

    # note, we'll keep the ability to directly parse dictionaries
    return _parse_fasta_all(content, 'list', header_parser, expect_unique)
    


####################################################################################################
#
#    
def _parse_fasta_all(content, mode, header_parser=None, expect_unique=True):
    
    # --------------------------------------------------------------
    ## Define local functions if we want to return a dictionary
    if mode == 'dict':
        def check_duplicate():
            if header in return_data:
                raise ProtFastaException('Found non-unique FASTA header [%s]'%(header))

        def update():
            return_data[header] = seq.upper()
        return_data={}

    # --------------------------------------------------------------
    # define local functions if we want to return a list
    elif mode == 'list':
        ## Define local functions if we want to return a dictionary        
        def check_duplicate():
            if expect_unique:
                if header in all_headers:
                    raise ProtFastaException('Found non-unique FASTA header [%s]'%(header))

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

    return return_data



# ------------------------------------------------------------------
#
def write_fasta_file(sequence_dict, filename, linelength=60):
    """
    Simple function that takes a dictionary of key to sequence values
    and writes out a valid FASTA file. 

    No return type, but writes a file to disk according to the location
    defined by filename .

    -----------------------------------------------------------------

    sequence_dict [dict]
    A dictionary for which keys are identifiers and the values are 
    amino acid sequences. When writing a file the '>' character is 
    automatically added to the header.

    filename [string]
    Filename to write to. Should end with .fasta but this is not 
    enforced.

    linelength [int] {60}
    Length of line to be written for sequence (note this does
    not effect the header line. 60 is default used by Uniprot.
    If set to 0, None or False no line-length limit is used.


    """

    # override line length for sane input. N
    if linelength == False or linelength == None or linelength < 1:
        linelength = False
        
    else:        
        # cast linelength to an integer here as a soft type checking
        linelength = int(linelength)


    # open the file handler - all
    with open(filename,'w') as fh:

        # for each entry
        for key in sequence_dict:

            # extract out the sequence
            seq = sequence_dict[key]

            # write the header line
            fh.write('>%s\n'%(key))

            
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
