"""
protfasta.py
A hyper robuts FASTA parser explicitly for protein sequences

Handles the primary functions
"""

from . import io, utilities
from .protfasta_exceptions import ProtfastaException

        
        

        

    
    

    
def read_fasta(filename, 
               expect_unique=True,
               header_parser=None,
               invalid_sequence='fail'
               return_list=False,
               verbose=False,
               write_final_sequences=False):

    """


    Paramaeters
    ---------------

    invalid_sequence :  string
        Selector that determines how to deal with invalid sequences. Options are as follows:

            'ignore'   : invalid sequences are completely ignored

            'fail'     : invalid sequence cause parsing to fail

            'exclude'  : invalid sequences are excluded

            'convert'  : invalid sequences are converted to valid sequences (assuming the invalid

    

    """
               
    
    raw = io.parse_fasta_file(filename, expect_unique=expect_unique, header_parser=header_parser, return_list=return_list, verbose=verbose)

    #  next deal with invalid amino acids
    updated = _deal_invalid_sequences(raw, invalid_sequence, verbose)



    # finally, if we wanted to write the final set of sequences we're going to use...:
    


    if return_list:
        return updated
    else:
        return utilities.convert_list_to_dictionary(updated)
        


# ------------------------------------------------------------------
#
def write_fasta_file(fasta_data, filename, linelength=60, verbose):
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

    if type(fasta_data) == list:
        def get_sequence():
            return (entry[0], entry[1]))

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
def _deal_invalid_sequences(raw, invalid_sequence='fail', verbose=False):

    # make sure case issues don't arise
    invalid_sequence = invalid_sequence.lower()


    # fail on an invalid sequence
    if invalid_sequence == 'fail':
        utilities.fail_on_invalid_sequence(raw)
        return raw

    # simply exclude invalid sequences
    if invalid_sequence == 'exclude':
        updated = utilities.exclude_invalid_sequence(raw)
        
        if verbose:
            print('Excluded %i of %i due to sequences with invalid characters' % (len(raw) - len(updated), len(raw)))

        return raw

    # convert invalid sequences
    if invalid_sequence == 'convert':
        (updated, count) = utilities.convert_invalid_sequence(raw)
        if verbose:
            print('Updated %i sequences')

        # note we then rescan in case there were still characters we couldn't deal with
        try:
            utilities.fail_on_invalid_sequence(updated)
        except ProtfastaException as e:
            print('******* Despite fixing fixable errors, additional problems remain with the sequence*********')
            raise e

        return updated

    # ignore invalid sequences
    if invalid_sequence == 'ignore':
        return raw
        
    raise ProtfastaException("Invalid option passed to the selector 'invalid_sequence': %s" %(invalid_sequence))


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pass
