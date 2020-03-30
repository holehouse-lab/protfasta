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


from . import io as _io
from . import utilities as _utilities
from .protfasta_exceptions import ProtfastaException

    

####################################################################################################
#
#
def _deal_with_invalid_sequences(raw, invalid_sequence_action='fail', verbose=False, correction_dictionary=None):
    """
    Function that defines how to deal with invalid amino acids in a nest-sequence list.

    Parameters
    -----------

    raw : list
        List of lists, where each sub-list has two elements (0 = header, 1 = sequence)

    invalid_sequence_action : {'fail','remove', 'ignore', 'convert', 'convert-ignore'}
        Keyword that defines the action taken upon encountering an invalid amino acid residue.
        default = 'fail'.

    verbose : bool
        Defines if the function provides extensive output to STDOUT (True) or not (False)
        default = False

    correction_dictionary : dict
        Dictionary that defines how specific residues will be converted. By default, the following conversions
        are used:
        
            * ``B -> N``
            * ``U -> C``
            * ``X -> G``
            * ``Z -> Q``
            * ``* -> <empty string>``
            * ``- -> <empty string>``

        HOWEVER, a customizable dictionary can be passed

    Returns
    --------
        Returns an updated nested list of lists where the appropriate action has been taken to alter entries
        with invalid sequences.

    """

    # sequencescode on an invalid sequence
    if invalid_sequence_action == 'fail':
        _utilities.fail_on_invalid_sequences(raw)
        return raw

    # simply remove invalid sequences
    if invalid_sequence_action == 'remove':
        updated = _utilities.remove_invalid_sequences(raw)
        
        if verbose:
            print('[INFO]: Removed %i of %i due to sequences with invalid characters' % (len(raw) - len(updated), len(raw)))

        return updated

    # convert invalid sequences
    if invalid_sequence_action == 'convert' or invalid_sequence_action == 'convert-ignore':
        (updated, count) = _utilities.convert_invalid_sequences(raw, correction_dictionary)
        if verbose:
            print('[INFO]: Converted %i sequences to valid sequences'%(count))

        # note we then rescan in case there were still characters we couldn't deal with
        if invalid_sequence_action == 'convert':
            try:
                _utilities.fail_on_invalid_sequences(updated)
            except ProtfastaException as e:
                raise ProtfastaException("\n\n******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s"%(str(e)))

        return updated

    # ignore invalid sequences
    if invalid_sequence_action == 'ignore':
        return raw
        
    raise ProtfastaException("Invalid option passed to the selector 'invalid_sequence_action': %s" %(invalid_sequence_action))



####################################################################################################
#
#
def _deal_with_duplicate_records(raw, duplicate_record_action='ignore', verbose=False):
    """
    Function that defines how to deal with duplicate FASTA records

    Parameters
    -----------

    raw : list
        List of lists, where each sub-list has two elements (0 = header, 1 = sequence)

    invalid_sequence_action : {'fail','remove', 'ignore'}
        Keyword that defines the action taken upon encountering a duplicate record.
        default = 'fail'.

    verbose : bool
        Defines if the function provides extensive output to STDOUT (True) or not (False)
        default = False

    Returns
    --------
        Returns an updated nested list of lists where the appropriate action has been taken 
        to alter duplicate records

    """

    if duplicate_record_action == 'ignore':
        pass

    if duplicate_record_action == 'fail':
        _utilities.fail_on_duplicates(raw)

    if duplicate_record_action == 'remove':
        updated = _utilities.remove_duplicates(raw)
        if verbose:
            print('[INFO]: Removed %i of %i due to duplicate records ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw



####################################################################################################
#
#
def _deal_with_duplicate_sequences(raw, duplicate_sequence_action='ignore', verbose=False):
    """
    Function that defines how to deal with duplicate sequences

    Parameters
    -----------

    raw : list
        List of lists, where each sub-list has two elements (0 = header, 1 = sequence)

    invalid_sequence_action : {'fail','remove', 'ignore'}
        Keyword that defines the action taken upon encountering a duplicate sequence.
        default = 'ignore'.

    verbose : bool
        Defines if the function provides extensive output to STDOUT (True) or not (False)
        default = False

    Returns
    --------
        Returns an updated nested list of lists where the appropriate action has been taken 
        to alter duplicate sequences

    """

    if duplicate_sequence_action == 'ignore':
        pass

    if duplicate_sequence_action == 'fail':
        _utilities.fail_on_duplicate_sequences(raw)
        
    if duplicate_sequence_action == 'remove':
        updated = _utilities.remove_duplicate_sequences(raw)
        if verbose:
            print('[INFO]: Removed %i of %i due to duplicate sequences ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw

