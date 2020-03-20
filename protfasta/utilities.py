from .protfasta_exceptions import ProtfastaException




STANDARD_CONVERSION = {'B':'N',
                       'U':'C',
                       'X':'G',
                       'Z':'Q',
                       '*':'',
                       '-':''}
STANDARD_AAS = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']



####################################################################################################
#
#    
def build_custom_dictionary(additional_dictionary):
    """

    """

    final_dict = {}
    for i in STANDARD_CONVERSION:        
        final_dict[i] = STANDARD_CONVERSION[i]

    # this deliberatly overwirtes the standard conversion
    # so a passed additional dictionary takes precedence! 
    for i in additional_dictionary:
        final_dict[i] = additional_dictionary[i]
        
    return final_dict



####################################################################################################
#
#    
def convert_to_valid(seq, additional_dictionary=None):
    """
    Function that converts non-standard amino acid residues to standard ones.
    Specifically:

    B -> N
    U -> C
    X -> G
    Z -> Q
    * -> <empty string>
    - -> <empty string>

    By default this alters the underlying sequence. If you wish to return a copy
    of the altered sequence instead set copy=True. Otherwise the underlying sequence
    is altered 


    Parameters
    ---------------
    seq : string 
        Amino acid sequence
    
    Returns
    ---------------
    string
        Returns the same string with non-standard residues converted according 
        
    """

    if additional_dictionary:
        converter = additional_dictionary
    else:
        converter = STANDARD_CONVERSION

    # cycle over each key in the converter dictionary and replace all values in
    # the sequence with the corresponding converted one 
    if i in converter:
        seq = seq.replace('%s'%(i), converter[i])

    return seq


def check_sequence_is_valid(seq):
    s = list(set(seq))
    for i in s:
        if i not in STANDARD_AAS:
            return False
    return True



def convert_sequences(dataset):

    for idx in range(0,len(dataset)):
        s = dataset[idx][1]
        dataset[idx][1] = convert_to_valid(s)
        
    return dataset


def ecluded_invalid_sequence(dataset):

    # this single line iterates through the dataset, and for each sequence only
    # only adds to to the growing new list IF that sequence is valid
    return [element for element in dataset if check_sequence_is_valid(element[1])]
            

def fail_on_invalid_sequence(dataset):

    # cycle over each entry and fail if a sequence is invald!
    for entry in dataset:
        if check_sequence_is_valid(entry[1]) is not True:
            ProtfastaException('Failed on invalid sequence:\n%s\n%s\n' % (entry[0],entry[1]))

            

def convert_list_to_dictionary(raw_list):
    return_dict={}
    for entry in raw_list:
        return_dict[entry[0]] = entry[1]

    return return_dict
