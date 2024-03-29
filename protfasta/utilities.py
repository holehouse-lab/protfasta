from .protfasta_exceptions import ProtfastaException
from ._configs import STANDARD_CONVERSION, STANDARD_CONVERSION_WITH_GAP, STANDARD_AAS, STANDARD_AAS_WITH_GAP

# To do: Need to fully annotate the docstring for these
# internal functions
#






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
def convert_to_valid(seq, correction_dictionary=None, alignment=False):
    """
    Function that converts non-standard amino acid residues to standard ones.
    Specifically:

    B -> N
    U -> C
    X -> G
    Z -> Q
    ' ' -> <empty string> (i.e. an empty space)
    * -> <empty string>
    - -> <empty string> (ONLY if alignment is False)

    By default this alters the underlying sequence. If you wish to return a copy
    of the altered sequence instead set copy=True. Otherwise the underlying sequence
    is altered 


    Parameters
    ---------------
    seq : string 
        Amino acid sequence

    alignment : bool
        If True then dashes are treated as OK and ignored.
    
    Returns
    ---------------
    string
        Returns the same string with non-standard residues converted according 
        
    """

    if correction_dictionary:
        converter = correction_dictionary
    else:
        if alignment:
            converter = STANDARD_CONVERSION_WITH_GAP
        else:
            converter = STANDARD_CONVERSION

    # cycle over each key in the converter dictionary and replace all values in
    # the sequence with the corresponding converted one 
    for i in converter:
        seq = seq.replace('%s'%(i), converter[i])

    return seq



####################################################################################################
#
#    
def check_sequence_is_valid(seq, alignment=False):
    """

    Parameters
    --------------
    seq : str
        Amino acid sequence

    alignment : bool
        Flag that defines if this alignment sequence rules
        should be applied or not.

    Returns
    ------------
    Tuple
        Returns a tuple of size 2 where 
        
        element 0 is a boolean (True or False) that flags if the 
        sequence was valid (True) or not (False).

        element 1 is a value that will return as the invalid amino acid
        (if sequence is invalid) OR if it's a valid sequence will be 0
    """

    if alignment == True:
        valid_AA_list = STANDARD_AAS_WITH_GAP
    else:
        valid_AA_list = STANDARD_AAS
    
    s = list(set(seq))
    for i in s:
        if i not in valid_AA_list:
            return (False, i)

    return (True, 0)



####################################################################################################
#
#    
def convert_invalid_sequences(dataset, correction_dictionary=None, alignment=False):

    count = 0
    for idx in range(0,len(dataset)):
        s = dataset[idx][1]
        dataset[idx][1] = convert_to_valid(s, correction_dictionary, alignment)
        
        if s != dataset[idx][1]:
            count = count + 1
        
    return (dataset, count)



####################################################################################################
#
#    
def remove_invalid_sequences(dataset, alignment=False):

    # this single line iterates through the dataset, and for each sequence only
    # only adds to to the growing new list IF that sequence is valid
    return [element for element in dataset if check_sequence_is_valid(element[1], alignment)[0]]
     

       
####################################################################################################
#
#    
def fail_on_invalid_sequences(dataset, alignment=False):

    # cycle over each entry and fail if a sequence is invald!
    for entry in dataset:

        (status, info) = check_sequence_is_valid(entry[1], alignment)

        if status is not True:
            raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (info, entry[0], entry[1]))



####################################################################################################
#
#    
def convert_list_to_dictionary(raw_list, verbose=False):

    if verbose:
        return_dict={}
        warning_count=0
        for entry in raw_list:
            if entry[0] in return_dict:
                warning_count=warning_count+1
                print('[WARNING]: Overwriting entry [count = %i]'%(warning_count))
            return_dict[entry[0]] = entry[1]
        if warning_count > 0:
            print('[INFO] If you want to avoid overwriting duplicate headers set return_list=True')
        else:
            print('[INFO]: All processed sequences uniquely added to the returning dictionary')
            
    else:
        return_dict={}
        for entry in raw_list:
            return_dict[entry[0]] = entry[1]

    return return_dict



####################################################################################################
#
#    
def fail_on_duplicates(dataset):
    lookup={}
    for entry in dataset:
        if entry[0] not in lookup:
            lookup[entry[0]] = entry[1]
        else:
            if lookup[entry[0]] == entry[1]:
                raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (entry[0], entry[1]))



####################################################################################################
#
#    
def remove_duplicates(dataset):


    # note - the lookup dictionary is ONLY used to keep track of what we have/have not
    # yet seen and is not returned. The updated list will be returned
    lookup={}
    updated = []

    # for each entry in the dataset
    for entry in dataset:

        # if we've never seen this header - add it to the lookup dictionary
        # as an entry in a list. This means multiple entries can have the same header (which
        # we're saying is OK)
        if entry[0] not in lookup:
            lookup[entry[0]] = [entry[1]]
            updated.append(entry)

        # however, if we HAVE seen this header before..
        else:
            found_dupe = False
            
            # for each sequence associated wth that header
            # ask if we found a duplicate 
            for d in lookup[entry[0]]:
                if d == entry[1]:
                    found_dupe=True

            # if we found no duplicates add
            if not found_dupe:
                lookup[entry[0]].append(entry[1])
                updated.append(entry)

    return updated

        

####################################################################################################
#
#    
def fail_on_duplicate_sequences(dataset):    
    seq_to_header={}
    for entry in dataset:
        if entry[1] in seq_to_header:
            raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (seq_to_header[entry[1]], entry[0]))
        seq_to_header[entry[1]] = entry[0]



####################################################################################################
#
#    
def remove_duplicate_sequences(dataset):    

    lookup = set([])
    updated=[]

    for entry in dataset:

        # if the sequence is in the lookup
        # set
        if entry[1] in lookup:
            continue

        # else its not in the lookup set so 
        # add it and then add to the updated list
        else:
            lookup.add(entry[1])
            updated.append(entry)
    return updated
            

    
                

                
