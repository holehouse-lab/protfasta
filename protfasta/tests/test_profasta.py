"""
Unit and regression test for the protfasta package.
"""

# Import package, test suite, and other packages as needed
import protfasta
from protfasta.protfasta_exceptions import ProtfastaException
import pytest
import sys

def test_protfasta_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "protfasta" in sys.modules

test_data={'test1':['sp|O00401|WASL_HUMAN Neural Wiskott-Aldrich syndrome protein OS=Homo sapiens OX=9606 GN=WASL PE=1 SV=2','MSSVQQQPPPPRRVTNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKKCSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALNFANEEEAKKFRKAVTDLLGRRQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVNNISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGWDPNTGFDLNNLDPELKNLFDMCGISEAQLKDRETSKVIYDFIEKTGGVEAVKNELRRQAPPPPPPSRGGPPPPPPPPHNSGPPPPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPPPPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLKKVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSKAIHSSDEDEDEDDEEDFEDDDEWED']}
def test_read_fasta_standard():

    
    test_data_dir = protfasta._get_data('test_data')
    
    simple_filename='%s/testset_1.fasta'%(test_data_dir)
    
    x = protfasta.read_fasta(simple_filename)
    assert len(x) == 9

    
    # check can read in a sequence correctly
    assert x[test_data['test1'][0]] == test_data['test1'][1]


def test_expect_unique_header_toggle():
    
    test_data_dir = protfasta._get_data('test_data')    
    simple_filename='%s/testset_1.fasta'%(test_data_dir)
    
    x = protfasta.read_fasta(simple_filename, expect_unique_header=False)
    assert len(x) == 9
    assert x[test_data['test1'][0]] == test_data['test1'][1]

    x = protfasta.read_fasta(simple_filename, expect_unique_header=True)
    assert len(x) == 9
    assert x[test_data['test1'][0]] == test_data['test1'][1]

    # bool only
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, expect_unique_header='dog')

    # bool only
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, expect_unique_header=1)


def test_header_parser():
    
    test_data_dir = protfasta._get_data('test_data')
    simple_filename='%s/testset_1.fasta'%(test_data_dir)

    def d(s):
        return s[0:10]

    def d_dumb(s):
        return "asas"

    def d_bad():
        return "asas"

    def d_bad2(s):
        return 1
    
    x = protfasta.read_fasta(simple_filename, header_parser=d)
    assert len(x) == 9
    assert x[test_data['test1'][0][0:10]] == test_data['test1'][1]

    # this dumb combination of settings means we overwrite the headers
    a = protfasta.read_fasta(simple_filename, header_parser=d_dumb, duplicate_sequence_action='ignore', expect_unique_header=False)
    assert len(a) == 1

    # now we at least avoid overwriting by setting the return type to be a list
    a = protfasta.read_fasta(simple_filename, header_parser=d_dumb, duplicate_sequence_action='ignore', expect_unique_header=False, return_list=True)
    assert len(a) == 9
                             
    # should fail because headers are duplicate
    with pytest.raises(ProtfastaException):
        assert  protfasta.read_fasta(simple_filename, header_parser=d_dumb)

    # bool only        
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, header_parser=d_bad)

    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, header_parser=d_bad2)



def test_duplicate_record_action():    

    test_data_dir = protfasta._get_data('test_data')
    duplicate_filename='%s/testset_duplicate.fasta'%(test_data_dir)
    simple_filename='%s/testset_1.fasta'%(test_data_dir)


    # this should be fine because simple_filename is valid
    assert len(protfasta.read_fasta(simple_filename, duplicate_record_action='fail')) == 9

    # this should fail because duplicate_filename has duplicates
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(duplicate_filename, duplicate_record_action='fail')

    # this should fail because this combination of options (i.e. implicit expect_unique=True)
    # will throw and error
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(duplicate_filename, duplicate_record_action='ignore')

    # THIS should fail because even though we've said remove, we are still expecting uniqe
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(duplicate_filename, duplicate_record_action='remove')

    x = protfasta.read_fasta(duplicate_filename, duplicate_record_action='remove', expect_unique_header=False)
    assert len(x) == 2


    x = protfasta.read_fasta(duplicate_filename, duplicate_record_action='ignore', expect_unique_header=False, return_list=True)
    assert len(x) == 3

    # this is not goood, BUT if we say expect uniuqe false, ignore duplicates and dont return a list we will use the first entry
    x = protfasta.read_fasta(duplicate_filename, duplicate_record_action='ignore', expect_unique_header=False)
    assert len(x) == 2


    # if we ignore or remove, same difference
    x = protfasta.read_fasta(duplicate_filename, duplicate_record_action='remove', expect_unique_header=False)
    assert len(x) == 2



def test_duplicate_sequence_action():    

    test_data_dir = protfasta._get_data('test_data')
    duplicate_filename='%s/testset_duplicate_seqs.fasta'%(test_data_dir)
    simple_filename='%s/testset_1.fasta'%(test_data_dir)

    # this should be fine because simple_filename is valid
    assert len(protfasta.read_fasta(simple_filename, duplicate_sequence_action='fail')) == 9

    # this should be fine because simple_filename is valid
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(duplicate_filename, duplicate_sequence_action='fail')

    # this should be fine because simple_filename is valid
    assert len(protfasta.read_fasta(duplicate_filename, duplicate_sequence_action='ignore')) == 3

    # remove duplciate sequence
    assert len(protfasta.read_fasta(duplicate_filename, duplicate_sequence_action='remove', verbose=True)) == 2

    # note only the sequences are duplicate, not the record
    assert len(protfasta.read_fasta(duplicate_filename, duplicate_record_action='remove', verbose=True)) == 3


def test_return_list():    
    test_data_dir = protfasta._get_data('test_data')
    duplicate_filename='%s/testset_duplicate_seqs.fasta'%(test_data_dir)
    duplicate_record='%s/testset_duplicate.fasta'%(test_data_dir)
    simple_filename='%s/testset_1.fasta'%(test_data_dir)

    x = protfasta.read_fasta(simple_filename, duplicate_sequence_action='fail')
    assert type(x) == dict


    x = protfasta.read_fasta(simple_filename, duplicate_sequence_action='fail', return_list=True)
    assert type(x) == list

    # show we can use return_list to read in a FASTA file with two identical records (note when we did this before and
    # return_list=False then len(x) == 2 because the dictionary overwrites
    x = protfasta.read_fasta(duplicate_record, duplicate_record_action='ignore', return_list=True, expect_unique_header=False)
    assert len(x) == 3

    x = protfasta.read_fasta(duplicate_record, duplicate_record_action='remove', return_list=True, expect_unique_header=False)
    assert len(x) == 2

