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

    
    test_data_dir = protfasta.get_data('test_data')
    
    simple_filename='%s/testset_1.fasta'%(test_data_dir)
    
    x = protfasta.read_fasta(simple_filename)
    assert len(x) == 9

    
    # check can read in a sequence correctly
    assert x[test_data['test1'][0]] == test_data['test1'][1]


def test_expect_unique_toggle():
    
    test_data_dir = protfasta.get_data('test_data')    
    simple_filename='%s/testset_1.fasta'%(test_data_dir)
    
    x = protfasta.read_fasta(simple_filename, expect_unique=False)
    assert len(x) == 9
    assert x[test_data['test1'][0]] == test_data['test1'][1]

    x = protfasta.read_fasta(simple_filename, expect_unique=True)
    assert len(x) == 9
    assert x[test_data['test1'][0]] == test_data['test1'][1]

    # bool only
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, expect_unique='dog')

    # bool only
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, expect_unique=1)


def test_header_parser():

    
    test_data_dir = protfasta.get_data('test_data')
    
    simple_filename='%s/testset_1.fasta'%(test_data_dir)

    def d(s):
        return s[0:10]

    def d_bad():
        return "asas"

    def d_bad2(s):
        return 1
    
    x = protfasta.read_fasta(simple_filename, header_parser=d)
    assert len(x) == 9
    assert x[test_data['test1'][0][0:10]] == test_data['test1'][1]

    # bool only        
    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, header_parser=d_bad)

    with pytest.raises(ProtfastaException):
        assert protfasta.read_fasta(simple_filename, header_parser=d_bad2)



    


