import protfasta

import pytest
import sys

TS1_FILE=['testset_1.fasta']
test_data_dir = shephard.get_data('test_data')



@pytest.fixture(scope='session', autouse=True)
def TS1(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    return TS1_proteome


@pytest.fixture(scope='session', autouse=True)
def TS1_domains(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    
    return TS1_proteome


@pytest.fixture(scope='session', autouse=True)
def TS1_domains2(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    
    return TS1_proteome

@pytest.fixture(scope='session', autouse=True)
def TS1_domains2_sites(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[3]))
    
    return TS1_proteome


@pytest.fixture(scope='session', autouse=True)
def TS1_domains2_sites_tracks(request):    
    TS1_proteome = build_proteome(TS1_FILE[0])
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[1]))
    si_domains.add_domains_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[2]))
    si_sites.add_sites_from_file(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[3]))
    si_tracks.read_in_tracks(TS1_proteome, '%s/%s' %(test_data_dir, TS1_FILE[4]))
    
    return TS1_proteome


#si_sites.add_sites(proteome, 'proteomescout_ptm_sites.csv')
"""
@pytest.fixture(scope='session', autouse=True)
def GS6_CP(request):
    GS6_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, GS6_FILES[1]), "%s/%s"%(test_data_dir, GS6_FILES[0])) 
    GS6_CP = GS6_CO.proteinTrajectoryList[0]    
    return GS6_CP


@pytest.fixture(scope='session', autouse=True)
def NTL9_CO(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    return NTL9_CO

@pytest.fixture(scope='session', autouse=True)
def NTL9_CP(request):    
    NTL9_CO = cttrajectory.CTTrajectory("%s/%s"%(test_data_dir, NTL9_FILES[1]), "%s/%s"%(test_data_dir, NTL9_FILES[0])) 
    NTL9_CP = NTL9_CO.proteinTrajectoryList[0]    
    return NTL9_CP
"""
