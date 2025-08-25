"""
protfasta - A simple but robuts FASTA parser explicitly for protein sequences.

This includes internal configuration stuff


.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other. 

"""

STANDARD_CONVERSION = {'B':'N',
                       'O':'K',
                       'U':'C',
                       'X':'G',
                       'Z':'Q',
                       '*':'',
                       '-':'',
                       ' ':''}

STANDARD_CONVERSION_WITH_GAP = {'B':'N',
                                'O':'K',
                                'U':'C',
                                'X':'G',
                                'Z':'Q',
                                ' ':'',
                                '*':''}

STANDARD_AAS = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
STANDARD_AAS_WITH_GAP = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
