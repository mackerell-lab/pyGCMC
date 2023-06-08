"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Abhishek A. Kognole, 
    	Aoxiang Tao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""


import numpy as np
import pkg_resources
import os


class GCMC:
    
    def __init__(self):

        # find the path of the mols files        
        self.mols_str_path = pkg_resources.resource_filename(__name__, 'mols') + '/'

        
        self.get_simulation()


        
        
    def get_fragment(self,fragment = None):
        if fragment is None:
            self.fragmentName = ['acey', 'benx', 'dmee', 'form', 'imia', 'mamy', 'meoh', 'prpx', 'sol']


            # fragmentPdb = [pmd.load_file('%s%s.pdb' % (self.mols_str_path,i)) for i in fragmentName]
            # fragmentPsf = [pmd.load_file('%s%s.psf' % (self.mols_str_path,i)) for i in fragmentName]
            
            # for i in range(len(fragmentName)):
            #     fragmentPsf[i].positions = fragmentPdb[i].positions
            # self.fragment = fragmentPsf

            # self.fragment = [(fragmentPdb[i],fragmentPsf[i]) for i in range(len(fragmentName))]


    def get_simulation(self):

        # NBFIX  rmin=<charmm_rmin>/2^(1/6), eps=4.184*<charmm_eps>
        # NB  rmin=<charmm_rmin>/2^(1/6)*2, eps=-4.184*<charmm_eps>

        # self.nbfix_dict = {i:[self.params.nbfix_types[i][1] / 2**(1./6) * 0.1,self.params.nbfix_types[i][0] * 4.184 ] 
        #             for i in self.params.nbfix_types}
        
        # self.nb_dict = {i:[self.params.atom_types_str[i].rmin / 2 ** (1./6) * 0.1 * 2,self.params.atom_types_str[i].epsilon *-4.184] 
        #                 for i in self.params.atom_types_str}
    
        
        
        self.get_fragment()

        

