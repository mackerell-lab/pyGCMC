"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import sys

class GCMCParameters:
    
    # Get and validate fragment excess chemical potential
    def get_fragmuex(self, fragmuex):
        fragmuex = [float(i) for i in fragmuex if len(i) > 0]
        if len(fragmuex) != len(self.fragmentName):
            print("Error: fragmuex number not match")
            sys.exit(1)
        self.fragmuex = fragmuex
    
    # Get and validate fragment concentration
    def get_fragconc(self, fragconc):
        fragconc = [float(i) for i in fragconc if len(i) > 0]
        if len(fragconc) != len(self.fragmentName):
            print("Error: fragconc number not match")
            sys.exit(1)
        self.fragconc = fragconc

    # Get and validate fragment configurational bias
    def get_fragconf(self, fragconf):
        fragconf = [int(i) for i in fragconf if len(i) > 0]

        for i in fragconf:
            if i <= 0:
                print("Error: fragconf number <= 0")
                sys.exit(1)

        if len(fragconf) == 1:
            self.fragconf = fragconf * len(self.fragmentName)
        elif len(fragconf) == len(self.fragmentName):
            self.fragconf = fragconf
        else:
            print("Error: fragconf number not match")
            sys.exit(1)
        
        self.configurational_bias = not all(x == 1 for x in self.fragconf)

    # Get and validate Monte Carlo time for each fragment
    def get_mctime(self, mctime):
        mctime = [float(i) for i in mctime if len(i) > 0]
        if len(mctime) != len(self.fragmentName):
            print("Error: mctime number not match")
            sys.exit(1)
        self.mctime = mctime

    # Display all parameters
    def show_parameters(self):
        print(f"MC steps: {self.mcsteps}")
        print("Solute Name: \t\t", '\t\t'.join(self.fragmentName))
        print("Solute Muex: \t\t", '\t\t'.join([str(i) for i in self.fragmuex]))
        print("Solute Conc: \t\t", '\t\t'.join([str(i) for i in self.fragconc]))
        print("Solute ConfB: \t", '\t\t'.join([str(i) for i in self.fragconf]))
        print("Solute mcTime: \t", '\t\t'.join([str(i) for i in self.mctime]))


