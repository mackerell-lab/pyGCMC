"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Abhishek A. Kognole, 
    	Aoxiang Tao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""

from .gcmc import GCMC
from .packages import *


def main():
    # file_output = open('Analyze_output.txt', 'w')
    # original_output = sys.stdout
    # sys.stdout = Tee(sys.stdout, file_output)


    starttime = time.time()
    
    


    parser = argparse.ArgumentParser(description='Perform GCMC Simulation')

    parser.add_argument('-p', '--paramfile', type=str, required=True, 
                        help='[Required] input parameter file')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='[Optional] verbose output')
    parser.add_argument('--logfile', type=str, 
                        help='[Optional] log file, if not specified, then output will be stdout')
    parser.add_argument('--debug', action='store_true', 
                        help='[Optional] for debug purpose')
    parser.add_argument('--version', action='version', version='GCMC version', 
                        help='Show program\'s version number and exit')
    args = parser.parse_args()



    out_file = args.logfile

    if out_file is not None:
        file_output = open(out_file, 'w')
        original_output = sys.stdout
        sys.stdout = Tee(sys.stdout, file_output)
        # print(f"Using output file: {out_file}")


    print('Start GCMC simulation at %s...' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))



    if out_file is not None:
        print(f"Using output file: {out_file}")

    gcmc = GCMC()

    try:
        gcmc.load_parameters(args.paramfile)
    except:
        print(f"Error reading parameter file: {args.paramfile}")
        sys.exit(1)
        


    endtime = time.time()
    print(f"Python time used: {endtime - starttime} s")


    
    print('GCMC simulation finished at %s...' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    endtime = time.time()
    print(f"Time used: {endtime - starttime} s")


    # if top_file is not None:
    #     try:
    #         # top = protein_data.read_top(top_file)
    #         protein_data.read_top(top_file)
    #     except:
    #         print(f"Error reading top file: {top_file}")
    #         sys.exit(1)
    #     # print(f"top atom number: {len(top)}")

    if out_file is not None:
        sys.stdout = original_output
        file_output.close()