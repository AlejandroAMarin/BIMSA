import os 
import sys

bindir    = "../tools/generate_dataset/" # Include path to generate_dataset binary
outputdir = "../inputs/"                 # Dir to store the files

def run():
    
    # The files will be created by the cartesian product of these parameters:
    NR_PAIRS    = [256, 512]         # Defines the numbers of sequences
    SEQ_LENGTHS = [100, 1000]        # Defines the lengths for the sequences
    ERRORS      = [0.02, 0.05, 0.1]  # Defines the errors (e.g., 0.1 =  10% error)

    num_files           = len(NR_PAIRS) * len(SEQ_LENGTHS) * len(ERRORS)
    num_processed_files = 0

    run_cmd = bindir + "generate_dataset --o #output --n #nr_patterns --l #seq_length --e #error"
    os.system("mkdir -p " + outputdir)    
    for p in NR_PAIRS:
        for s in SEQ_LENGTHS:
            for e in ERRORS:
                out_file_name = outputdir + "n" + str(p) + "_l" + str(s) + "_e" + str(int(e*100)) + ".seq"          
                r_cmd = run_cmd.replace("#output",    out_file_name)
                r_cmd = r_cmd.replace("#nr_patterns", str(p))
                r_cmd = r_cmd.replace("#seq_length",  str(s))
                r_cmd = r_cmd.replace("#error",       str(e))
                        
                print ("["+  "{0:0=3d}".format(int((num_processed_files / num_files) * 100))+ "%] Generating [NR_PAIRS=" + str(p) + ", SEQ_LENGTH=" + str(s) + ", ERROR=" + str(int(e*100)) + "%] at " +  out_file_name)
                num_processed_files = num_processed_files + 1
                try:
                    os.system(r_cmd) 
                except:  
                    pass 

    print ("["+  "{0:0=3d}".format(int((num_processed_files / num_files) * 100))+ "%] Done!")
def main():
    run()

if __name__ == "__main__":
    main()
