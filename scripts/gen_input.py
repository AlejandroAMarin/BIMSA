import os
import sys
import argparse

bindir = "../tools/generate_dataset/"  # Include path to generate_dataset binary
outputdir = "../inputs/"  # Dir to store the files

def remove_first_line(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()  # Read all lines
        
        # Write the lines back, excluding the first line
        with open(filename, 'w') as file:
            file.writelines(lines[1:])
        print(f"First line removed from {filename}")
    except Exception as e:
        print(f"Error removing first line from {filename}: {e}")

def run(nr_pairs, seq_length, error):

    # The single values are passed as arguments, no need to iterate over lists
    run_cmd = bindir + "generate_dataset --o #output --n #nr_patterns --l #seq_length --e #error"
    os.system("mkdir -p " + outputdir)
    
    out_file_name = outputdir + "n" + str(nr_pairs) + "_l" + str(seq_length) + "_e" + str(int(error*100)) + ".seq"
    r_cmd = run_cmd.replace("#output", out_file_name)
    r_cmd = r_cmd.replace("#nr_patterns", str(nr_pairs))
    r_cmd = r_cmd.replace("#seq_length", str(seq_length))
    r_cmd = r_cmd.replace("#error", str(error))

    print(f"Generating [NR_PAIRS={nr_pairs}, SEQ_LENGTH={seq_length}, ERROR={int(error*100)}%] at {out_file_name}")
    try:
        os.system(r_cmd)
        remove_first_line(out_file_name)  # Remove the first line after generation
    except:
        pass

    print("Done!")


def main():
    parser = argparse.ArgumentParser(description="Generate datasets with specified parameters")
    parser.add_argument('-n', '--nr_pairs', type=int, required=True, help='Number of pairs (NR_PAIRS)')
    parser.add_argument('-l', '--seq_length', type=int, required=True, help='Sequence length (SEQ_LENGTH)')
    parser.add_argument('-e', '--error', type=float, required=True, help='Error rate (ERROR)')

    args = parser.parse_args()

    run(args.nr_pairs, args.seq_length, args.error)


if __name__ == "__main__":
    main()

