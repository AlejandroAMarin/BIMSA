import argparse
import os
import subprocess

def run_make_command(code_directory, args):
    exec_call = ['make', f'NR_TASKLETS={args.tasklets}', f'NR_DPUS={args.dpus}', 
        f'WFT={args.wf_trans}', f'SEQT={args.seq_trans}', 
        f'CIGART={args.cigar_trans}', 
        f'MAX_DISTANCE_THRESHOLD={args.max_distance}', f'MAX_ERROR={round((args.max_distance * 2)/int(args.size), 2)}',
        f'DYNAMIC={args.dynamic}', f'BATCH_SIZE={args.batch_size}',
        f'PRINT={args.print}', f'DEBUG={args.debug}']
    if args.banded:
        exec_call.append("BANDED=1")
    if args.adaptive:
        exec_call.append("ADAPTIVE=1")

    print(exec_call)
    try:
        subprocess.check_call(['make', 'clean'], cwd=code_directory)
        subprocess.check_call(exec_call, cwd=code_directory)
    except subprocess.CalledProcessError as e:
        print(f"make command failed with error code {e.returncode}")
        exit(1)

def run_c_program(code_directory, args):
    exec_call = ['./bin/bimsa_host', '-i', f"../{args.file}", '-s', f'{args.size}']
    try:
        subprocess.check_call(exec_call,cwd=code_directory)
    except subprocess.CalledProcessError as e:
        print(f"C program execution failed with error code {e.returncode}")
        exit(1)

def main():
    parser = argparse.ArgumentParser(description='Python program to run make and execute a C program with an input file.')
    parser.add_argument('-f', '--file', help='Path to the input file', required=True)
    parser.add_argument('-t', '--tasklets', help='Number of tasklets (default 12)', required=False, default=12)
    parser.add_argument('-d', '--dpus', help='Number of DPUs (dafault 2500)', required=False, default=2500)
    parser.add_argument('-wt', '--wf_trans', help='Wavefront transfer and cache size in power of 2', required=False, default=9)
    parser.add_argument('-st', '--seq_trans', help='Sequence and CIGAR transfer and cache size in power of 2', required=False, default=3)
    parser.add_argument('-ct', '--cigar_trans', help='Sequence and CIGAR transfer and cache size in power of 2', required=False, default=9)
    parser.add_argument('-s', '--size', help='Define maximum sequence size', required=True)
    parser.add_argument('-m', '--max_distance', help='Maximum wavefront distance (default 5000)', required=False, default=5000)
    parser.add_argument('-b', '--batch_size', help='Number of pairs to be executed in batches (default 0)', required=False, default=0)
    parser.add_argument('-dn', '--dynamic', help='Asign pairs dynamically to threads', required=False, action='store_true')
    parser.add_argument('-p', '--print', help='Print extra dpu execution information', required=False, action='store_true')
    parser.add_argument('-db', '--debug', help='Print extra validation information', required=False, action='store_true')
    parser.add_argument('-bd', '--banded', help='Use banded heuristic', required=False, action='store_true')
    parser.add_argument('-ad', '--adaptive', help='Use adaptive heuristic', required=False, action='store_true')

    args = parser.parse_args()
    code_directory = "upmem/"
    run_make_command(code_directory, args)
    run_c_program(code_directory, args)

if __name__ == '__main__':
    main()
