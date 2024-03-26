import argparse
import os
import subprocess

def run_make_command(code_directory, args):
    exec_call = ['make', f'NR_TASKLETS={args.tasklets}', f'NR_DPUS={args.dpus}', 
        f'BL={args.block_size}', f'BLI={args.block_size_inputs}']
    if args.perf_instructions and not args.perf_cycles:
        exec_call.append("PERF=1")
        exec_call.append("PERF_INSTRUCTIONS=1")
    if args.perf_cycles and not args.perf_instructions:
        exec_call.append("PERF=1")
        exec_call.append("PERF_CYLES=1")
    if args.perf_instructions and args.perf_cycles:
        print("Can not profile cycles and instructions at the same time. Using INSTRUCTIONS as dafault...")
        exec_call.append("PERF=1")
        exec_call.append("PERF_INSTRUCTIONS=1")
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
    exec_call = ['./bin/ulsapim_host', '-i', f"{args.file}", '-s', f'{args.sets}']
    try:
        subprocess.check_call(exec_call,cwd=code_directory)
    except subprocess.CalledProcessError as e:
        print(f"C program execution failed with error code {e.returncode}")
        exit(1)

def main():
    parser = argparse.ArgumentParser(description='Python program to run make and execute a C program with an input file.')
    parser.add_argument('-f', '--file', help='Path to the input file', required=True)
    parser.add_argument('-t', '--tasklets', help='Number of tasklets (default 1)', required=False, default=1)
    parser.add_argument('-d', '--dpus', help='Number of DPUs (dafault 1)', required=False, default=1)
    parser.add_argument('-bl', '--block_size', help='Wavefront transfer and cache size in power of 2', required=False, default=8)
    parser.add_argument('-bli', '--block_size_inputs', help='Sequence and CIGAR transfer and cache size in power of 2', required=False, default=5)
    parser.add_argument('-bd', '--banded', help='Use banded heuristic', required=False, action='store_true')
    parser.add_argument('-ad', '--adaptive', help='Use adaptive heuristic', required=False, action='store_true')
    parser.add_argument('-s', '--sets', help='Define different sequence size sets', required=True)
    parser.add_argument('-pi', '--perf_instructions', help='Profile instructions', action='store_true')
    parser.add_argument('-pc', '--perf_cycles', help='Profile cycles', action='store_true')

    args = parser.parse_args()
    code_directory = "upmem/"
    run_make_command(code_directory, args)
    run_c_program(code_directory, args)

if __name__ == '__main__':
    main()
