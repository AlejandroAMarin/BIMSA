import argparse

def read_file_by_index(file_path, index):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        index = (index) * 2
        if index < len(lines):
            #print("\n")
            print(lines[index] + lines[index+1],end="")
        else:
            print(f"Pair {int(index/2)} is out of range for file {file_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read a file by index')
    parser.add_argument('-f', '--file', type=str, help='Path to the file to be read')
    parser.add_argument('-i', '--index', type=int, help='Index of the line to be read')
    args = parser.parse_args()

    read_file_by_index(args.file, args.index)
