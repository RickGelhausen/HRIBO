#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='Minimal value of a input filelist')
    parser.add_argument('files', type=argparse.FileType('r'), nargs='+')
    args = parser.parse_args()
    minimal_value=int(100000000)
    for file in args.files:
        for line in file:
            processed_line=int(line.rstrip("\r\n"))
            if processed_line < minimal_value:
                minimal_value=processed_line
    print(str(minimal_value))
if __name__ == '__main__':
    main()
