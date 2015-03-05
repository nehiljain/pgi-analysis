__author__ = 'nehil'



import sys
import argparse


parser = argparse.ArgumentParser(description="Takes the shell command and a directory name and creates a shell script with all scripts")
parser.add_argument("-cmd","--command", help="Tested Shell Command example string.")
parser.add_argument("-extn","--extension", help="The unique extension of the file. Can be .txt, .sorted.bam etc")
args = parser.parse_args()



if args.command:
    print(args.command)
else:

if args.extension:
    print(args.extension)


