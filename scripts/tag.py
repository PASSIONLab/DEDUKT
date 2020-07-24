#!/usr/bin/python
from sys import argv
from fileinput import FileInput

if len(argv) < 3:
  print("2 arguments are required: (1) a text file containing the list of read tags, 1 per line, and (2) a file containing the corresponding diBELLA untagged alignment output")
  quit()

tag_map=dict()
filename=argv[1]
with open(filename) as file:
  tag_map={k:v for (k,v) in enumerate(file)}

#now open alignments file and replace first two words with of each line with value from dictionary
space=' '
filename=argv[2]
with FileInput(filename, inplace=True, backup='.bak') as file:
  for line in file:
    linesplit=line.split()
    linesplit[0]=tag_map[(int(linesplit[0])-1)].rstrip()
    linesplit[1]=tag_map[(int(linesplit[1])-1)].rstrip()
    lineout=' '.join(x for x in linesplit)
    print(lineout)
