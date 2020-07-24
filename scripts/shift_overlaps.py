#!/usr/bin/python
from sys import argv
from fileinput import FileInput

if len(argv) < 3:
  print("2 arguments are required: (1) a file of alignments (2) a file containing the corresponding diBELLA overlaps")

prev_set=set()
filename=argv[1]
with open(filename) as file:
  for line in file:
    linesplit=line.split()
    prev_set.add( (linesplit[0],linesplit[1]) )

filename=argv[2]
keep=[]
with FileInput(filename, inplace=False, backup='.old') as file:
  for line in file:
    linesplit=line.split()
    if (linesplit[0],linesplit[1]) not in prev_set:
      keep.append(line)

outputfile=filename+".new"
with open(outputfile, 'w+') as ofile:
  for item in keep:
    ofile.write("%s" % item)

