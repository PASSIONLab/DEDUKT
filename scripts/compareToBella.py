#!/usr/bin/python
#
# @author:mme
#
from sys import argv

if len(argv) < 2:
  print("2 arguments required: the path to a bella output file and another output file in bella format")
  quit()

bellas={}
tot_bella=0
with open(argv[1], 'r') as bella_file:
  for line in bella_file:
    sortedline=sorted( line.split()[:3] )
    bellas[ (sortedline[1], sortedline[2]) ] = sortedline[0]
    tot_bella+=1

hits=0
misses=0
mismatches=0
with open(argv[2], 'r') as file:
  for line in file:
    sortedline=sorted( line.split()[:3] )
    key=(sortedline[1],sortedline[2])
    if key in bellas:
      hits+=1
      if sortedline[0] != bellas[key]: 
        mismatches+=1
    else: misses+=1

tot_dibella=hits+misses

print("Total BELLA: ", tot_bella)
print("Total diBELLA: ", tot_dibella)
print("Total hits: ", hits)
print("Total misses: ", misses)
print("Total mismatches (w.r.t. number of shared k-mers): ", mismatches)