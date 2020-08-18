# Python script for merging path names lists

import os

refsoc = '/u/home/m/mmckinne/hab_proj/Isca/src/extra/model/socrates/path_names'
refreg = '/u/home/m/mmckinne/hab_proj/Isca/src/extra/model/isca/path_names'

regconts = open(refreg,'r')
reglines = regconts.readlines()
regconts.close()

socconts = open(refsoc,'r')
soclines = socconts.readlines() 
socconts.close()

newlines = []

for line in reglines:
 if line in soclines:
  continue
 else:
  newlines.append(line)

socedit = open(refsoc,'w+')
for line in newlines:
 socedit.write(line)
for line in soclines:
 socedit.write(line)
socedit.close()

