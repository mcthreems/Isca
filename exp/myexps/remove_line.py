# Python script for removing lines from text files using a search term

import glob


pyfiles = glob.glob('*.py')

for pyfile in pyfiles:
 fid = open(pyfile,'r')
 ftext = fid.read()
 fid.close()
 
 if sterm not in ftext:
  continue
 
 tind = ftext.index(sterm)
 tind1 = tind
 tind2 = tind
 
 while ftext[tind1] != '\n':
  tind1 += -1
 while ftext[tind2] != '\n':
  tind2 += 1
  
 newt1 = ftext[:tind1]
 newt2 = ftext[tind2:]
 newt = newt1 + newt2
 
 fid = open(pyfile,'w')
 fid.write(newt)
 fid.close()

