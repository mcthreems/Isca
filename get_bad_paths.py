# Python script for finding a bad path name

import os

fname = './src/extra/model/socrates/path_names'

fid = open(fname,'r')
fcont = fid.readlines()
fid.close()

for pname in fcont:
 try:
  os.system('ls ./src/'+pname)
 except:
  print(pname)

