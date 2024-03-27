# Python script for adding line to many files

import glob

fcard = 'fandry_qmax_higheq' #wildcard to get files to add line to
refline = "diag.add_field('atmosphere','dt_qg_diffusion', time_avg=True)" #line after which to place new line
newline = ("diag.add_field('atmosphere', 'fnq_term', time_avg=True)\n"
           "diag.add_field('atmosphere', 'enq_term', time_avg=True)")

fnames = glob.glob(fcard+'*.py')
for fname in fnames:
 fid = open(fname,'r')
 ftx = fid.read()
 fid.close()

 if refline in ftx and newline not in ftx:
  getloc = ftx.index(refline) + len(refline)
  newtx = ftx[:getloc] + '\n' + newline + ftx[getloc:]
  fid = open(fname,'w+')
  fid.write(newtx)
  fid.close()
  print('New line added to '+fname)
 elif newline in ftx:
  print('New line already present for '+fname)
 else:
  print('Reference line not found in '+fname)

