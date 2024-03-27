# Python script for slimming data down to last 5 years

import os, glob

datdir = '/u/scratch/m/mmckinne/isca_data'
keeplast = 10

exnames = glob.glob(datdir+'/*')

for exname in exnames:
 runs = glob.glob(exname+'/run*')
 if len(runs) == 0:
  continue
 else:
  rnums = [int(irun[-4:]) for irun in runs]
  if len(rnums) > 1:
   rnums.sort()
   gaps = [1]+[rnums[ii] - rnums[ii-1] for ii in range(1,len(rnums))]
   discont = [igap > 1 for igap in gaps]
   if True in discont:
    dloc = discont.index(True) - 1
    maxnum = rnums[dloc]
   else:
    maxnum = max(rnums)
  else:
   maxnum = max(rnums)
  
  del1 = maxnum - keeplast
  delnums = [del1 - ii for ii in range(0,del1)]
  for delnum in delnums:
   if os.path.isdir(exname+'/run'+str(delnum).zfill(4)):
    os.system('rm -r '+exname+'/run'+str(delnum).zfill(4))

