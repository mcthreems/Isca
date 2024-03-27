# Python script for getting restart point for an experiment file

import os, glob

def get_start(expname):
 fnames = glob.glob('/u/scratch/m/mmckinne/isca_data/'+expname+'/restarts/res*')
 if len(fnames) == 0:
  return 1
 else:
  rnums = [int(fname[-11:-7]) for fname in fnames]
  rmax = max(rnums)
  srun = rmax+1
  return srun
