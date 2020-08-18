# Python script for submitting an experiment file

import os

def func_sub(expname):
 i_exp_temp = open('exp_sub','r')
 exp_temp = i_exp_temp.read()
 i_exp_temp.close

 edit_loc1 = exp_temp.index('set exp_file  = ')+len('set exp_file  = ')
 edit_loc2 = edit_loc1+exp_temp[edit_loc1:].index('\n')

 exp_new = exp_temp[:edit_loc1]+expname+'.py'+exp_temp[edit_loc2:]

 i_exp_new = open('sub_'+expname,'w+')
 i_exp_new.write(exp_new)
 i_exp_new.close()

