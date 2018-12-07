#!usr/bin/env python

from subprocess import call
import sys
import nibabel as nib
import numpy as np
import pandas as pd

csv = sys.argv[1]

df = pd.read_csv(csv)

for i in df.index:
	T1 = df.t1[i]
	T2 = df.t2[i]

	#extract brain from T1 and T2
	t1_brain = T1[:-7]+'_bet.nii.gz'
	call(['bet', T1, t1_brain])

	t2_brain = T2[:-7]+'_bet.nii.gz'
	call(['bet', T2, t2_brain])

	#rigid registration of brain from T2 to T1
	reg_brain_file =  T1[:-7]+'_brain_rigid_T1_reg.nii.gz'
	rigid_transform = T1[:-7]+'_brain_rigid_T1_reg.mat'
	call(['flirt', '-cost', 'normmi', 'dof', '6', '-in', t2_brain, '-ref', t1_brain, '-out', reg_brain_file, '-omat', 			rigid_transform]) 

	#apply transform to T1 masks to register them to T2 
