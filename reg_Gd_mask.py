#!/usr/bin/env python

from subprocess import call
import sys
import nibabel as nib
import numpy as np
import pandas as pd

csv = sys.argv[1]

df = pd.read_csv(csv)

for i in df.index:
	T1_pre = df.pre[i]
	T1_post = df.post[i]
	DTI_FA = df.FA[i]
	mask = df.Gd_mask[i]

	#extract brain from T1 pre and post, rather use low fractional intensity threshold to make sure all brain voxels are included, as images are subtracted later
	pre_brain = T1_pre[:-7]+'_bet.nii.gz'
	call(['bet', T1_pre, pre_brain, '-B', '-f', '0.20'])

	post_brain = T1_post[:-7]+'_bet.nii.gz'
	call(['bet', T1_post, post_brain, '-B', '-f', '0.20'])

	#register postGdT1 to preGdT1
	reg_brain = post_brain[:-7]+'reg.nii.gz'
	transform = post_brain[:-7]+'transform.mat
	call(['flirt','-in', post_brain, '-ref', pre_brain, '-out', reg_brain, '-cost', 'normmi', '-dof', '6', '-omat', transform]) 
	#'-interp', 'nearestneighbour', 

	#register DTI to preGdT1
	reg_FA = DTI_FA[:-7]+'reg.nii.gz'
	transform_FA = DTI_FA[:-7]='transform.mat'
	call(['flirt', '-in', DTI_FA, '-ref', pre_brain, '-out', reg_FA, '-cost', 'normmi', '-dof', '6', '-omat', transform_FA])

	#registering Gd_mask to original DTI by registering inverse DTI matrix on Gd_mask
	mask_reg = mask[:-7]+'_reg.nii.gz'
	call(['flirt', '-in', 'mask', '-ref', 'DTI_FA', '-out', mask_reg, '-interp', 'nearestneighbour', '-applyxfm', '-init', transform, '-concat', '-convert_xfm', '-inverse', transform_FA])

	'''
	#subtract brain extracted preT1 from postT1
	sub_file = pre_brain[:-7]+_'sub.nii.gz'
	call(['fslmaths', reg_brain, '-sub', pre_brain, sub_file])'''

	
