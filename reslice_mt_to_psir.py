#!/usr/bin/env python
from subprocess import call
import os
import sys
import numpy as np
import nibabel as nib
import glob
import argparse

def psir_slices(folder):
    # process MT
    mt = os.path.join(folder, "MT_Map.nii.gz")
    mt_thr = mt[:-7]+'_thr.nii.gz'
    mt_psir_slices = mt_thr[:-7]+'_psir.nii.gz'
    call(['fslmaths', mt, '-thr', '0', '-uthr', '1', mt_thr])
    call(['fslroi', mt_thr, mt_psir_slices, '0', '-1', '0', '-1', '31', '2'])
    average_slices(mt_psir_slices)

    # process MTR
    mtr = os.path.join(folder, "MTR.nii.gz")
    mtr_psir_slices = mtr[:-7]+'_psir.nii.gz' 
    call(['fslroi', mtr, mtr_psir_slices, '0', '-1', '0', '-1', '31', '2'])
    average_slices(mtr_psir_slices)

    # process R
    r = os.path.join(folder, "R.nii.gz")
    r_psir_slices = r[:-7]+'_psir.nii.gz'
    r_psir_slices_thr = r_psir_slices[:-7]+'_thr.nii.gz'
    r_psir_slices_avg = r_psir_slices[:-7]+'_thr_avg.nii.gz'
    call(['fslroi', r, r_psir_slices, '0', '-1', '0', '-1', '31', '2'])

    rpsir = nib.load(r_psir_slices)
    rpsir_affine = rpsir.get_affine()
    rpsir_data = np.array(rpsir.get_data())
    inverse_data = np.ones(np.shape(rpsir_data)) / rpsir_data
    inverse_data = np.nan_to_num(inverse_data)
    data_thr0 = inverse_data > 0
    data_thr5000 = inverse_data < 5000
    thr_data = inverse_data * data_thr0 * data_thr5000
    avg_data = np.average(thr_data, 2)

    nib.save(nib.Nifti1Image(thr_data, rpsir_affine), r_psir_slices_thr)
    nib.save(nib.Nifti1Image(avg_data, rpsir_affine), r_psir_slices_avg)
    
    return

def average_slices(roi_nii):
    img = nib.load(roi_nii)
    affine = img.get_affine()
    data = img.get_data()
    avg_data = np.average(data, 2)
    nib.save(nib.Nifti1Image(avg_data,affine), roi_nii[:-7]+"_avg.nii.gz")
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mt_folders", nargs="+")
    args = parser.parse_args()

    for folder in args.mt_folders:
        print folder
        psir_slices(folder)
