from subprocess import check_call
import os
from glob import glob

imgs = glob('/data/henry7/PBR/subjects/mse*/nii/tca_roi/*_i1_filled_cerebellum_ci_1.nii.gz')
for img in imgs:
    if not os.path.exists(img.split('_i1_filled_cerebellum_ci_1.nii.gz')[0] + '.roi'):
        cmd = ['Jim', img]
        print(cmd)
        check_call(cmd)
