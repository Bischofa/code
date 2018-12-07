from glob import glob
from subprocess import check_call
import os
import argparse
from pbr.base import _get_output

if __name__ == '__main__':

 # code to grab command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('MSEID', help='txt file of the mseids you want to get')
    args = parser.parse_args()
    
 # open the file passed as args.MSEID as f, strips off empty lines in the txt file of mseids
    with open(args.MSEID, 'r') as f:
        mseid_list = [x.strip() for x in f.readlines()]

 # iterate through each mseid in the file
    for mseid in mseid_list:
        img = glob('{0}/{1}/nii/*_IRSPGR_new.nii.gz'.format(_get_output(mseid), mseid))
        
        if img:
            cmd = ['Jim', img[0]]
            print(cmd)
            check_call(cmd)
        else:
            print("missing {}" .format(mseid))


    


