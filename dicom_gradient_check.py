from glob import glob
from os.path import join as opj
import argparse
from pydicom import dcmread 
import pandas as pd
from nipype.utils.filemanip import load_json

if __name__ == '__main__':

    # code to grab command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('MSEID', help='txt file of the mseids you want to get')
    parser.add_argument('OUT', help='output name of the output csv')
    args = parser.parse_args()
    
    # load the heuristic
    heuristic = load_json('/data/henry2/arajesh/PBR/pbr/heuristic.json')['filetype_mapper']

    # make a pandas dataframe where we will store output
    df = pd.DataFrame()
 
    # open the file passed as args.MSEID as f
    with open(args.MSEID, 'r') as f:
        mseid_list = [x.strip() for x in f.readlines()]
  
    # iterate through each mseid in the file, and also get a numerical number (idx) for how many times we go thru the loop
    for idx,mseid in enumerate(mseid_list):
    
        # grabs a list of all files that fufill criteria of the glob function
        dicom_series_glob = glob(opj('/working/henry_temp/PBR/dicoms/', '*{}*/E*/*'.format(mseid))) 
 
        # go through each series for tha tgiven mseid
        for series in dicom_series_glob:
 
            print('looking at series {}'.format(series.split('/')[-1]))
   
            # grab first dicm in each series folder
            first_dcm_glob = glob(opj(series, '*.DCM'))
            if not first_dcm_glob:
                continue
            else:
                first_dcm = first_dcm_glob[0]
            first_dcm_data = dcmread(first_dcm)
 
            # grab the modality of the dicom
            try:
                modality = first_dcm_data[0x08, 0x103e].value
                modality = modality.replace(' ', '_')
                modality = modality.replace('-', '_')
            except:
                modality = ''            
          
            # try to find the filetype in the heuristic
            try:
               # if the heuristic cotains the modality this will run
               filetype = heuristic[modality] 
            except:
               # if modlaity is not in the heuristic then the try loop will erroir causing an exception and making filetype be equal to 'NA' 
               filetype='NA'
   
            print('modality', modality)
            # this means we have the correct dicom series
            if filetype == 'T1':
              
               print('found t1, with modality {}'.format(modality))
 
               # do same logic with try/except loop as above but with checking the  
               try:
                   gradient = first_dcm_data[0x43,0x102D].value.strip()
               except:
                   gradient = '' 
            
               # set the row of idx and column of 'mseid' to be the mseid value
               df.loc[idx, 'mseid'] = mseid

               print('gradient', gradient)
               try:
                   if 'w' in gradient:
                       df.loc[idx, 'gradient'] = '3D'
                   else:
                       df.loc[idx, 'gradient'] = '2D'
                   df.loc[idx, 'gradient_val'] = gradient
               except:
                   print('binary file')               
                   try:
                       if b'w' in gradient:
                           df.loc[idx, 'gradient'] = '3D'
                       else:
                           df.loc[idx, 'gradient'] = '2D'
                       df.loc[idx, 'gradient_val'] = gradient.decode('ascii')
                   except:
                       print('idk broooo')
    
               break
     
    print('saving data to {}'.format(args.OUT)) 
    df.to_csv(args.OUT, index=False)
