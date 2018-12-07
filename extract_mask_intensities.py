''' This script allows you to extract the mean intensities of an image that is covered by a mask '''
import argparse
import numpy as np 
import nibabel as nib 

if __name__ == '__main__':

    parser = argparse.ArgumentParser('This script allows you to extract the mean intensities of an image that is covered by a mask ')
    parser.add_argument('img', help='Image that intensities will be extracted from')
    parser.add_argument('mask', help='Binary mask of an ROI where you want to extract the intensities')
    args = parser.parse_args()

    img = nib.load(args.img).get_data()
    mask = nib.load(args.mask).get_data()

    mean = np.mean(img[np.where(mask == 1)]) 
    std = np.std(img[np.where(mask == 1)])
    print('\nMean Intensity of Image = {} from ROI = {}\nMean= {}\nStd= {}\nNumber of Voxels= {}\n'.format(args.img, args.mask, mean, std, len(img[np.where(mask == 1)])))
