from nipype import config
import nibabel as nib
import nipype.pipeline.engine as pe
from nipype import Workflow, Node
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.utility import Function
import os
import argparse
from glob import glob
import nibabel as nib
import numpy as np
from pbr.base import _get_output



def pca_data(img):
    #This node takes in the FIRST segmented branstem and then runs sklearn PCA on the segmentation to find the axis along which there is greatest variance in the coordinates of the segmentation. I.E, the main axis along the brainstem.
 
    import nibabel as nib
    import numpy as np
    from sklearn.decomposition import PCA
    from math import pi
    
    # PCA FIT
    img_data = nib.load(img).get_data()
    bs = np.where(img_data != 0)
    X = list(zip(bs[0], bs[1], bs[2]))
    pca = PCA(n_components=3)
    pca.fit(X)
    print('\npca variance', pca.explained_variance_ratio_)
    print('pca componenets', pca.components_)
    component = pca.components_[np.where(pca.explained_variance_ratio_ == np.max(pca.explained_variance_ratio_))][0]
    component[2] = abs(component[2]) #have to correct the z component to be positive
    print('\nselected component', component)
    center = pca.mean_
    
    # find cross product between the z-axis [0,0,1] and the PCA brainstem axis to find the axis of rotation
    X_pdct = np.cross(component, [0, 0,1])
    print('\ncross product of vectors', X_pdct)
    
    # find the angle between the PCA axis and the z-axis
    # Got 3d vector angle from this website http://www.analyzemath.com/stepbystep_mathworksheets/vectors/vector3D_angle.html
    theta = np.arccos(np.dot(component,[0,0,1])/(np.linalg.norm([0, 0,1]) * np.linalg.norm(component)))/pi * 360
    print('Theta', theta)
    
    return center, X_pdct, theta




def makerot(theta, in_img, out_mat, axis, center=[]):
    #wrapper fsl makerot, makes a rotation matrix if you give an angle, a center of rotation and an axis
    from subprocess import check_call
    
    cmd = ['makerot', '--theta={}'.format(theta), '--cov={}'.format(in_img), 
       '--axis={},{},{}'.format(axis[0], axis[1], axis[2]), 
        '--out={}'.format(out_mat)]
    if len(center) == 3:
        cmd.append('--centre={},{},{}'.format(center[0], center[1], center[2]))
    print(cmd)
    check_call(cmd)

    return out_mat, theta

def run_flirt(in_img, ref, mat, out_img, theta=0, degrees=False):
    #wrapper for fsl flirt

    from subprocess import check_call
    def print_shape(img):
        import nibabel as nib
        print('\n{} shape\n{}\n'.format(img, nib.load(img).get_data().shape))
    print_shape(in_img)    

    if degrees == True:
        out_img = out_img.split('.nii')[0]+'_{}_d.nii.gz'.format(theta)
    
    cmd = ['flirt', '-in', in_img, '-ref', ref, '-applyxfm', '-init', mat, '-out', out_img]
    print(cmd)
    check_call(cmd)
    return out_img

def JIMMY(img_pth):
    #Allows user to specify the location of the initial axis of rotation    
    from subprocess import Popen 
    import nibabel as nib

    #get img header
    zooms = nib.load(img_pth).header.get_zooms()[:3]
    
    # open img with Jim
    cmd = ['/data/henry2/arajesh/py_scripts/henrylab_utils/Jimbo', img_pth]
    print('\n*************\nTHE INPUT FOR THE COORDS WILL COME AFTER YOU CLOSE JIM, MAKE SURE YOU REMEMBER X,Y,Z COORDS OF MOST CRANIAL MARKER!!!\n**************')
    print('\nOpening {}'.format(img_pth))
    print('\nvoxel dimensions of img: {}'.format(zooms))
    proc = Popen(cmd)
    
        
    # ask for user input of x,y,z coords for most cranial marker placement
    jim_coord = str(input('Jim coords for the most cranial marker of the cord!\nPlease use commas to parse like this: 3.44,12.1,3  for X,Y,Z coords\nAlso you only need to enter to the .001 decimal place for marker placement\nType here: ')).split(',')
    jim_coord = [float(x) for x in jim_coord]
    
    img_data, img_aff = nib.load(img_pth).get_data(), nib.load(img_pth).affine
    x_dim, y_dim = img_data.shape[:-1]
    
    
    # convert from jim coords to nifti coords
    nii_x = int(x_dim/2.0 * (1.0+jim_coord[0]/100.0))
    nii_y = int(y_dim/2.0 * (1.0-jim_coord[1]/100.0))
    center = [nii_x, nii_y, int(jim_coord[2])]
    
    abs_dist = 5 #ABS DIST IS THE DISTANCE FROM THE MARKER INPUT Z COORD AND THE CUTOFF FOR THE IMG
    z_vox = nib.load(img_pth).header.get_zooms()[2]
    rel_dist = abs_dist/z_vox
    slice_pos = int(jim_coord[2] - rel_dist) #slice_pos is the the cutoff for the image when we resample it

    #center is the center of the rotation that is returned for the cord, it is 2mm down from the most cranial marker placement
    center = (nii_x, nii_y, int(jim_coord[2] - 2/z_vox))
    print('CENTER WHICH WE ARE GOING TO ROTATE AROUND', center)

    #Make centering translation matrix
    trans_aff = '/'.join(img_pth.split('/')[:-1])+'/translation.mat'
    f = open(trans_aff, 'w')
    f.write('1  0  0  {}\n'.format(int((img_data.shape[0]/2-center[0])/2)))
    f.write('0  1  0  {}\n'.format(int((img_data.shape[1]/2-center[1])/2)))
    f.write('0  0  1  {}\n'.format(int((img_data.shape[2]/2-center[2])/2)))
    f.write('0  0  0  1')
    f.close()


    #Make centering translation inverse matrix
    inv_trans_aff = '/'.join(img_pth.split('/')[:-1])+'/inv_translation.mat'
    f = open(inv_trans_aff, 'w')
    f.write('1  0  0  {}\n'.format(-int((img_data.shape[0]/2-center[0])/2)))
    f.write('0  1  0  {}\n'.format(-int((img_data.shape[1]/2-center[1])/2)))
    f.write('0  0  1  {}\n'.format(-int((img_data.shape[2]/2-center[2])/2)))
    f.write('0  0  0  1')
    f.close()

    out_img = img_pth

    return out_img, center, slice_pos, inv_trans_aff, trans_aff


def open_jim(img_pth):
    #not currently used. 
    from subprocess import check_call
    
    cmd = ['/data/henry2/arajesh/py_scripts/henrylab_utils/Jimbo', img_pth]
    print('\nOpening {}'.format(img_pth))
    check_call(cmd)
    slice_pos = int(input('\nWhat is the z axis location of the pontomedullary junction: '))
    return img_pth, slice_pos


def resample(img_pth):
    #resamples img to 1mm^3 voxel size 
    def print_shape(img):
        import nibabel as nib
        print('\n{} shape\n{}\n'.format(img, nib.load(img).get_data().shape))
    print_shape(img_pth)
    import nibabel as nib
    from dipy.align.reslice import reslice
    
    img = nib.load(img_pth)
    zooms = img.header.get_zooms()[:3]
    new_zooms = [1.0, 1.0, 1.0] #size of the new voxel dimensions
    newdata, newaff = reslice(img.get_data(), img.affine, zooms, new_zooms)

    #out_img = '/'.join(img_pth.split('/')[:-1]) + '/img_resample.nii.gz'
    out_img = img_pth
    nib.save(nib.Nifti1Image(newdata, newaff), out_img) 

    return out_img 


def crop(slice_pos, img_pth):
    #crops an image from a given slice_pos in the z-dimension
    import nibabel as nib
    def print_shape(img):
        import nibabel as nib
        print('\n{} shape\n{}\n'.format(img, nib.load(img).get_data().shape))
    print_shape(img_pth)
#    import numpy as np
#    import math
#    abs_dist = 25.0 # in mm
#    z_vox = nib.load(img_pth).header.get_zooms()[2]
#    rel_dist = abs_dist/z_vox
#    angle = abs(float(img_pth.split('/')[-1].split('_')[2]))
#    crop_dis = int(rel_dist * math.cos(math.radians(angle))) # rounding to closest integer
    
    img = nib.load(img_pth)
    img_data, img_aff = img.get_data(), img.affine
   
    #out_img =  '/'.join(img_pth.split('/')[:-1]) + '/img_cropped.nii.gz'
    out_img = img_pth
    nib.save(nib.Nifti1Image(img_data[:,:,slice_pos:], img_aff), out_img )
    
    return out_img
    
def reshape_img(img_pth, grow=True):

    import nibabel as nib
    import numpy as np
    from os.path import join as opj
    def print_shape(img):
        import nibabel as nib
        print('\n{} shape\n{}\n'.format(img, nib.load(img).get_data().shape))
    print_shape(img_pth)   
    img = nib.load(img_pth)
    img_data = img.get_data()
    growth = 100 #factor to grow the image by

    if grow == True:
        print('\n ********************* \n GROW true\n***************************\n')
        new_img = np.zeros([x+2*growth for x in img_data.shape])
        new_img[growth:img_data.shape[0]+growth, growth:img_data.shape[1]+growth, growth:img_data.shape[2]+growth] = img_data
        out_img = img_pth
        nib.save(nib.Nifti1Image(new_img, img.affine), out_img) 

    elif grow == False:
        print('\n ********************* \n GROW FALSE\n***************************\n')
        new_img = img_data[growth:img_data.shape[0]-growth, growth:img_data.shape[1]-growth, growth:img_data.shape[2]-growth]
        out_img = img_pth
        #out_img = '/'.join(img_pth.split('/')[:-1]) + '/img_reshaped.nii.gz'
        nib.save(nib.Nifti1Image(new_img, img.affine), out_img) 

    return out_img 


def Rotate_wf(bstem_reo, T1_reo, mseid):
    
    import nipype.interfaces.io as nio
    import nipype.interfaces.fsl as fsl 

    theta_array = np.arange(0, 38.5, 3.5) #this line changes the angles of rotation 
    print('\n*******\nThe angles you are rotating arounda are {}\n*******\n'.format(theta_array))
 
    working_dir = '/working/henry_temp/keshavan'
    wf = Workflow(name='{}_rotate'.format(mseid), base_dir=working_dir)

    inputspec = Node(IdentityInterface(fields=['brstm_reo', 't1_reo']),
                    name='inputspec')
    
    inputspec.inputs.brstm_reo = bstem_reo
    inputspec.inputs.t1_reo = T1_reo

    pca_node = Node(name='pca_node',
                   interface=Function(input_names=['img'],
                                     output_names=['center', 'X_pdct', 'theta'],
                                     function=pca_data))
                  
    wf.connect(inputspec, 'brstm_reo', pca_node, 'img')

    xyz_search = Node(name='xyz_search',
                      interface=Function(input_names=['img_pth'],
                                         output_names=['out_img', 'center', 'slice_pos', 'inv_trans_aff', 'trans_aff'],
                                         function=JIMMY))


    rot = Node(name='rot',
                 interface=Function(input_names=['theta', 'in_img', 'out_mat', 'axis', 'center'],
                                   output_names=['out_mat', 'theta'],
                                   function=makerot))
    
    rot.inputs.out_mat = '/'.join(bstem_reo.split('/')[:-1]) + '/{}_rot.mat'.format(mseid)
    wf.connect(pca_node, 'center', rot, 'center')
    wf.connect(pca_node, 'theta', rot, 'theta')
    wf.connect(pca_node, 'X_pdct', rot, 'axis')
    wf.connect(inputspec, 't1_reo', rot, 'in_img')
    
    flirt = Node(name='flirt',
               interface=Function(input_names=['in_img', 'ref', 'mat', 'out_img'],
                                 output_names=['out_img'],
                                 function=run_flirt))
    
    wf.connect(inputspec, 't1_reo', flirt, 'in_img')
    wf.connect(inputspec, 't1_reo', flirt, 'ref')
    wf.connect(rot, 'out_mat', flirt, 'mat')
    flirt.inputs.out_img = '/'.join(bstem_reo.split('/')[:-1]) + '/{}_PCA_align.nii.gz'.format(mseid)

    wf.connect(flirt, 'out_img', xyz_search, 'img_pth')

    img_reshape1 = Node(name='img_reshape',
                       interface=Function(input_names=['img_pth'],
                                          output_names=['out_img'],
                                          function=reshape_img))

    wf.connect(xyz_search, 'out_img', img_reshape1, 'img_pth')

    translate = Node(name='translate',
                     interface=Function(input_names=['in_img', 'ref', 'mat', 'out_img'],
                                        output_names=['out_img'],
                                        function=run_flirt))

    wf.connect(img_reshape1, 'out_img', translate, 'in_img')
    wf.connect(img_reshape1, 'out_img', translate, 'ref')
    wf.connect(xyz_search, 'trans_aff', translate, 'mat')
    translate.inputs.out_img = '/'.join(bstem_reo.split('/')[:-1]) + '/{}_PCA_align_transl.nii.gz'.format(mseid)

    small_rot = Node(name='small_rot',
             interface=Function(input_names=['theta', 'in_img', 'out_mat', 'axis', 'center'],
                              output_names=['out_mat', 'theta'],
                               function=makerot))
    small_rot.iterables = ('theta', theta_array) 
    small_rot.inputs.axis = [1,0,0]
    small_rot.inputs.out_mat = '/'.join(bstem_reo.split('/')[:-1]) + '/{}_rot_M.mat'.format(mseid)

    wf.connect(translate, 'out_img', small_rot, 'in_img')
    wf.connect(xyz_search, 'center', small_rot, 'center')
    
    m_flirt =  Node(name='m_flirt',
                   interface=Function(input_names=['in_img', 'ref', 'mat', 'out_img', 'degrees', 'theta'],
                                     output_names=['out_img'],
                                     function=run_flirt))
    m_flirt.inputs.out_img = '/'.join(T1_reo.split('/')[:-1]) +  '/{}_PCA.nii.gz'.format(mseid)
    m_flirt.inputs.degrees = True
    wf.connect(small_rot, 'theta', m_flirt, 'theta')
    wf.connect(small_rot, 'out_mat', m_flirt, 'mat')
    wf.connect(translate, 'out_img', m_flirt, 'in_img')
    wf.connect(translate, 'out_img', m_flirt, 'ref')


    inv_translate = Node(name='inv_translate',
                     interface=Function(input_names=['in_img', 'ref', 'mat', 'out_img'],
                                        output_names=['out_img'],
                                        function=run_flirt))

    wf.connect(m_flirt, 'out_img', inv_translate, 'in_img')
    wf.connect(m_flirt, 'out_img', inv_translate, 'ref')
    wf.connect(xyz_search, 'inv_trans_aff', inv_translate, 'mat')
    wf.connect(m_flirt, 'out_img', inv_translate, 'out_img')    

    img_reshape2 = Node(name='img_reshape2',
                       interface=Function(input_names=['img_pth', 'grow'],
                                          output_names=['out_img'],
                                          function=reshape_img))
    img_reshape2.inputs.grow= False #shrinks the image back to small size
    wf.connect(inv_translate, 'out_img', img_reshape2, 'img_pth') 
       
    cropper = Node(name='cropper',
               interface=Function(input_names=['slice_pos', 'img_pth'],
                                 output_names=['out_img'],
                                 function=crop))

    wf.connect(img_reshape2, 'out_img', cropper, 'img_pth')
    wf.connect(xyz_search, 'slice_pos', cropper, 'slice_pos')
    
    resamp = Node(name='resamp',
                 interface=Function(input_names=['img_pth'],
                                   output_names=['out_img'],
                                   function=resample))
    
    wf.connect(cropper, 'out_img', resamp, 'img_pth')   

#    FMresamp = Node(name='FMresamp',
#                 interface=Function(input_names=['img_pth'],
#                                   output_names=['out_img'],
#                                   function=resample))
#
#    wf.connect(img_reshape2, 'out_img', FMresamp, 'img_pth')
 
#    sinker = pe.Node(nio.DataSink(), name='sinker', try_hard_link_dataset=False,
#                     run_without_submitting = True)
#    sinker.inputs.base_directory = '/'.join(T1_reo.split('/')[:-1])
#    print('\nBASE DIR {}'.format('/'.join(T1_reo.split('/')[:-1])))
#    wf.connect(resamp, 'out_img', sinker, '@resamp')
    #wf.connect(FMresamp, 'out_img', sinker, 'FM')
    #wf.write_graph(graph2use='orig', simple_form=False)
    wf.run()

def main(mseid):
    from glob import glob

    for t1 in glob(os.path.join(_get_output(mseid), mseid, 'C1_auto', '*', '{}_BrStm_reo.nii.gz'.format(mseid))):
        bs_reo = '/'.join(t1.split('/')[:-1]) + '/{}_bs_reo.nii.gz'.format(mseid) 
        print('\nt1', t1, 'bs_reo' ,bs_reo)
        Rotate_wf( t1, bs_reo, mseid)    

if __name__ == '__main__':

    parser = argparse.ArgumentParser('This script aligns the images for you')
    parser.add_argument('mseid', help='which mseid to run on')
    args = parser.parse_args()
    main(args.mseid)
