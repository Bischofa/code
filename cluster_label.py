#!/usr/bin/env python

import numpy as np
import os
import nibabel as nb
from collections import defaultdict
import argparse

def has_adjacent(i,j):
	r_list=[]
	try:
		if im[i-1,j+1]!=0 and tabl[(i-1,j+1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i-1,j+1))
	except:
		pass
		
	try:
		if im[i,j+1]!=0 and tabl[(i,j+1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i,j+1))
	except:
		pass
		
	try:
		if im[i+1,j+1]!=0 and tabl[(i+1,j+1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i+1,j+1))
	except:
		pass
		
	try:
		if im[i+1,j]!=0 and tabl[(i+1,j)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i+1,j))
	except:
		pass
	try:
		if im[i+1,j-1]!=0 and tabl[(i+1,j-1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i+1,j-1))
	except:
		pass
	try:
		if im[i,j-1]!=0 and tabl[(i,j-1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i,j-1))
	except:
		pass
	try:
		if im[i-1,j-1]!=0 and tabl[(i-1,j-1)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i-1,j-1))
	except:
	        pass
	try:
		if im[i-1,j]!=0 and tabl[(i-1,j)]==0 and abs(i)==i and abs(j)==j:
			r_list.append((i-1,j))
	except:
		pass
	return r_list
def cluster_find(i,j,num):
	#recursive function
	global tabl
	#print('first {}:{}'.format(i,j))
	if  im[i,j]==0 or tabl[(i,j)] or abs(i)!=i or abs(j)!=j:
		return
	else:
		#print('{}:{}'.format(i,j))
		tabl[(i,j)]=num
		for kernel in has_adjacent(i,j):
			cluster_find(kernel[0],kernel[1],num)
def lesion_fill(cluster_num,im):
	for i in cluster_num:
		
		lesion=sorted(cluster_num[i],key=lambda x:x[0])
		lesion=sorted(lesion,key=lambda x:x[1])
		line_dict=defaultdict(list)
		for j in lesion:
			line_dict[j[1]].append(j)
		for k in line_dict:
			#print(line_dict[k])
			#input('line_dict')
			for filler in range(line_dict[k][-1][0]-line_dict[k][0][0]):
				#print(filler)
				im[filler+line_dict[k][0][0],k]=1
	return im

#2d lesion finding/labeling and filling
def main_func(img_path):
	
	print('\nRunning cluster label...')	
	image=nb.load(img_path)
	imall=image.get_data()
	sh=imall.shape
	num=0
	global tabl
	tabl=defaultdict(int)        
	global im
	for x in range(sh[2]):
		tabl=defaultdict(int)
		im=imall[:,:,x]	
		for i in range(sh[0]):
			for j in range(sh[1]):
				if im[i,j]==0 or tabl[(i,j)]:
					continue
				num=num+1
				#print('{}:{}'.format(i,j))
				cluster_find(i,j,num)
		clusters=defaultdict(int)
		cluster_num=defaultdict(list)
		for coords in tabl:
			clusters[tabl[coords]]+=1
			cluster_num[tabl[coords]].append(coords)
			#im[coords[0],coords[1]]=tabl[coords]
		imall[:,:,x]=lesion_fill(cluster_num,im)

	return imall, image.affine	

if __name__ == '__main__':
    parser = argparse.ArgumentParser('This code fills a jim roi that is converted into a nifti file into a filled binary mask - useful for lesion filling and segmentation filling')
    parser.add_argument('img', help='path of nifti file to be filled')
    parser.add_argument('out', help='filename of output filled lesion mask')
    args = parser.parse_args()
    img, imgaff = cluster_label(args.img)
    nb.save(nb.Nifti1Image(img,imgaff),args.out)
