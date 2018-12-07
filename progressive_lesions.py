#!/usr/bin/env python
from subprocess import check_call
import numpy as np
import os
try:
    import pandas as pd
except:
    os.environ["PYTHONPATH"] = "/data/henry7/software/"
    import pandas as pd
import argparse


def get_thr(date):
    if date < 20050402:
        threshold = 300
    elif (date >= 20090301) and (date <= 20091001):
        threshold = 435
    elif (date > 20091001) or (np.isnan(date)):
        threshold = 424
    else:
        threshold = 400
    return threshold

def find_previous(tp, df):
    checking = range(1, int(tp))[::-1]
    for c in checking:
        prev = str(c)
        if not np.isnan(df[prev].values[0]):
            break
    return prev

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument("-msid", metavar="ms####",
                        help="msid (Example: ms0056)")
    parser.add_argument("-tp", metavar="yr#",
                        help="time point (Example: yr1)")
    args = parser.parse_args()

    subject = args.msid
    tp = str(int(args.tp[-1]) + 1)

    dates = "/data/henry6/alyssa/ratio_hist_norm/cb_scan_dates.csv"
    df = pd.read_csv(dates)
    df = df[df.msid == subject]
    logfile = open("processing_log.txt", "w")

    ratio_paths = {"6": "/data/henry6/esha/Old_Files/T2T1_Ratio_Zscores/Yr5/calculating_ratios/"+subject+"_yr5_ratio_brain.nii.gz",
                   "5": "/data/henry6/esha/Old_Files/T2T1_Ratio_Zscores/Yr4/calculating_ratios/"+subject+"_yr4_ratio_brain.nii.gz",
                   "4": "/data/henry6/esha/Old_Files/T2T1_Ratio_Zscores/Yr3/calculating_ratios/"+subject+"_yr3_ratio_brain.nii.gz",
                   "3": '/data/henry6/esha/Old_Files/EPIC_Yr2/LST_files/ratio/'+subject+'_yr2_ratio_brain.nii.gz',
                   "2": '/data/henry6/esha/Old_Files/T2T1_Ratio_Zscores/Yr1/calculating_ratios/'+subject+'_yr1_ratio_brain.nii.gz',
                   "1": '/data/henry6/esha/Old_Files/EPIC_Baseline/LST_files/ratio_brain/'+subject+'_yr0_ratio_brain.nii.gz'
                  }

    lesion_paths = {"1": '/data/henry6/alyssa/lesions_reg/manual_lesions/baseline/'+subject+'_yr0_lesions.nii.gz',
                    "2": "/data/henry8/alyssa/lesions_reg/get_yr1_lesions/"+subject+"_yr1_lesions.nii.gz",
                    "3": "/data/henry8/alyssa/lesions_reg/get_yr2_lesions/"+subject+"_yr2_lesions.nii.gz",
                    "4": "/data/henry8/alyssa/lesions_reg/get_yr3_lesions/"+subject+"_yr3_lesions.nii.gz",
                    "5": "/data/henry8/alyssa/lesions_reg/get_yr4_lesions/"+subject+"_yr4_lesions.nii.gz",
                    "6": "/data/henry8/alyssa/lesions_reg/get_yr5_lesions/"+subject+"_yr5_lesions.nii.gz"
                   }

    prev = find_previous(tp, df)
    scandate = df[tp].values[0]
    prevdate = df[prev].values[0]
    threshold = get_thr(scandate)
    print "\n\nThe threshold that should be used for %s %s is %i\n\n" %(subject, args.tp, threshold)

    if not os.path.exists(ratio_paths[prev]):
        print "Couldn't find the previous ratio"
    elif not os.path.exists(ratio_paths[tp]):
        print "Couldn't find the current ratio"
    else:
        prev_lesions = os.path.join(os.path.dirname(lesion_paths[tp]), 
                                    subject+'_prevTO%s_lesions.nii.gz' %args.tp)
        prev_ratio = os.path.join(os.path.dirname(lesion_paths[tp]), 
                                  subject+'_prevTO%s_ratio_brain.nii.gz' %args.tp)
        if args.tp in ["yr2", "yr5"]:
            manual = "/data/henry6/alyssa/lesions_reg/manual_lesions"
            if args.tp == "yr2":
                orig = os.path.join(manual, "yr2", 
                                      "_".join([subject, args.tp, 
                                                "lesions.nii.gz"]))
                lesion = os.path.join(manual, "yr2", 
                                      "_".join([subject, args.tp, 
                                                "lesions1.nii.gz"]))
                check_call(["cp", orig, lesion])
            else:
                orig = os.path.join(manual, "yr5_flair", 
                                      "_".join([subject, args.tp, 
                                                "lesions.nii.gz"]))
                lesion = os.path.join(manual, "yr5_flair", 
                                      "_".join([subject, args.tp, 
                                                "lesions1.nii.gz"]))
                check_call(["cp", orig, lesion])
            check_call(["freeview", "-v", prev_ratio+":grayscale=0,1000", 
                        "-v", ratio_paths[tp]+":grayscale=0,1000",
                        "-v", prev_lesions+":colormap=heat",
                        "-v", lesion+":colormap=jet:opacity=0.5"])
            check_call(["mv", lesion, lesion_paths[tp]])
        else:
            check_call(["freeview", "-v", prev_ratio+":grayscale=0,1000", 
                        "-v", ratio_paths[tp]+":grayscale=0,1000",
                        "-v", prev_lesions+":colormap=heat",
                        "-v", lesion_paths[tp]+":colormap=jet:opacity=0.5"])
