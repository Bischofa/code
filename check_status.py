#!/usr/bin/env python
from subprocess import check_output
from glob import glob
import pandas as pd
import argparse
import csv
import os

#usage: check_status.py <filename csv spreadsheet with accession numbers>

# get command line input
parser = argparse.ArgumentParser()
parser.add_argument("csv", help="csv spreadsheet with ACCESSION column")
args = parser.parse_args()

# initiate output spreadsheet
writer = open(args.csv[:-4]+"_status.csv", "w")
spreadsheet = csv.DictWriter(writer, fieldnames=["Accession", "upload", "msID", "dicoms", "niftis", "freesurfer"])
spreadsheet.writeheader()

# load input spreadsheet
df = pd.read_csv(args.csv)

# processing
for i in df.index:
    row = {}
    row["Accession"] = df.ACCESSION[i]
    output = check_output(["ms_get_exam_id", "-a", str(df.ACCESSION[i])])
    for line in output.split("\n"):
        if "mse_number:" in line:
            mse = "mse%s" % line.split()[-1]
            row["upload"] = mse
	if "msID:" in line:
            msID = "msID%s" % line.split() [-1]
            row["msID"] = msID
    if glob("/data/henry6/PBR/dicoms/"+mse):
        row["dicoms"] = "x"
    if glob("/data/henry[78]/PBR/subjects/"+mse):
        row["niftis"] = "x"
    if glob("/data/henry6/PBR/surfaces/*%s*" %mse):
        row["freesurfer"] = "x"
        row["msID"] = os.path.basename(glob("/data/henry6/PBR/surfaces/*%s*" %mse)[0]).split("-")[0]
    spreadsheet.writerow(row)

writer.close()
