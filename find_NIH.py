#!/usr/bin/env python
from subprocess import Popen, PIPE, check_call
from getpass import getpass
from glob import glob
import argparse
import dicom
import csv
import os




def main(path, sequence):
    exams = sorted(glob(path+"/*/E*"))
    if exams:
        for exam in exams:
            if os.path.islink(os.path.dirname(exam)):
                continue
            firsts = sorted(glob(exam+"/[0-9]*/*I1.DCM*"))
            if firsts:
                a = 0
                b = 0
                for first in firsts:
                    try:
                        if "gz" in first:
                            check_call(["gunzip", first])
                            first = first[:-3]
                        dcm = dicom.read_file(first)
                    except:
                        break
                    series = dcm.SeriesDescription
                    if sequence in series:
                        b = 1
                    if ("NIH FLAIR" in series):
                        a = 1

                    if a and b:
                        row = {}
                        row["msid"] = os.path.basename(path)
                        row["name"], row["date"] = get_info(dcm)
                        row["path"] = exam
                        spreadsheet.writerow(row)
                        break

def get_info(dcm):
    name = dcm.PatientName
    if name[:2] == "ms":
        proc = Popen(["ms_get_phi", "--examID", dcm.AccessionNumber, "-p", password], stdout=PIPE)
        for line in proc.stdout.readlines():
            line = line.strip()
            if "PatientName" in line:
                name = line.split(" = ")[-1]
            if "StudyDate" in line:
                date = line.split(" = ")[-1]
    else:
        date = dcm.StudyDate

    return name, date

if __name__ == "__main__":

    # add command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", help="the keyword for the sequence you'd like to search for")
    parser.add_argument("spreadsheet", help="spreadsheet name to write to")
    parser.add_argument("paths", nargs="+", help="paths to msid folder; e.g. /data/henry[0-9]/EPIC*/ms*")
    args = parser.parse_args()
    
    # initiate spreadsheet
    to_write = open(args.spreadsheet, "wb")
    spreadsheet = csv.DictWriter(to_write, 
                                 fieldnames=["msid", "name", "date", "path"])
    spreadsheet.writeheader()

    # connecting to database
    password = getpass("mspacman password: ")
    check_call(["ms_dcm_echo", "-p", password])

    # searching data
    for i, path in enumerate(args.paths):
        print "querying %i/%i" %(i+1, len(args.paths))
        if os.path.islink(os.path.dirname(path)):
            continue
        main(path, args.sequence)

    # finishing up
    to_write.close()
