#!/usr/bin/env python
from subprocess import call, check_call
from getpass import getpass
import csv
import os
import sys
import argparse

from nipype.utils.filemanip import load_json, save_json

## load study keys ##
study_to_directory = load_json('/data/henry6/alyssa/scripts/studies.json')
keys = ""
for i, key in enumerate(sorted(study_to_directory.keys())):
    if i != (len(study_to_directory.keys()) - 1):
        keys += "{}, ".format(key)
    else:
        keys += key
# controls currently being sent to NOT YET SCREENED folder
# what is study_opera2 ???

## help ##
usage = """ms_relocate_scans.py <csv spreadsheet>"""
parser = argparse.ArgumentParser(description="ms_relocate_scans.py grabs directory and ACC# from csv spreadsheet and then copies the data accordingly. Included studies are %s" %keys)
parser.add_argument("spreadsheet", help="csv with first column titled STUDY and the second titled ACCESSION")

if len(sys.argv) != 2:
    parser.print_help()
    sys.exit(1)
else:
    args = parser.parse_args()

## let's get started ##
password = getpass("password: ")

data = csv.DictReader(open(args.spreadsheet, 'r'))
# fns = csv.fieldnames

for row in data:
    try:
        directory = study_to_directory[row['STUDY']]
    except KeyError:
        new_study = raw_input("Is %s a new study?"
                              "\n(Reminder: new study tags need to be initialized in the DB)\n"
                              "yes/no: " % row["STUDY"])
        if (new_study == "yes") or (new_study == "y"):
            directory = raw_input("What is the study directory path? \n")
            print("Thanks! Adding this new study to the dictionary...")
            if not os.path.exists(directory):
                os.mkdir(directory)
            study_to_directory[row["STUDY"]] = directory
            save_json("/data/henry6/alyssa/scripts/studies.json", study_to_directory)
        elif (new_study == "no") or (new_study == "n"):
            print("Your input (%s) was not found in our dictionary. What should the study name have been?"
                  % row["STUDY"])
            print(keys)
            row['STUDY'] = raw_input("What should the study name have been?"
                                     "\n(Reminder: study keys are case sensitive)\n"
                                     "Study name: ")
            directory = study_to_directory[row['STUDY']]
    print("changing directory to %s" % directory)
    os.chdir(directory)
    print("copying exam to study directory")
    call(['ms_get_scanner_data', '-a', row['ACCESSION'], '--study', row['STUDY'], '--all_raw', '-p', password]) #previously check_call


