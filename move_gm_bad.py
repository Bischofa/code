import os
import pandas as pd
import glob as glob
from shutil import move

working_directory = '/data/henry6/antje/PMD/subjects'
df = pd.read_excel('/data/henry6/antje/PMD/scripts/PMD_GMA.xlsx')

for mse in df['not_reliable'].values:
    mseid = 'mse' + str(mse)
    src = os.path.join(working_directory, mseid)
    print(src)
    dst = os.path.join(working_directory, 'unreliable', mseid)
    print(dst)
    move(src, dst)
