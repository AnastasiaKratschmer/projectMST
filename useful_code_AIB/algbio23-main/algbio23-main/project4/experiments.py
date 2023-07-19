import csv

from rfdist import rfdist_from_files
import os
from fnmatch import fnmatch

files_normal = []
files_permuted = []
for file in os.listdir('experimentsdata'):
    if fnmatch(file, '*permuted'):
        files_permuted.append('experimentsdata/' + file)
    elif fnmatch(file, 'kalign*') or fnmatch(file, 'CLUSTAL*') or fnmatch(file, 'MUSCLE*') and not fnmatch(file,'*permuted'):
        files_normal.append('experimentsdata/' + file)
files_normal.sort()
files_permuted.sort()

res_normal = [[0] * 6 for i in range(0,6)]
res_permuted = [[0] * 6 for i in range(0,6)]
for i in range(0,6):
    for j in range(0, 6):
        res_normal[i][j] = rfdist_from_files(files_normal[i], files_normal[j])
        res_permuted[i][j] = rfdist_from_files(files_permuted[i], files_permuted[j])
res_3 = [0 for i in range(0,6)]
for k in range(0,6):
    res_3[k] = rfdist_from_files(files_normal[k], files_permuted[k])

with open('normaldist', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(files_normal)
    writer.writerows(res_normal)

with open('permuteddist', 'w') as fp:
    writer = csv.writer(fp)
    writer.writerow(files_permuted)
    writer.writerows(res_permuted)

with open('exp3', 'w') as fh:
    writer = csv.writer(fh)
    writer.writerow(files_normal)
    writer.writerow(res_3)