import csv
import os

fileCount = 0

directory = r'C:\Users\patmcc\Desktop\CSV_CLEANER'

if not os.path.exists(directory + '\cleaned'):
    os.makedirs(directory + '\cleaned')

for root,dirs,files in os.walk(directory):
    for file in files:
       if file.endswith(".csv"):
            with open(directory + '\\'+ file, 'r') as inp, open(directory + '\\cleaned' + '\\'+ file[:-4] + '_LC.csv', 'w') as out:
                writer = csv.writer(out, lineterminator = '\n')
                for row in csv.reader(inp):
                    if row[3] != '':
                        writer.writerow(row)
            with open(directory + '\\'+ file, 'r') as inp, open(directory + '\\cleaned' + '\\'+ file[:-4] + '_ACCEL.csv', 'w') as out:
                writer = csv.writer(out, lineterminator = '\n')
                for row in csv.reader(inp):
                    if row[15] != '':
                        writer.writerow(row)
       fileCount += 1
    break  # prevent decending into subfolders
print('Cleaned ' + str(fileCount) + ' files.')