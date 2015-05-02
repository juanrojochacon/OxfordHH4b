#!/bin/python

import os
import sys
import math
regimes = ['res', 'inter','boost']
folder = sys.argv[1]
#folder = "oxford_combined/signal"
print folder

recursive_folders = ['signal', 'background', 'SHERPA_QCD4b', 'SHERPA_QCD2b2j', 'SHERPA_QCD4j']
file_types = ['res', 'inter', 'boost']

new_data = {}

for type in file_types:
    new_data[type] = []

for k, recursive_folder in enumerate(recursive_folders):

    for j, file_type in enumerate(file_types):

        filename = 'histo_CF_' + file_type + '.dat'
        filePath = os.path.join(folder,recursive_folder,filename)

        with open(filePath) as f:

            Cutflow = f.readlines()[4:8]

            for l, line in enumerate(Cutflow):
                if k == 0:
                    new_line = []
                    split_line = line.split('\t')
                    new_line.append(float(split_line[0]))
                    new_line.append(split_line[2])           
                    new_data[file_type].append(new_line)
                else:
                    new_data[file_type][l].append(line.split('\t')[2])

for file_type in file_types:
    list_of_rows = new_data[file_type]
    #print list_of_rows[0][0]
    outfile = 'cutflow_' + folder + '_' + file_type + '.dat'
    with open(outfile, 'w') as f:
        # generate header line
        f.write(file_type+'\t')
        for recursive_folder in recursive_folders:
            f.write(recursive_folder+'\t')
        f.write('\n')
        for h, row in enumerate(list_of_rows):
            for i, item in enumerate(row):
                if i == 0:
                    f.write('C%d' % int(item))
                else:
                    f.write(str(item))
                    if h>0:
                        #f.write(' ')
                        #f.write(str(list_of_rows[h-1][i]))
                        if float(list_of_rows[h-1][i]) != 0:
                            cuteff = 100*float(item)/float(list_of_rows[h-1][i])
                        else:
                            cuteff = 0
                        f.write(' ')
                        #f.write(str(cuteff))
                        f.write('(%.2f' % cuteff)
                        f.write('%)')
                        #cuteff = float(item)/float(list_of_rows[h-1][i])
                        #f.write(str(cuteff) )
                if i == len(row) - 1:
                    f.write('\n')
                else:
                    f.write('\t')


summary_outfile = 'cutflow_' + folder + '_summary.dat'
with open(summary_outfile, 'w') as f:
    for i in range(0,4):
        #print i
        # header line
        if i == 0 or i == 1:
            f.write(recursive_folders[i]+'\t')
        elif i == 2:
            f.write('S/B'+'\t')
        elif i == 3:
            f.write('S/sqrt(B)'+'\t')
        for file_type in file_types:
            f.write(file_type+'\t')
        f.write('\n')
        # write main lines
        for j in range(len(new_data["res"])):
            for k, file_type in enumerate(file_types):
                if k == 0:
                    f.write('C%d' % (new_data[file_type][j][0])+'\t')
                # switch for tables, first sig+bkg, then S/B and sqrt.
                if i == 0 or i==1:
                    f.write(new_data[file_type][j][i+1]+'\t')
                elif i == 2:
                    if float(new_data[file_type][j][2]) == 0.0:
                        SoverB = 0
                    else:
                        SoverB = float(new_data[file_type][j][1])/float(new_data[file_type][j][2])
                    f.write('%e'% SoverB + '\t')
                elif i == 3:
                    if float(new_data[file_type][j][2]) == 0.0:
                        SoverRootB = 0
                    else:
                        SoverRootB = float(new_data[file_type][j][1])/math.sqrt(float(new_data[file_type][j][2]))
                    f.write('%e'% SoverRootB + '\t')
            f.write('\n')
        f.write('\n\n')

    
    
