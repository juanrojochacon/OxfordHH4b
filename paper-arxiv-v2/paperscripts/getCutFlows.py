#!/usr/bin/python2.7

import os
import sys
import math
regimes = ['res', 'inter','boost']
folder = "../../paper/plotdata/results_noPU/"
print folder

recursive_folders = ['diHiggs', 'background', 'SHERPA_QCD4b', 'SHERPA_QCD2b2j', 'SHERPA_QCD4j', 'SHERPA_QCDttbar']
file_types = ['res', 'inter', 'boost']

Cnum = ['C0', 'C1a', 'C1b', 'C1c', 'C1d', 'C1e', 'C2' ]

new_data = {}
luminosity = 3000.0 # used for S/sqrt(B) only.

for type in file_types:
    new_data[type] = []

for k, recursive_folder in enumerate(recursive_folders):

    for j, file_type in enumerate(file_types):

        filename = 'histo_CF_' + file_type + '.dat'
        filePath = os.path.join(folder,recursive_folder,filename)
        print(filePath)

        with open(filePath) as f:

            Cutflow = f.readlines()[4:11]

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
    outfile = 'cutflow_' + file_type + '.dat'
    with open(outfile, 'w') as f:
        # generate header line
        f.write(file_type+'\t')
        for recursive_folder in recursive_folders:
            f.write(recursive_folder+'\t')
        f.write('\n')
        for h, row in enumerate(list_of_rows):
            for i, item in enumerate(row):
                if i == 0:
                    f.write(Cnum[int(item)])
                else:
                    f.write(str(item))
                    if h>0:
                        if float(list_of_rows[h-1][i]) != 0:
                            cuteff = 100*float(item)/float(list_of_rows[h-1][i])
                        else:
                            cuteff = 0
                        f.write(' ')
                        f.write('(%.2f' % cuteff)
                        f.write('%)')
                if i == len(row) - 1:
                    f.write('\n')
                else:
                    f.write('\t')


summary_outfile = 'cutflow' + '_summary.dat'
with open(summary_outfile, 'w') as f:
    for i in range(0,10):
        if i == 0 or i == 1 or i == 2:
            f.write(''+recursive_folders[i]+'\t')
        elif i == 3:
            f.write(''+'S/B'+'\t')
        elif i == 4:
            f.write(''+'S/sqrt(B)'+'\t')
        elif i == 5:
            f.write(''+'Nev_signal'+'\t')
        elif i == 6:
            f.write(''+'Nev_back'+'\t')
        elif i == 7:
            f.write(''+'Nev_back_4b'+'\t')
        elif i == 8:
            f.write(''+'S/B_4b'+'\t')
        elif i == 9:
            f.write(''+'S/sqrt(B)_4b'+'\t')
        for file_type in file_types:
            f.write(file_type+'\t')
        f.write('\n')
        # write main lines
        for j in range(len(new_data["res"])):
            for k, file_type in enumerate(file_types):
                if k == 0:
                    f.write(Cnum[j] +'\t')
                # switch for tables, first sig+bkg, then S/B and sqrt.
                if i == 0 or i==1 or i==2:
                    f.write(new_data[file_type][j][i+1]+'\t')
                elif i == 3:
                    if float(new_data[file_type][j][2]) == 0.0:
                        SoverB = 0
                    else:
                        SoverB = float(new_data[file_type][j][1])/float(new_data[file_type][j][2])
                    f.write('%e'% SoverB + '\t')
                elif i == 4:
                    if float(new_data[file_type][j][2]) == 0.0:
                        SoverRootB = 0
                    else:
                        SoverRootB = math.sqrt(luminosity)*float(new_data[file_type][j][1])/math.sqrt(float(new_data[file_type][j][2]))
                    f.write('%e'% SoverRootB + '\t')
                elif i == 5:
                    Nev = luminosity*float(new_data[file_type][j][1])
                    f.write('%e'% Nev + '\t')
                elif i == 6:
                    Nev_back = luminosity*float(new_data[file_type][j][2])
                    f.write('%e'% Nev_back + '\t')
                elif i == 7:
                    Nev_back_4b = luminosity*float(new_data[file_type][j][3])
                    f.write('%e'% Nev_back_4b + '\t')
                elif i == 8:
                    if float(new_data[file_type][j][3]) == 0.0:
                        SoverB = 0
                    else:
                        SoverB = float(new_data[file_type][j][1])/float(new_data[file_type][j][3])
                    f.write('%e'% SoverB + '\t')
                elif i == 9:
                    if float(new_data[file_type][j][3]) == 0.0:
                        SoverRootB = 0
                    else:
                        SoverRootB = math.sqrt(luminosity)*float(new_data[file_type][j][1])/math.sqrt(float(new_data[file_type][j][3]))
                    f.write('%e'% SoverRootB + '\t')
            f.write('\n')
        f.write('\n\n')


         #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
    
