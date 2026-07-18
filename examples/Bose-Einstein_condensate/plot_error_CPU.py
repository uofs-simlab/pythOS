import numpy as np
from matplotlib import pyplot as plt
import os
import sys

datafile = "CPU_vs_error_omega1.0_dx0.04_t1.0.txt"

if len(sys.argv) > 1: # file from command line
    datafile = sys.argv[1]

dx_str = '' # get dx
i = 0
while datafile[i] != 'd' or datafile[i+1] != 'x':
    i += 1
i += 2 # skip 'dx'
while datafile[i] != '_':
    dx_str += datafile[i]
    i += 1
dx = float(dx_str)

t_end_str = '' # get t_end
i += 1
if datafile[i] == 't':
    i += 1
    while datafile[i+1] != 't' and datafile[i] != '_':
        t_end_str += datafile[i]
        i += 1
    t_end = float(t_end_str)




fig, ax = plt.subplots(figsize=(8,4.8))


abbrevs = ['LT', 'Str-2-nl', 'Str-2-diffR', 'Str-3-R', 'Str-3-nl', 'Str-3-diff', 'Str-4-R', 'Str-4-nl', 'Str-4-x', 'Str-4-y', 'Ruth1', 'Ruth2', 'OS2-1', 'OS2-2', 'AKS3-1', 'AKS3-2']

leglist = ['Lie Trotter', 'Strang, nonl (2-split)', 'Strang, diff+Rabi (2-split)', 'Strang, Rabi (3-split)', 'Strang, nonl (3-split)', 'Strang, diff (3-split)', 'Strang, Rabi (4-split)', 'Strang, nonl (4-split)', 'Strang, diff x (4-split)', 'Strang, diff y (4-split)', 'Ruth-1', 'Ruth-2', 'OS2-1', 'OS2-2', 'AKS3-1', 'AKS3-2']

colors = ['tab:blue', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:orange', 'mediumslateblue', 'darkslateblue', 'tab:green', 'darkgreen', 'tab:olive', 'olivedrab','tab:red', 'darkred']

ls = ['-', '-', '-', '--', '--', '--', ':', ':', ':', ':', '-', '-', '-', '-', '-', '-']

markers = ['o', 's', 's', 's', 's', 's', 's', 's', 's', 's', '^', '^', 'v', 'v', '*', '*']


times = [[] for _ in range(len(leglist))]
err = [[] for _ in range(len(leglist))]
my_indices = [] # which plots we have data for

with open(datafile, 'r') as file:
    count = 0 

    lines = file.readlines()

    for i in range(len(lines)):
        words = lines[i].split()
        
        if len(words) == 0:
            count = -1 # reset count at linebreak
            name = ''

        if count % 2 == 0:
            name = words[0] # OS method name
            try:
                err_str = lines[i+1].split()
            except:
                err_str = ''

        try:
            index = abbrevs.index(name) # find index of method
            if index not in my_indices:
                my_indices.append(index)
        except:
            if name != '':
                print(f"Method name {name} not in list for plotting.")
            count += 1
            continue

        if count % 2 == 0: # times
            for j in range(1,len(words)):
                mytime = float(words[j])
                if len(times[index]) < j:
                    times[index].append(-1) # add el to times[index]
                try:
                    myerror = err[index][j-1]
                except:
                    myerror = -1
                if myerror-1e-08 < float(err_str[j-1]) < myerror+1e-08 or myerror < 0 or float(err_str[j-1]) < 0: # errors match, as they should
                    if mytime > 0:
                        if mytime < times[index][j-1] or times[index][j-1] < 1e-08:
                            times[index][j-1] = mytime # get minimum time
                else:
                    print(f"Error: time {times[index][j-1]} not aligned in file. Errors do not match.")
        else: # errors
            for j in range(len(words)):
                myerror = float(words[j])
                if len(err[index]) < j+1:
                    err[index].append(-1) # add el to err[index]
                if myerror > 0:
                    err[index][j] = myerror # get error

        count += 1

my_indices.sort()
plots = [0 for i in range(len(times))]
handles = []
for ind in my_indices: #range(len(times)):
    plots[ind], = plt.plot(times[ind], err[ind], color=colors[ind], linestyle=ls[ind], marker=markers[ind], label=leglist[ind]) # plot
    handles.append(plots[ind])
    #print(times[ind])

            
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel("CPU time (s)")
plt.ylabel("relative norm of error")
#plt.title(f"Total error vs. CPU time (dx = {dx}, t = {t_end})")

leg = plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5))
fig.add_artist(leg)
plt.tight_layout()

plt.show()




