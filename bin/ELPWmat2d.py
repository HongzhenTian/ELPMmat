#!/usr/bin/env python

#Program:
#This program is toused to calculate elastic constants, in-plane Young`s,\
#shear and in-plane Poisson`s ratio
#History:
#20180620 HongzhenTian First release

import numpy as np
import generate2d
import os
import subprocess

dir0 = os.getcwd()
os.system("generate2d.py")
strains = [-0.018, -0.015, -0.012, -0.009, -0.006, -0.003, 0.000,\
           0.003, 0.006, 0.009, 0.012, 0.015, 0.018]

def grep(filename,arg):
    process = subprocess.Popen(['grep','-i',arg,filename],stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout

pspnumbers = len(np.unique(generate2d.numbers))
pspfile = []
for i in range(1,pspnumbers+1):
    psp = grep("etot.input", "in.psp" + str(i))
    pspfile.append(psp.split("=")[1].strip())

dirname = ''
if generate2d.csystem == 1:
    dirname = 'oblique'
elif generate2d.csystem == 2:
    dirname = 'rectanglar'
elif generate2d.csystem == 3:
    dirname = 'square'
elif generate2d.csystem == 4:
    dirname = 'hexagon'
for strain in strains:
    for j in range(generate2d.enumber):
        os.system("mkdir dir" + "_" + str(strain) + "_" + str(j) )
        os.system("cp " + dirname + "/atom.config" + "_" + str(strain) + "_" + str(j) + \
                  " dir" + "_" + str(strain) + "_" + str(j) +"/atom.config")
        os.system("cp etot.input " + " dir" + "_" + str(strain) + "_" + str(j))
        for psp in pspfile:
            os.system("cp " + psp + " dir" + "_" + str(strain) + "_" + str(j))
        os.system("cp job.pbs" + " dir" + "_" + str(strain) + "_" + str(j))
        os.chdir(dir0 + "/dir" + "_" + str(strain) + "_" + str(j))
        os.system("qsub job.pbs")
        os.system("sleep 10")
        os.chdir(dir0)
os.system("rm -rf " + dirname)


