#!/usr/bin/env python

import numpy as np
import generate2d
import os
from scipy import optimize
import math

dir0 = os.getcwd()
with open(dir0 + "/dir_0.0_0/atom.config") as f:
	lines = f.readlines()
	lines = lines[2:5]
	a = []
	b = []
	c = []
        d = [0.0, 0.0, 0.0]
	for i in range(0,len(lines)):
		lines[i] = lines[i].split()[0:3]
	for j in lines[0]:
		a.append(float(j))
	for j in lines[1]:
		b.append(float(j))
	for j in lines[2]:
		c.append(float(j))
	a = np.array(a)
	b = np.array(b)
	c = np.array(c)
        d[0] = a[0]*b[1] - a[1]*b[0]
        d[1] = a[0]*b[2] - a[2]*b[0]
        d[2] = a[1]*b[2] - a[2]*b[1]
        area = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2)
#print volume


#print generate2d.enumber
#print generate2d.strains
#print generate2d.csystem
osawk = '''awk '/END/ {print $5}' RELAXSTEPS | tail -1'''
#osawk = '''awk '/E_tot\(eV\)/ {print $3}' REPORT | tail -1'''
E0 = []
for k in range(generate2d.enumber):
    locals()['energy'+str(k)] = []
    E0tmp = 0.0
    for strain in generate2d.strains:
        os.chdir(dir0 + "/dir_" + str(strain) + "_" + str(k))
        tmp = float(os.popen(osawk).read().strip())
        locals()['energy'+str(k)].append(tmp)
    #print locals()['energy'+str(generate.enumber)]
    E0.append(min(locals()['energy'+str(k)]))
#print E0
for i in range(generate2d.enumber):
    locals()['dE' + str(i)] = []
    locals()['strain' + str(i)] = []
    for energy in locals()['energy'+str(i)]:
        dEtmp = (float(energy) - float(E0[i])) / area * 16.02000
        if dEtmp >= 0 and dEtmp <= 2:
            locals()['dE' + str(i)].append(dEtmp)
            index = (locals()['energy'+str(i)]).index(energy)
            locals()['strain' + str(i)].append(generate2d.strains[index])

parameter = []
for i in range(generate2d.enumber):
    x = []
    y = []
    x = locals()['strain' + str(i)]
    y = locals()['dE' + str(i)]
    x1 = np.arange(min(x), max(x), 0.001)
    def fmax(x,a):
        return a*x**2
    fita, fitb = optimize.curve_fit(fmax,x,y)
    parameter.append(fita[0])



def elastic_oblique(parameter):

    elastic_constant=[[0 for i in range(3)]for i in range(3)]

    elastic_constant[0][0] = 2*parameter[0]
    elastic_constant[1][1] = 2*parameter[1]
    elastic_constant[2][2] = 2*parameter[2]

    elastic_constant[0][1] = parameter[3] - parameter[0] - parameter[1]
    elastic_constant[0][2] = parameter[4] - parameter[0] - parameter[2]
    elastic_constant[1][2] = parameter[5] - parameter[1] - parameter[2]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    return elastic_constant

def elastic_rectanglar(parameter):

    elastic_constant=[[0 for i in range(3)]for i in range(3)]

    elastic_constant[0][0] = 2*parameter[0]
    elastic_constant[1][1] = 2*parameter[1]
    elastic_constant[2][2] = 2*parameter[2]

    elastic_constant[0][1] = parameter[3] - parameter[0] - parameter[1]

    elastic_constant[1][0] = elastic_constant[0][1]
    return elastic_constant


def elastic_square(parameter):

    elastic_constant=[[0 for i in range(3)]for i in range(3)]

    elastic_constant[0][0] = 2*parameter[2]
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2*parameter[1]

    elastic_constant[0][1] = parameter[1] - 2*parameter[2]

    elastic_constant[1][0] = elastic_constant[0][1]
    return elastic_constant



def elastic_hexagon(parameter):

    elastic_constant=[[0 for i in range(3)]for i in range(3)]

    elastic_constant[0][0] = (parameter[0] + 4*parameter[1])/2.00
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2*parameter[1]

    elastic_constant[0][1] = (parameter[0] - 4*parameter[1])/2.00

    elastic_constant[1][0] = elastic_constant[0][1]
    return elastic_constant



if generate2d.csystem == 1:
    elastic_constant=elastic_oblique(parameter)
elif generate2d.csystem == 2:
    elastic_constant=elastic_rectanglar(parameter)
elif generate2d.csystem == 3:
    elastic_constant=elastic_square(parameter)
elif generate2d.csystem == 4:
    elastic_constant=elastic_hexagon(parameter)
else:
    print "space group error!"


Ex=(elastic_constant[0][0]*elastic_constant[1][1]-elastic_constant[1][0]*elastic_constant[0][1])/elastic_constant[1][1]
Ey=(elastic_constant[0][0]*elastic_constant[1][1]-elastic_constant[1][0]*elastic_constant[0][1])/elastic_constant[0][0]

Muxy=elastic_constant[1][0]/elastic_constant[1][1]
Muyx=elastic_constant[0][1]/elastic_constant[0][0]

compliance_constant=np.linalg.inv(elastic_constant)

with open(dir0+"/elastic.txt",'w+') as f:
    print >> f, "Elastic constants (GPa):"
    print >> f, "c%d%d=%-8.2f"%(1, 1, elastic_constant[0][0]),
    print >> f, "c%d%d=%-8.2f"%(1, 2, elastic_constant[0][1]),
    print >> f, "c%d%d=%-8.2f"%(1, 6, elastic_constant[0][2]),
    print >> f, "\n",
    print >> f, "c%d%d=%-8.2f"%(2, 1, elastic_constant[1][0]),
    print >> f, "c%d%d=%-8.2f"%(2, 2, elastic_constant[1][1]),
    print >> f, "c%d%d=%-8.2f"%(2, 6, elastic_constant[1][2]),
    print >> f, "\n",
    print >> f, "c%d%d=%-8.2f"%(1, 6, elastic_constant[2][0]),
    print >> f, "c%d%d=%-8.2f"%(2, 6, elastic_constant[2][1]),
    print >> f, "c%d%d=%-8.2f"%(6, 6, elastic_constant[2][2]),
    print >> f, "\n"
    print >> f, "Compliance constants:"
    print >> f, "s%d%d=%-8.4f"%(1, 1, compliance_constant[0][0]),
    print >> f, "s%d%d=%-8.4f"%(1, 2, compliance_constant[0][1]),
    print >> f, "s%d%d=%-8.4f"%(1, 6, compliance_constant[0][2]),
    print >> f, "\n",
    print >> f, "s%d%d=%-8.4f"%(2, 1, compliance_constant[1][0]),
    print >> f, "s%d%d=%-8.4f"%(2, 2, compliance_constant[1][1]),
    print >> f, "s%d%d=%-8.4f"%(2, 6, compliance_constant[1][2]),
    print >> f, "\n",
    print >> f, "s%d%d=%-8.4f"%(1, 6, compliance_constant[2][0]),
    print >> f, "s%d%d=%-8.4f"%(2, 6, compliance_constant[2][1]),
    print >> f, "s%d%d=%-8.4f"%(6, 6, compliance_constant[2][2]),
    print >> f, "\n"
    print >> f, "Ex         Ey         Poisson(xy)     Poisson(yx)"
    print >> f, "------     ------     ----------      ----------"
    print >> f, "%-8.2f   %-8.2f   %-8.2f        %-8.2f"%(Ex, Ey, Muxy, Muyx)

