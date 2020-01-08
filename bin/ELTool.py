#!/usr/bin/env python

import numpy as np
import generate
import os
from scipy import optimize

dir0 = os.getcwd()
with open(dir0 + "/dir_0.0_0/atom.config") as f:
	lines = f.readlines()
	lines = lines[2:5]
	a = []
	b = []
	c = []
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
	s = np.cross(a,b)
	volume = np.dot(s,c)
#print volume


#print generate.enumber
#print generate.strains
#print generate.csystem
osawk = '''awk '/END/ {print $5}' RELAXSTEPS | tail -1'''
#osawk = '''awk '/E_tot\(eV\)/ {print $3}' REPORT | tail -1'''
E0 = []
for k in range(generate.enumber):
    locals()['energy'+str(k)] = []
    E0tmp = 0.0
    for strain in generate.strains:
        os.chdir(dir0 + "/dir_" + str(strain) + "_" + str(k))
        tmp = float(os.popen(osawk).read().strip())
        locals()['energy'+str(k)].append(tmp)
    #print locals()['energy'+str(generate.enumber)]
    E0.append(min(locals()['energy'+str(k)]))
#print E0
for i in range(generate.enumber):
    locals()['dE' + str(i)] = []
    locals()['strain' + str(i)] = []
    for energy in locals()['energy'+str(i)]:
        dEtmp = (float(energy) - float(E0[i])) / volume * 160.2000
        if dEtmp >= 0 and dEtmp <= 2:
            locals()['dE' + str(i)].append(dEtmp)
            index = (locals()['energy'+str(i)]).index(energy)
            locals()['strain' + str(i)].append(generate.strains[index])

parameter = []
for i in range(generate.enumber):
    x = []
    y = []
    x = locals()['strain' + str(i)]
    y = locals()['dE' + str(i)]
    x1 = np.arange(min(x), max(x), 0.001)
    def fmax(x,a):
        return a*x**2
    fita, fitb = optimize.curve_fit(fmax,x,y)
    parameter.append(fita[0])



def elastic_triclinic(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = 2*parameter[0]
    elastic_constant[1][1] = 2*parameter[1]
    elastic_constant[2][2] = 2*parameter[2]
    elastic_constant[3][3] = 2*parameter[3]
    elastic_constant[4][4] = 2*parameter[4]
    elastic_constant[5][5] = 2*parameter[5]

    elastic_constant[0][1] = parameter[6] - parameter[0] - parameter[1]
    elastic_constant[0][2] = parameter[7] - parameter[0] - parameter[2]
    elastic_constant[0][3] = parameter[8] - parameter[0] - parameter[3]
    elastic_constant[0][4] = parameter[9] - parameter[0] - parameter[4]
    elastic_constant[0][5] = parameter[10] - parameter[0] - parameter[5]
    elastic_constant[1][2] = parameter[11] - parameter[1] - parameter[2]
    elastic_constant[1][3] = parameter[12] - parameter[1] - parameter[3]
    elastic_constant[1][4] = parameter[13] - parameter[1] - parameter[4]
    elastic_constant[1][5] = parameter[14] - parameter[1] - parameter[5]
    elastic_constant[2][3] = parameter[15] - parameter[2] - parameter[3]
    elastic_constant[2][4] = parameter[16] - parameter[2] - parameter[4]
    elastic_constant[2][5] = parameter[17] - parameter[2] - parameter[5]
    elastic_constant[3][4] = parameter[18] - parameter[3] - parameter[4]
    elastic_constant[3][5] = parameter[19] - parameter[3] - parameter[5]
    elastic_constant[4][5] = parameter[20] - parameter[4] - parameter[5]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][0] = elastic_constant[0][2]
    elastic_constant[3][0] = elastic_constant[0][3]
    elastic_constant[4][0] = elastic_constant[0][4]
    elastic_constant[5][0] = elastic_constant[0][5]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[3][1] = elastic_constant[1][3]
    elastic_constant[4][1] = elastic_constant[1][4]
    elastic_constant[5][1] = elastic_constant[1][5]
    elastic_constant[3][2] = elastic_constant[2][3]
    elastic_constant[4][2] = elastic_constant[2][4]
    elastic_constant[5][2] = elastic_constant[2][5]
    elastic_constant[4][3] = elastic_constant[3][4]
    elastic_constant[5][3] = elastic_constant[0][2]
    elastic_constant[5][4] = elastic_constant[4][5]
    return elastic_constant

def elastic_monoclinic(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = 2*parameter[0]
    elastic_constant[1][1] = 2*parameter[1]
    elastic_constant[2][2] = 2*parameter[2]
    elastic_constant[3][3] = 2*parameter[3]
    elastic_constant[4][4] = 2*parameter[4]
    elastic_constant[5][5] = 2*parameter[5]

    elastic_constant[0][1] = parameter[6] - parameter[0] - parameter[1]
    elastic_constant[1][2] = parameter[7] - parameter[1] - parameter[2]
    elastic_constant[0][2] = parameter[8] - parameter[0] - parameter[2]
    elastic_constant[0][4] = parameter[9] - parameter[0] - parameter[4]
    elastic_constant[1][4] = parameter[10] - parameter[1] - parameter[4]
    elastic_constant[2][4] = parameter[11] - parameter[2] - parameter[4]
    elastic_constant[3][5] = parameter[12] - parameter[3] - parameter[5]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    elastic_constant[4][0] = elastic_constant[0][4]
    elastic_constant[4][1] = elastic_constant[1][4]
    elastic_constant[4][2] = elastic_constant[2][4]
    elastic_constant[5][3] = elastic_constant[3][5]
    return elastic_constant


def elastic_orthorhombic(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = 2*parameter[0]
    elastic_constant[1][1] = 2*parameter[1]
    elastic_constant[2][2] = 2*parameter[2]
    elastic_constant[3][3] = 2*parameter[3]
    elastic_constant[4][4] = 2*parameter[4]
    elastic_constant[5][5] = 2*parameter[5]

    elastic_constant[0][1] = parameter[6] - parameter[0] - parameter[1]
    elastic_constant[1][2] = parameter[7] - parameter[1] - parameter[2]
    elastic_constant[0][2] = parameter[8] - parameter[0] - parameter[2]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    return elastic_constant



def elastic_tetragonalII(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = parameter[0] - (parameter(4) - 2.00*parameter(5) + parameter(2))
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2*parameter[2]
    elastic_constant[3][3] = parameter[3]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = 2*parameter[1]

    elastic_constant[0][1] = parameter[4] - 2*parameter[5] + parameter[2]
    elastic_constant[1][2] = parameter[4] - parameter[2] - parameter[0]/2.00
    elastic_constant[0][5] = parameter[6] - elastic_constant[0][0]/2.00 - elastic_constant[5][5]/2.00
    elastic_constant[1][2] = elastic_constant[0][2]
    elastic_constant[1][5] = elastic_constant[0][5]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    elastic_constant[5][0] = elastic_constant[0][5]
    elastic_constant[5][1] = elastic_constant[1][5]
    return elastic_constant



def elastic_tetragonalI(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = parameter[0] - (parameter[4] - 2.00*parameter[5] + parameter[2])
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = parameter[2]*2.00
    elastic_constant[3][3] = parameter[3]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = parameter[1]*2.00

    elastic_constant[0][1] = parameter[4] - 2.00*parameter[5] + parameter[2]
    elastic_constant[0][2] = (parameter[4] -parameter[2] -parameter[0])/2.00
    elastic_constant[1][2] = elastic_constant[0][2]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    return elastic_constant




def elastic_rhombohedralII(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = (parameter[0]  + 4.00*parameter[1])/2.00
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2.00*parameter[2]
    elastic_constant[3][3] = parameter[3]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = 2.00*parameter[1]

    elastic_constant[0][1] = (parameter[0] - 4*parameter[1])/2.00
    elastic_constant[0][2] = (parameter[4] - parameter[0] - parameter[2])/2.00
    elastic_constant[1][2] = elastic_constant[0][2]
    elastic_constant[0][3] = parameter[5] - parameter[3]/2.00 - parameter[1]
    elastic_constant[0][4] = -(parameter[6] - elastic_constant[0][0]/2.00 - elastic_constant[3][3]/2.00)
    elastic_constant[3][5] = -elastic_constant[0][4]
    elastic_constant[1][3] = -elastic_constant[0][3]
    elastic_constant[1][4] = -elastic_constant[0][4]
    elastic_constant[4][5] = elastic_constant[0][3]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    elastic_constant[3][0] = elastic_constant[0][3]
    elastic_constant[4][0] = elastic_constant[0][4]
    elastic_constant[3][2] = elastic_constant[2][3]
    elastic_constant[4][1] = elastic_constant[1][4]
    elastic_constant[5][4] = elastic_constant[4][5]
    elastic_constant[5][3] = elastic_constant[3][5]
    return elastic_constant




def elastic_rhombohedralI(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = (parameter[0]  + 4.00*parameter[1])/2.00
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2.00*parameter[2]
    elastic_constant[3][3] = parameter[3]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = 2.00*parameter[1]

    elastic_constant[0][1] = (parameter[0] - 4*parameter[1])/2.00
    elastic_constant[0][2] = (parameter[4] - parameter[0] - parameter[2])/2.00
    elastic_constant[0][3] = parameter[5] - parameter[3]/2.00 - parameter[1]
    elastic_constant[1][2] = elastic_constant[0][2]
    elastic_constant[1][3] = -elastic_constant[0][3]
    elastic_constant[4][5] = elastic_constant[0][3]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    elastic_constant[3][0] = elastic_constant[0][3]
    elastic_constant[3][1] = elastic_constant[1][3]
    elastic_constant[5][4] = elastic_constant[4][5]
    return elastic_constant


def elastic_hexagonal(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = (parameter[0]  + 4.00*parameter[1])/2.00
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = 2.00*parameter[2]
    elastic_constant[3][3] = parameter[3]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = 2.00*parameter[1]

    elastic_constant[0][1] = (parameter[0] - 4*parameter[1])/2.00
    elastic_constant[0][2] = (parameter[4] - parameter[0] - parameter[2])/2.00
    elastic_constant[1][2] = elastic_constant[0][2]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]

    return elastic_constant

def elastic_cubic(parameter):

    elastic_constant=[[0 for i in range(6)]for i in range(6)]

    elastic_constant[0][0] = 2.00*parameter[1] - 2.00/3.00*parameter[2]
    elastic_constant[1][1] = elastic_constant[0][0]
    elastic_constant[2][2] = elastic_constant[0][0]
    elastic_constant[3][3] = 2.00/3.00*parameter[0]
    elastic_constant[4][4] = elastic_constant[3][3]
    elastic_constant[5][5] = elastic_constant[3][3]

    elastic_constant[0][1] = 2.00/3.00*parameter[2] - parameter[1]
    elastic_constant[1][2] = elastic_constant[0][1]
    elastic_constant[0][2] = elastic_constant[0][1]

    elastic_constant[1][0] = elastic_constant[0][1]
    elastic_constant[2][1] = elastic_constant[1][2]
    elastic_constant[2][0] = elastic_constant[0][2]
    return elastic_constant



if generate.csystem == 1:
    elastic_constant=elastic_triclinic(parameter)
elif generate.csystem == 2:
    elastic_constant=elastic_monoclinic(parameter)
elif generate.csystem == 3:
    elastic_constant=elastic_orthorhombic(parameter)
elif generate.csystem == 4:
    elastic_constant=elastic_tetragonalII(parameter)
elif generate.csystem == 5:
    elastic_constant=elastic_tetragonalI
elif generate.csystem == 6:
    elastic_constant=elastic_rhombohedralII(parameter)
elif generate.csystem == 7:
    elastic_constant=elastic_rhombohedralI
elif generate.csystem == 8:
    elastic_constant=elastic_hexagonal(parameter)
elif generate.csystem == 9:
    elastic_constant=elastic_cubic(parameter)


#print elastic_constant
#Young`s moduli: E
#shear moduli: G
#bulk moduli: K
#Poisson ratio: Mu

#E=9KG/(3K+G)
#Mu=(3K-2G)/(2(3K+G))

#Voigt
#9Kv=(c11+c22+c33)+2(c12+c23+c13)
#15Gv=(c11+c22+c33)-(c12+c23+c13)+3(c44+c55+c66)

Pv= elastic_constant[0][0] + elastic_constant[1][1] + elastic_constant[2][2]
Qv= elastic_constant[0][1] + elastic_constant[0][2] + elastic_constant[1][2]
Rv= elastic_constant[3][3] + elastic_constant[4][4] + elastic_constant[5][5]

Ev = ((Pv+2.00*Qv)*(Pv-Qv+3.00*Rv))/(3.00*(2.00*Pv+3.00*Qv+Rv))
Gv = (Pv-Qv+3.00*Rv)/15.00
Kv = Ev*Gv/(3.00*(3.00*Gv-Ev))
Muv = (Ev/(2.00*Gv))-1.00


#Reuss
#1/Kr=(s11+s22+s33)+2(s12+s23+s13)
#15/Gr=4(s11+s22+s33)-4(s12+s23+s13)+3(s44+s55+s66)

compliance_constant=np.linalg.inv(elastic_constant)
Pr = compliance_constant[0][0] + compliance_constant[1][1] + compliance_constant[2][2]
Qr = compliance_constant[0][1] + compliance_constant[0][2] + compliance_constant[1][2]
Rr = compliance_constant[3][3] + compliance_constant[4][4] + compliance_constant[5][5]

Er = 15.00/(3.00*Pr+2.00*Qr+Rr)
Gr = 15.00/(4.00*(Pr-Qr)+3.00*Rr)
Kr=Er*Gr/(3.00*(3.00*Gr-Er))
Mur=(Er/(2.00*Gr))-1.00


#Hill

Eh=(Ev+Er)/2.00
Gh=(Gv+Gr)/2.00
Kh=(Kv+Kr)/2.00
Muh=(Muv+Mur)/2.00






with open(dir0+"/elastic.txt",'w+') as f:
    print >> f, "Elastic constants (GPa):"
    print >> f, "-----------------------"
    for i in range(6):
        for j in range(6):
            if i > j:
                print >> f, "c%d%d=%-8.2f"%(j+1, i+1, elastic_constant[i][j]),
            else:
                print >> f, "c%d%d=%-8.2f"%(i+1, j+1, elastic_constant[i][j]),
        print >> f, "\n",
    print >> f, "\n"
    print >> f, "Compliance constants:"
    print >> f, "-----------------------"
    for i in range(6):
        for j in range(6):
            if i > j:
                print >> f, "s%d%d=%-8.4f"%(j+1, i+1, compliance_constant[i][j]),
            else:
                print >> f, "s%d%d=%-8.4f"%(i+1, j+1, compliance_constant[i][j]),
        print >> f, "\n",
    print >> f, "\n"
    print >> f, "Modulus           Voigt       Reuss       Hill"
    print >> f, "------------      ------      ------      ------"
    print >> f, "Young`s (GPa)     %-6.2f      %-6.2f      %-6.2f"%(Ev, Er, Eh)
    print >> f, "Shear (GPa)       %-6.2f      %-6.2f      %-6.2f"%(Gv, Gr, Gh)
    print >> f, "Bulk (GPa)        %-6.2f      %-6.2f      %-6.2f"%(Kv, Kr, Kh)
    print >> f, "Poisson           %-6.2f      %-6.2f      %-6.2f"%(Muv, Mur, Muh)
