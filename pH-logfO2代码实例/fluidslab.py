from numpy import log10,abs
from scipy.optimize import root,fsolve,leastsq
import DEW_HL
import pandas as pd
from tqdm import tqdm
import numpy as np

T = 600 #C
P = 50  #kbar

def read_aqueous_data():
	aqueous_species_databasae = pd.read_excel('logk_DEW2019_HP.xlsx',sheet_name='aqueous_species')
	return aqueous_species_databasae
def read_minerals_data():
	minerals_databasae = pd.read_excel('logk_DEW2019_HP.xlsx',sheet_name='minerals')
	return minerals_databasae
def read_gamma_data():
	gamma_databasae = pd.read_excel('gamma.xlsx',sheet_name='gamma')
	return gamma_databasae
# main_concentration_list [CO3]
# all_concentration_list [CO3,HCO3,CO2,CH4]
#                         0   1    2  3  

aqueous_sheet = read_aqueous_data()
P_list = aqueous_sheet['P/kbar']
T_list = aqueous_sheet['T/C']
logk_OH_list = aqueous_sheet['OH-']
logk_HCO3_list = aqueous_sheet['HCO3-']
logk_CO2_list = aqueous_sheet['CO2']
logk_CH4_list = aqueous_sheet['CH4']
logk_CH3COO_list = aqueous_sheet['CH3COO-']
gamma_sheet = read_gamma_data()
a_gamma = gamma_sheet['a_gamma']
b_gamma = gamma_sheet['b_gamma']
# print(T_list[75])
# print(P_list[75])
for i in range(len(P_list)):
	if T == T_list[i] and P == P_list[i]:
		k = i


Z_list = [-2,-1,0,0,-1]
def fx(x,k):
	eqs = []
	eqs.append(x[0]+x[1]+x[2]+x[3]-main_concentration)  # C mole
	eqs.append((-pH)+log10(x[0])+logr[0]-log10(x[1])-logr[1]-(logk_HCO3_list[k]))  #  HCO3=H + CO3
	eqs.append((-pH)+log10(x[1])+logr[1]-log10(x[2])-logr[2]-log10(aH2O)-(logk_CO2_list[k])) # CO2+ H2O = H + HCO3
	eqs.append(log10(x[0])+logr[0]+2*log10(aH2O)+2*(-pH)-2*logfO2-log10(x[3])-logr[3]-logk_CH4_list[k]) # CH4+2O2 = 2H + 2H2O +CO3
	eqs.append(2*log10(x[0])+2*logr[0]-3*pH-2*logfO2-log10(x[4])-logr[4]-logk_CH3COO_list[k])  # CH3COO- + 2O2 = 2CO32-  + 3H+
	return eqs


def equation_solver(k):
	x0 = [0.1]*5
	solution = list(leastsq(fx,x0,k)[0])
	# print(solution)
	# print(solution[0]+solution[4]+solution[7]+solution[9])
	# print(solution[1]+solution[2]+solution[3]+solution[4]+solution[7]+solution[8])
	return solution

pH_list = np.arange(0,10,1)
logfO2_list = np.arange(-24,-12,2)
result_list = []
tip = 0
for i in range(len(pH_list)):
	for ii in range(len(logfO2_list)):
		pH = pH_list[i]
		logfO2 = logfO2_list[ii]
		main_concentration = 0.2
		logr = [1]*5
		s = [0.1]*5
		for iii in range(6):
			aH2O = 55.56/(55.56+sum(s))
			logr = DEW_HL.all_Debye_huckel(P*1000,T,s,Z_list,a_gamma[k],b_gamma[k])
			s = equation_solver(k)
		result_list.append(s)


from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()  #定义新的三维坐标轴
ax3 = fig.gca(projection='3d')

CO3_list = []
HCO3_list = []
CO2_list = []
CH4_list = []
CH3COO_list = []

for i in range(len(pH_list)*len(logfO2_list)):
	CO3_list.append(np.log10(result_list[i][0]))
	HCO3_list.append(np.log10(result_list[i][1]))
	CO2_list.append(np.log10(result_list[i][2]))
	CH4_list.append(np.log10(result_list[i][3]))
	CH3COO_list.append(np.log10(result_list[i][4]))
CO3_array=np.array(CO3_list).reshape((len(pH_list),len(logfO2_list)))
HCO3_array=np.array(HCO3_list).reshape((len(pH_list),len(logfO2_list)))
CO2_array=np.array(CO2_list).reshape((len(pH_list),len(logfO2_list)))
CH4_array=np.array(CH4_list).reshape((len(pH_list),len(logfO2_list)))
CH3COO_array=np.array(CH3COO_list).reshape((len(pH_list),len(logfO2_list)))

X, Y = np.meshgrid(logfO2_list,pH_list)
surf5 = ax3.plot_surface(X,Y,CH3COO_array,color='grey',alpha=0.8,label='CH3COO')
surf1 = ax3.plot_surface(X,Y,CO3_array,color='b',alpha=0.8,label='CO3')
surf2 = ax3.plot_surface(X,Y,HCO3_array,color='r',alpha=0.8,label='HCO3')
surf3 = ax3.plot_surface(X,Y,CO2_array,color='y',alpha=0.8,label='CO2')
surf4 = ax3.plot_surface(X,Y,CH4_array,color='g',alpha=0.8,label='CH4')

ax3.set_xlabel('logfO2')
ax3.set_ylabel('pH')
ax3.set_zlabel('concentration')
surf1._facecolors2d=surf1._facecolor3d
surf1._edgecolors2d=surf1._edgecolor3d
surf2._facecolors2d=surf2._facecolor3d
surf2._edgecolors2d=surf2._edgecolor3d
surf3._facecolors2d=surf3._facecolor3d
surf3._edgecolors2d=surf3._edgecolor3d
surf4._facecolors2d=surf4._facecolor3d
surf4._edgecolors2d=surf4._edgecolor3d
surf5._facecolors2d=surf5._facecolor3d
surf5._edgecolors2d=surf5._edgecolor3d
ax3.legend()
plt.show()



fw = open('./result.txt','w')
for i in range(len(pH_list)):
	for ii in range(len(logfO2_list)):
		fw.write('{} {} {} {} {} {} {}\n'.format(pH_list[i],logfO2_list[ii],CO3_array[i,ii],HCO3_array[i,ii],CO2_array[i,ii],CH4_array[i,ii],CH3COO_array[i,ii]))
