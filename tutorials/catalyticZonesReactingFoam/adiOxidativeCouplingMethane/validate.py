import cantera as ct
import numpy as np

import os
import sys
catchy_etc = os.environ.get('CATCHY_ETC') 
sys.path.append(os.path.join(catchy_etc, "validationTools"))

from reactormodels_thermo import surfacePFR

#######################################################################
# Input Parameters
#######################################################################

T0 = 973.0   # Temperature 
P0 = 1.1e5   # Pressure
Y0 = {'CH4':0.0977562, 'O2':0.048746, 'N2':0.853498} # Composition

epsB = 0.4             # Voidage
epsC = 0.0             # Catalyst porosity
velocity = 1.5         # Interstitial inlet velocity
cat_area_per_vol = 1e6 # Catalyst surface area per unit reactor(!) volume

ID = 0.005
L = 0.01

# input file containing the surface reaction mechanism
cti_file = 'constant/mechanism/ocm_lct_srla2o3.cti'
gas = ct.Solution(cti_file, 'gas')
surf = ct.Interface(cti_file,'surface1', [gas])

# output files
output_filename_gas = 'pfr_gas.csv'
output_filename_surf = 'pfr_surf.csv'

#####################################################################

pfr = surfacePFR(gas, surf, 1.0, ID, epsB, cat_area_per_vol/(1.0-epsB), epsC = epsC, isothermal=False)
pfr.gas.TPY = T0, P0, Y0
pfr.mdot = pfr.A*pfr.gas.density*velocity
y0 = np.hstack((pfr.gas.T, pfr.gas.P, pfr.gas.Y))
states, surfstates = pfr.solve(y0, L, 1000)

states.write_csv('pfr_gas.csv', cols=('extra','T','density','P','Y'))
print("Gas phase results saved to '{0}'".format(output_filename_gas))
surfstates.write_csv('pfr_surf.csv', cols=('extra','T','X'))
print("Surface results saved to '{0}'".format(output_filename_surf))

#####################################################################

import matplotlib
import matplotlib.pyplot as plt
import glob

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['font.family'] = "serif"

latestTime = os.popen("foamListTimes | tail -1").read().split('\n')[0]
pfrdata = np.genfromtxt('pfr_gas.csv', dtype=float, delimiter = ',', names=True)
pfrdataS = np.genfromtxt('pfr_surf.csv', dtype=float, delimiter = ',', names=True)
offiles = glob.glob('postProcessing/sample/{}/*.csv'.format(latestTime))
ofdata = np.genfromtxt(offiles[0], dtype=float, delimiter = ',', names=True)

fig, [[ax1,ax2,ax3],[ax4,ax5,ax6]] = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
ax1.set_xlabel('Axial position [m]')
ax1.set_ylabel('CH4 mass fraction [-]')
ax1.plot(pfrdata['z'], pfrdata['Y_CH4'], color='crimson', label='Cantera')
ax1.plot(ofdata['y'], ofdata['CH4'], '--', color='crimson', label='catchyFOAM')
handles, labels = ax1.get_legend_handles_labels()
ax2.set_xlabel('Axial position [m]')
ax2.set_ylabel('O2 mass fraction [-]')
ax2.plot(pfrdata['z'], pfrdata['Y_O2'], color='crimson', label='Cantera')
ax2.plot(ofdata['y'], ofdata['O2'], '--', color='crimson', label='catchyFOAM')
ax3.set_xlabel('Axial position [m]')
ax3.set_ylabel('Temperature [K]')
ax3.plot(pfrdata['z'], pfrdata['T'], color='crimson', label='Cantera')
ax3.plot(ofdata['y'], ofdata['T'], '--', color='crimson', label='catchyFOAM')
ax6.set_xlabel('Axial position [m]')
ax6.set_ylabel('Density [kg/m3]')
ax6.plot(pfrdata['z'], pfrdata['density'], color='crimson', label='Cantera')
ax6.plot(ofdata['y'], ofdata['thermorho'], '--', color='crimson', label='catchyFOAM')
ax4.set_xlabel('Axial position [m]')
ax4.set_ylabel('Mass fraction [-]')
ax4.plot(pfrdata['z'], pfrdata['Y_C2H2'], color='limegreen', label='C2H2')
ax4.plot(ofdata['y'], ofdata['C2H2'], '--', color='limegreen')
ax4.plot(pfrdata['z'], pfrdata['Y_C2H4'], color='navy', label='C2H4')
ax4.plot(ofdata['y'], ofdata['C2H4'], '--', color='navy')
ax4.plot(pfrdata['z'], pfrdata['Y_C2H6'], color='dodgerblue', label='C2H6')
ax4.plot(ofdata['y'], ofdata['C2H6'], '--', color='dodgerblue')
ax4.plot(pfrdata['z'], pfrdata['Y_CO'], color='darkorange', label='CO')
ax4.plot(ofdata['y'], ofdata['CO2'], '--', color='darkorange')
ax4.plot(pfrdata['z'], pfrdata['Y_CO2'], color='crimson', label='CO2')
ax4.plot(ofdata['y'], ofdata['CO2'], '--', color='crimson')
ax4.legend(loc=0, frameon=False)
ax5.set_xlabel('Axial position [m]')
ax5.set_ylabel('Mass fraction [-]')
ax5.plot(pfrdataS['z'], pfrdataS['X__S_'], color='crimson', label='S_')
ax5.plot(ofdata['y'], ofdata['_S_catalyst'], '--', color='crimson')
ax5.plot(pfrdataS['z'], pfrdataS['X_O_S'], color='darkorange', label='O_S')
ax5.plot(ofdata['y'], ofdata['O_Scatalyst'], '--', color='darkorange')
ax5.plot(pfrdataS['z'], pfrdataS['X_OH_S'], color='dodgerblue', label='OH_S')
ax5.plot(ofdata['y'], ofdata['OH_Scatalyst'], '--', color='dodgerblue')
ax5.plot(pfrdataS['z'], pfrdataS['X_CH3O_S'], color='navy', label='CH3O_S')
ax5.plot(ofdata['y'], ofdata['CH3O_Scatalyst'], '--', color='navy')
ax5.plot(pfrdataS['z'], pfrdataS['X_CH2O_S'], color='limegreen', label='CH2O_S')
ax5.plot(ofdata['y'], ofdata['CH2O_Scatalyst'], '--', color='limegreen')
ax5.legend(loc=0, frameon=False)
plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')
print("Plot saved to 'results_comparison.png'")

#####################################################################
