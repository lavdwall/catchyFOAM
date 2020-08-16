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

T0 = 900.0   # Temperature 
P0 = 1.3e5   # Pressure
Y0 = {'CH4':0.55, 'O2':0.45} # Composition

epsB = 0.4             # Voidage
epsC = 0.0             # Catalyst porosity
velocity = 1.5         # Interstitial inlet velocity
cat_area_per_vol = 3e2 # Catalyst surface area per unit reactor(!) volume

ID = 0.005
L = 0.01

# input file containing the surface reaction mechanism
cti_file = 'constant/mechanism/methane_pox_on_pt.cti'
gas = ct.Solution(cti_file, 'gas')
surf = ct.Interface(cti_file,'Pt_surf', [gas])

# output files
output_filename_gas = 'pfr_gas.csv'
output_filename_surf = 'pfr_surf.csv'

#####################################################################

pfr = surfacePFR(gas, surf, 1.0, ID, epsB, cat_area_per_vol/(1.0-epsB), epsC = epsC, isothermal=True)
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
ax4.plot(pfrdata['z'], pfrdata['Y_H2O'], color='limegreen', label='H2O')
ax4.plot(ofdata['y'], ofdata['H2O'], '--', color='limegreen')
ax4.plot(pfrdata['z'], pfrdata['Y_H2'], color='navy', label='H2')
ax4.plot(ofdata['y'], ofdata['H2'], '--', color='navy')
ax4.plot(pfrdata['z'], pfrdata['Y_CO'], color='darkorange', label='CO')
ax4.plot(ofdata['y'], ofdata['CO2'], '--', color='darkorange')
ax4.plot(pfrdata['z'], pfrdata['Y_CO2'], color='crimson', label='CO2')
ax4.plot(ofdata['y'], ofdata['CO2'], '--', color='crimson')
ax4.legend(loc=0, frameon=False)
ax5.set_xlabel('Axial position [m]')
ax5.set_ylabel('Mass fraction [-]')
ax5.plot(pfrdataS['z'], pfrdataS['X_PTS'], color='crimson', label='PT(S)')
ax5.plot(ofdata['y'], ofdata['PTScatalyst'], '--', color='crimson')
ax5.plot(pfrdataS['z'], pfrdataS['X_OHS'], color='darkorange', label='OH(S)')
ax5.plot(ofdata['y'], ofdata['OHScatalyst'], '--', color='darkorange')
ax5.plot(pfrdataS['z'], pfrdataS['X_OS'], color='dodgerblue', label='O(S)')
ax5.plot(ofdata['y'], ofdata['OScatalyst'], '--', color='dodgerblue')
ax5.plot(pfrdataS['z'], pfrdataS['X_COS'], color='navy', label='CO(S)')
ax5.plot(ofdata['y'], ofdata['COScatalyst'], '--', color='navy')
ax5.legend(loc=0, frameon=False)
plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')
print("Plot saved to 'results_comparison.png'")

#####################################################################
