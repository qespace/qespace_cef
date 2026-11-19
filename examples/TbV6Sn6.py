import sys

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add the parent directory (project root) to sys.path
root_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(root_dir))


import crysfipy
from crysfipy.observables import spectral_weight


# Analysis of the crystal field effects in TbV6Sn6
#
# Based on:
# Phys Rev Res 7 L042045

calculateTAS = False
calculateTOF = False
calculateFormFactor = False
calculateMagnetism = True

# Bdict = {'B20':-0.1036, 'B40':-5.91e-4, 'B60':2.31e-6, 'B66':3.3e-6}
Bdict = {'B20':-0.1151, 'B40':-5.57e-4, 'B60':2.31e-6, 'B66':23e-6}
CEFpars = crysfipy.CEFpars(Bdict=Bdict, units='meV', pointGroup='D6')
cefion = crysfipy.CEFion(ion=crysfipy.Ion('Tb'), cfp=CEFpars, Hfield=[0,0,1e-3])



print(cefion.cfp)
print(cefion.__str__(evecs=True, precision=5))




e = np.linspace(0.01,2,1000)
temperature = 5


De, Dint = spectral_weight(cefion,temperature, Q=2, scheme='powder', Ei=1e6)

### Jx Jy Jz J
print(np.square(np.abs(cefion.Jx))[0,1]    )
print(np.square(np.abs(cefion.Jy))[0,1]    )
print(np.square(np.abs(cefion.Jz))[0,1]    )
print(np.einsum('ijk->jk',np.square(np.abs(cefion.J)))[0,1]    )

# print(De)
# print(Dint)
    
spectrum = np.zeros(e.shape)
for Etr, Itr in zip(De.flatten(), Dint.flatten()):
    spectrum += Itr*np.exp(-(Etr-e)**2/0.00001)


fig, ax = plt.subplots()
#label = f'T={temperature} K'
ax.plot(e, spectrum)

ax.legend()
fig.savefig(Path(__file__).parent / 'TbV6Sn6-INS-spectra.png')