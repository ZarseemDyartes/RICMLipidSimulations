# II: Bulk Simulation, Base Adlayer
import numpy as np
import matplotlib.pyplot as plt

# IRM Light Wavelength
lambda1 = 450e-9
lambda2 = 532e-9
lambda3 = 610e-9

# Refractive Indices
n0 = 1.52  # Substrate refractive index
n1 = 1.48  # Lipid bilayer refractive index - DOPC
n2 = 1.33  # PBS top solution refractive index

# Distances and Wavelength
m = np.arange(31)
npoint = len(m)
d1 = 5e-9  # Lipid bilayer thickness
d2 = 1e-9  # Adlayer thickness

# Fresnel Reflectance Coefficients
r02 = (n0 - n2) / (n0 + n2)  # Substrate to adlayer
r12 = (n1 - n2) / (n1 + n2)  # Bilayer to top solution
r21 = -r12  # Adlayer to bilayer

# Fresnel Transmittance Coefficients
t02 = 2 * n0 / (n0 + n2)  # Substrate to adlayer
t21 = 2 * n2 / (n1 + n2)  # Adlayer to bilayer
t12 = 2 * n1 / (n1 + n2)  # Bilayer to top solution

# Transfer Matrices between layers
M02 = np.array([[1, r02], [r02, 1]]) / t02  # Substrate to adlayer
M21 = np.array([[1, r21], [r21, 1]]) / t21  # Adlayer to bilayer
M12 = np.array([[1, r12], [r12, 1]]) / t12  # Bilayer to top solution

# Transfer Matrices within layers - 450 nm
phi1_l1 = 2 * np.pi * n1 * d1 / lambda1  # Bilayer
M1_l1 = np.array([[np.exp(1j * phi1_l1), 0], [0, np.exp(-1 * 1j * phi1_l1)]])  # Bilayer
phi2_l1 = 2 * np.pi * n2 * d2 / lambda1  # Adlayer
M2_l1 = np.array([[np.exp(1j * phi2_l1), 0], [0, np.exp(-1 * 1j * phi2_l1)]])  # Adlayer

# Transfer Matrices within layers - 532 nm
phi1_l2 = 2 * np.pi * n1 * d1 / lambda2  # Bilayer
M1_l2 = np.array([[np.exp(1j * phi1_l2), 0], [0, np.exp(-1 * 1j * phi1_l2)]])  # Bilayer
phi2_l2 = 2 * np.pi * n2 * d2 / lambda2  # Adlayer
M2_l2 = np.array([[np.exp(1j * phi2_l2), 0], [0, np.exp(-1 * 1j * phi2_l2)]])  # Adlayer

# Transfer Matrices within layers - 610 nm
phi1_l3 = 2 * np.pi * n1 * d1 / lambda3  # Bilayer
M1_l3 = np.array([[np.exp(1j * phi1_l3), 0], [0, np.exp(-1 * 1j * phi1_l3)]])  # Bilayer
phi2_l3 = 2 * np.pi * n2 * d2 / lambda3  # Adlayer
M2_l3 = np.array([[np.exp(1j * phi2_l3), 0], [0, np.exp(-1 * 1j * phi2_l3)]])  # Adlayer


# Light Pathway: m02 - m2 - m21 - j(m1) - m12
r_l1 = np.zeros(npoint, dtype=complex)
t_l1 = np.zeros(npoint, dtype=complex)
r_l2 = np.zeros(npoint, dtype=complex)
t_l2 = np.zeros(npoint, dtype=complex)
r_l3 = np.zeros(npoint, dtype=complex)
t_l3 = np.zeros(npoint, dtype=complex)

for j in range(1,npoint):
    M_l1 = M02 @ M2_l1 @ M21 @ (M1_l1 ** j) @ M12
    r_l1[j] = M_l1[1, 0] / M_l1[0, 0]
    t_l1[j] = 1 / M_l1[0, 0]

    M_l2 = M02 @ M2_l2 @ M21 @ (M1_l2 ** j) @ M12
    r_l2[j] = M_l2[1, 0] / M_l2[0, 0]
    t_l2[j] = 1 / M_l2[0, 0]

    M_l3 = M02 @ M2_l3 @ M21 @ (M1_l3 ** j) @ M12
    r_l3[j] = M_l3[1, 0] / M_l3[0, 0]
    t_l3[j] = 1 / M_l3[0, 0]

R_l1 = np.abs(r_l1) ** 2
T_l1 = np.abs(t_l1) ** 2 * n2 / n0
R_l2 = np.abs(r_l2) ** 2
T_l2 = np.abs(t_l2) ** 2 * n2 / n0
R_l3 = np.abs(r_l3) ** 2
T_l3 = np.abs(t_l3) ** 2 * n2 / n0

Rbg = ((n0 - n2) / (n0 + n2)) ** 2
CIRM_l1 = (R_l1 / Rbg)*100
CIRM_l2 = (R_l2 / Rbg)*100
CIRM_l3 = (R_l3 / Rbg)*100
CIRM_l1[0] = CIRM_l2[0] = CIRM_l3[0] = 100

# Creating the graph
plt.figure(2)
print("IRM450:")
print (*CIRM_l1,sep=',')
print('\n'+"IRM532:")
print (*CIRM_l2,sep=',')
print('\n'+"IRM650:")
print (*CIRM_l3,sep=",")
plt.plot(m, CIRM_l1, 'b--o', label = "IRM450")
plt.plot(m, CIRM_l2, 'g--o', label = "IRM532")
plt.plot(m, CIRM_l3, 'r--o', label = "IRM610")
plt.xlabel('Lipid Layers')
plt.ylabel('IRM I/I_0')
plt.grid()
plt.ylim([0, 110])
plt.xlim([min(m) - 0.5, max(m) + 0.5])
textct = f'n_Sol={n2:.3f}\nn_Lipid={n1:.3f}\nn_Glass={n0:.3f}'
plt.text(0, 40, textct)
plt.legend()
plt.savefig(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, BulkBaseAdlayer.png')
np.savetxt(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, BulkBaseAdlayer.dat', np.column_stack((m, CIRM_l1)), delimiter='\t')

I_glass = 4000
Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l1[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l1[3]
plt.imsave(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, BulkBaseAdlayer.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l2, 'g--o')
plt.savefig(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, BulkBaseAdlayer.png')
np.savetxt(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, BulkBaseAdlayer.dat', np.column_stack((m, CIRM_l2)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l2[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l2[3]
plt.imsave(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, BulkBaseAdlayer.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l3, 'r--o')
plt.savefig(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, BulkBaseAdlayer.png')
np.savetxt(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, BulkBaseAdlayer.dat', np.column_stack((m, CIRM_l3)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l3[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l3[3]
plt.imsave(f'nSol={n2:.2f},nLipid={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, BulkBaseAdlayer.tif', Mimg, format='tiff')
