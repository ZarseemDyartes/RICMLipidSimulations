# IV: Head-Tail Simulation, Base Adlayer
import numpy as np
import matplotlib.pyplot as plt

# IRM Light Wavelength
lambda1 = 450e-9
lambda2 = 532e-9
lambda3 = 610e-9

# Refractive Indices
n0 = 1.52  # Substrate refractive index
n1 = 1.49  # Phosphocholine lipid head group refractive index
n2 = 1.44  # Lipid tail group refractive index (octadecane used in model)
n3 = 1.33  # PBS top solution refractive index

# Distances and Wavelength
m = np.arange(31)
npoint = len(m)
d1 = 0.5e-9  # Phosphocholine head thickness
d2 = 4e-9  # Tail thickness
d3 = 1e-9    # Adlayer thickness

# Fresnel Reflectance Coefficients
r03 = (n0 - n3) / (n0 + n3)  # Substrate to adlayer
r12 = (n1 - n2) / (n1 + n2)  # Head to tail
r21 = -r12  # Tail to head
r13 = (n1 - n3) / (n1 + n3)  # Head to top solution
r31 = -r13  # Adlayer to head

# Fresnel Transmittance Coefficients
t03 = 2 * n0 / (n0 + n3)  # Substrate to adlayer
t12 = 2 * n1 / (n1 + n2)  # Head to tail
t21 = 2 * n2 / (n1 + n2)  # Tail to head
t13 = 2 * n1 / (n1 + n3)  # Head to top solution
t31 = 2 * n3 / (n1 + n3)  # Adlayer to head

# Transfer Matrices between layers
M03 = np.array([[1, r03], [r03, 1]]) / t03  # Substrate to adlayer
M12 = np.array([[1, r12], [r12, 1]]) / t12  # Head to tail
M21 = np.array([[1, r21], [r21, 1]]) / t21  # Tail to head
M13 = np.array([[1, r13], [r13, 1]]) / t13  # Head to top solution
M31 = np.array([[1, r31], [r31, 1]]) / t31  # Adlayer to head

# Transfer Matrices within layers - 450 nm
phi1_l1 = 2 * np.pi * n1 * d1 / lambda1  # Head
phi2_l1 = 2 * np.pi * n2 * d2 / lambda1  # Tail
phi3_l1 = 2 * np.pi * n3 * d3 / lambda1  # Adlayer
M1_l1 = np.array([[np.exp(1j * phi1_l1), 0], [0, np.exp(-1 * 1j * phi1_l1)]])  # Head
M2_l1 = np.array([[np.exp(1j * phi2_l1), 0], [0, np.exp(-1 * 1j * phi2_l1)]])  # Tail
M3_l1 = np.array([[np.exp(1j * phi3_l1), 0], [0, np.exp(-1 * 1j * phi3_l1)]])  # Adlayer

# Transfer Matrices within layers - 532 nm
phi1_l2 = 2 * np.pi * n1 * d1 / lambda2  # Head
phi2_l2 = 2 * np.pi * n2 * d2 / lambda2  # Tail
phi3_l2 = 2 * np.pi * n3 * d3 / lambda2  # Adlayer
M1_l2 = np.array([[np.exp(1j * phi1_l2), 0], [0, np.exp(-1 * 1j * phi1_l2)]])  # Head
M2_l2 = np.array([[np.exp(1j * phi2_l2), 0], [0, np.exp(-1 * 1j * phi2_l2)]])  # Tail
M3_l2 = np.array([[np.exp(1j * phi3_l2), 0], [0, np.exp(-1 * 1j * phi3_l2)]])  # Adlayer

# Transfer Matrices within layers - 610 nm
phi1_l3 = 2 * np.pi * n1 * d1 / lambda3  # Head
phi2_l3 = 2 * np.pi * n2 * d2 / lambda3  # Tail
phi3_l3 = 2 * np.pi * n3 * d3 / lambda3  # Adlayer
M1_l3 = np.array([[np.exp(1j * phi1_l3), 0], [0, np.exp(-1 * 1j * phi1_l3)]])  # Head
M2_l3 = np.array([[np.exp(1j * phi2_l3), 0], [0, np.exp(-1 * 1j * phi2_l3)]])  # Tail
M3_l3 = np.array([[np.exp(1j * phi3_l3), 0], [0, np.exp(-1 * 1j * phi3_l3)]])  # Adlayer

# Light Pathway: m03 - m3 - m31 - j(m1 - m12 - m2 - m21 - m1) - m13
r_l1 = np.zeros(npoint, dtype=complex)
t_l1 = np.zeros(npoint, dtype=complex)
r_l2 = np.zeros(npoint, dtype=complex)
t_l2 = np.zeros(npoint, dtype=complex)
r_l3 = np.zeros(npoint, dtype=complex)
t_l3 = np.zeros(npoint, dtype=complex)

for j in range(1,npoint):
    M_l1 = M03 @ M3_l1 @ M31 @ np.linalg.matrix_power(M1_l1 @ M12 @ M2_l1 @ M21 @ M1_l1, j) @ M13
    r_l1[j] = M_l1[1, 0] / M_l1[0, 0]
    t_l1[j] = 1 / M_l1[0, 0]

    M_l2 = M03 @ M3_l2 @ M31 @ np.linalg.matrix_power(M1_l2 @ M12 @ M2_l2 @ M21 @ M1_l2, j) @ M13
    r_l2[j] = M_l2[1, 0] / M_l2[0, 0]
    t_l2[j] = 1 / M_l2[0, 0]

    M_l3 = M03 @ M3_l3 @ M31 @ np.linalg.matrix_power(M1_l3 @ M12 @ M2_l3 @ M21 @ M1_l3, j) @ M13
    r_l3[j] = M_l3[1, 0] / M_l3[0, 0]
    t_l3[j] = 1 / M_l3[0, 0]
    
R_l1 = np.abs(r_l1) ** 2
T_l1 = np.abs(t_l1) ** 2 * n3 / n0
R_l2 = np.abs(r_l2) ** 2
T_l2 = np.abs(t_l2) ** 2 * n3 / n0
R_l3 = np.abs(r_l3) ** 2
T_l3 = np.abs(t_l3) ** 2 * n3 / n0

Rbg = ((n0 - n3) / (n0 + n3)) ** 2
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
textct = f'n_Sol={n3:.3f}\nn_Head={n1:.3f}\nn_Tail={n2:.3f}\nn_Glass={n0:.3f}'
plt.text(0, 40, textct)
plt.legend()
plt.savefig(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailBaseAdlayer.png')
np.savetxt(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailBaseAdlayer.dat', np.column_stack((m, CIRM_l1)), delimiter='\t')

I_glass = 4000
Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l1[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l1[3]
plt.imsave(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailBaseAdlayer.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l2, 'g--o')
plt.savefig(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailBaseAdlayer.png')
np.savetxt(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailBaseAdlayer.dat', np.column_stack((m, CIRM_l2)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l2[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l2[3]
plt.imsave(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailBaseAdlayer.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l3, 'r--o')
plt.savefig(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailBaseAdlayer.png')
np.savetxt(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailBaseAdlayer.dat', np.column_stack((m, CIRM_l3)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l3[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l3[3]
plt.imsave(f'nSol={n3:.2f},nHead={n1:.2f},nTail={n2:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailBaseAdlayer.tif', Mimg, format='tiff')
