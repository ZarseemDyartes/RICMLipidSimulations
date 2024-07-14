# V: Head-Tail Simulation, Interstitial Adlayers
import numpy as np
import matplotlib.pyplot as plt

# IRM Light Wavelength
lambda1 = 450e-9
lambda2 = 532e-9
lambda3 = 610e-9

# Refractive Indices
n0 = 1.52  # Substrate refractive index
n1 = 1.33  # Base adlayer refractive index
n2 = 1.49  # Phosphocholine lipid head group refractive index
n3 = 1.44  # Lipid tail group refractive index (octadecane used in model)
n4 = 1.33  # Interstitial adlayer refractive index
n5 = 1.33  # Top solution refractive index

# Distances and Wavelength
m = np.arange(31)
npoint = len(m)
d1 = 1e-9  # Base adlayer thickness
d2 = 0.5e-9  # Phosphocholine head thickness
d3 = 4e-9  # Tail thickness
d4 = 0.75e-9  # Interstitial adlayer thickness

# Fresnel Reflectance Coefficients
r01 = (n0-n1)/(n0+n1)  # Substrate to base adlayer
r12 = (n1-n2)/(n1+n2)  # Base adlayer to head
r23 = (n2-n3)/(n2+n3)  # Head to tail
r32 = -r23  # Tail to head
r24 = (n2-n4)/(n2+n4)  # Head to interstitial adlayer
r42 = -r24  # Interstitial adlayer to head
r25 = (n2-n5)/(n2+n5)  # Head to top solution

# Fresnel Transmittance Coefficients
t01 = 2*n0/(n0+n1)  # Substrate to base adlayer
t12 = 2*n1/(n1+n2)  # Base adlayer to head
t23 = 2*n2/(n2+n3)  # Head to tail
t32 = 2*n3/(n2+n3)  # Tail to head
t24 = 2*n2/(n2+n4)  # Head to interstitial adlayer
t42 = 2*n4/(n2+n4)  # Interstitial adlayer to head
t25 = 2*n2/(n2+n5)  # Head to top solution

# Transfer Matrices between layers
M01 = np.array([[1, r01], [r01, 1]]) / t01  # Substrate to base adlayer
M12 = np.array([[1, r12], [r12, 1]]) / t12  # Base adlayer to head
M23 = np.array([[1, r23], [r23, 1]]) / t23  # Head to tail
M32 = np.array([[1, r32], [r32, 1]]) / t32  # Tail to head
M24 = np.array([[1, r24], [r24, 1]]) / t24  # Head to interstitial adlayer
M42 = np.array([[1, r42], [r42, 1]]) / t42  # Interstitial adlayer to head
M25 = np.array([[1, r25], [r25, 1]]) / t25  # Head to top solution

# Transfer Matrices within layers - 450 nm
phi1_l1 = 2 * np.pi * n1 * d1 / lambda1 # Base adlayer
phi2_l1 = 2 * np.pi * n2 * d2 / lambda1 # Head
phi3_l1 = 2 * np.pi * n3 * d3 / lambda1 # Tail
phi4_l1 = 2 * np.pi * n4 * d4 / lambda1 # Interstitial adlayer
M1_l1 = np.array([[np.exp(1j * phi1_l1), 0], [0, np.exp(-1 * 1j * phi1_l1)]])  # Base adlayer
M2_l1 = np.array([[np.exp(1j * phi2_l1), 0], [0, np.exp(-1 * 1j * phi2_l1)]])  # Head
M3_l1 = np.array([[np.exp(1j * phi3_l1), 0], [0, np.exp(-1 * 1j * phi3_l1)]])  # Tail
M4_l1 = np.array([[np.exp(1j * phi4_l1), 0], [0, np.exp(-1 * 1j * phi4_l1)]])  # Interstitial Adlayer

# Transfer Matrices within layers - 532 nm
phi1_l2 = 2 * np.pi * n1 * d1 / lambda2 # Base adlayer
phi2_l2 = 2 * np.pi * n2 * d2 / lambda2 # Head
phi3_l2 = 2 * np.pi * n3 * d3 / lambda2 # Tail
phi4_l2 = 2 * np.pi * n4 * d4 / lambda2 # Interstitial adlayer
M1_l2 = np.array([[np.exp(1j * phi1_l2), 0], [0, np.exp(-1 * 1j * phi1_l2)]])  # Base adlayer
M2_l2 = np.array([[np.exp(1j * phi2_l2), 0], [0, np.exp(-1 * 1j * phi2_l2)]])  # Head
M3_l2 = np.array([[np.exp(1j * phi3_l2), 0], [0, np.exp(-1 * 1j * phi3_l2)]])  # Tail
M4_l2 = np.array([[np.exp(1j * phi4_l2), 0], [0, np.exp(-1 * 1j * phi4_l2)]])  # Interstitial Adlayer

# Transfer Matrices within layers - 610 nm
phi1_l3 = 2 * np.pi * n1 * d1 / lambda3 # Base adlayer
phi2_l3 = 2 * np.pi * n2 * d2 / lambda3 # Head
phi3_l3 = 2 * np.pi * n3 * d3 / lambda3 # Tail
phi4_l3 = 2 * np.pi * n4 * d4 / lambda3 # Interstitial adlayer
M1_l3 = np.array([[np.exp(1j * phi1_l3), 0], [0, np.exp(-1 * 1j * phi1_l3)]])  # Base adlayer
M2_l3 = np.array([[np.exp(1j * phi2_l3), 0], [0, np.exp(-1 * 1j * phi2_l3)]])  # Head
M3_l3 = np.array([[np.exp(1j * phi3_l3), 0], [0, np.exp(-1 * 1j * phi3_l3)]])  # Tail
M4_l3 = np.array([[np.exp(1j * phi4_l3), 0], [0, np.exp(-1 * 1j * phi4_l3)]])  # Interstitial Adlayer

# Light Pathway: m01 - m1 - m12 - (j-1)[m2 - m23 - m3 - m32 - m2 - m24 - m4- m42] - m2 - m23 - m3 - m32 - m2 - m25
r_l1 = np.zeros(npoint, dtype=complex)
t_l1 = np.zeros(npoint, dtype=complex)
r_l2 = np.zeros(npoint, dtype=complex)
t_l2 = np.zeros(npoint, dtype=complex)
r_l3 = np.zeros(npoint, dtype=complex)
t_l3 = np.zeros(npoint, dtype=complex)

for j in range(1,npoint):
    M_l1 = M01 @ M1_l1 @ M12 @ np.linalg.matrix_power(M2_l1 @ M23 @ M3_l1 @ M32 @ M2_l1 @ M24 @ M4_l1 @ M42, j-1) @ M2_l1 @ M23 @ M3_l1 @ M32 @ M2_l1 @ M25
    r_l1[j] = M_l1[1, 0] / M_l1[0, 0]
    t_l1[j] = 1 / M_l1[0, 0]

    M_l2 = M01 @ M1_l2 @ M12 @ np.linalg.matrix_power(M2_l2 @ M23 @ M3_l2 @ M32 @ M2_l2 @ M24 @ M4_l2 @ M42, j-1) @ M2_l2 @ M23 @ M3_l2 @ M32 @ M2_l2 @ M25
    r_l2[j] = M_l2[1, 0] / M_l2[0, 0]
    t_l2[j] = 1 / M_l2[0, 0]

    M_l3 = M01 @ M1_l3 @ M12 @ np.linalg.matrix_power(M2_l3 @ M23 @ M3_l3 @ M32 @ M2_l3 @ M24 @ M4_l3 @ M42, j-1) @ M2_l3 @ M23 @ M3_l3 @ M32 @ M2_l3 @ M25
    r_l3[j] = M_l3[1, 0] / M_l3[0, 0]
    t_l3[j] = 1 / M_l3[0, 0]
    
R_l1 = np.abs(r_l1) ** 2
T_l1 = np.abs(t_l1) ** 2 * n5 / n0
R_l2 = np.abs(r_l2) ** 2
T_l2 = np.abs(t_l2) ** 2 * n5 / n0
R_l3 = np.abs(r_l3) ** 2
T_l3 = np.abs(t_l3) ** 2 * n5 / n0

Rbg = ((n0 - n5) / (n0 + n5)) ** 2
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
textct = f'n_Sol={n5:.3f}\nn_IAd={n4:.3f}\nn_Tail={n3:.3f}\nn_Head={n2:.3f}\nn_BAd={n1:.3f}'
plt.text(10, 80, textct)
plt.legend()
plt.savefig(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailInterstitialAdlayers.png')
np.savetxt(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailInterstitialAdlayers.dat', np.column_stack((m, CIRM_l1)), delimiter='\t')

I_glass = 4000
Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l1[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l1[3]
plt.imsave(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda1*1e9:.1f}nm, HeadTailInterstitialAdlayers.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l2, 'g--o')
plt.savefig(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailInterstitialAdlayers.png')
np.savetxt(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailInterstitialAdlayers.dat', np.column_stack((m, CIRM_l2)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l2[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l2[3]
plt.imsave(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda2*1e9:.1f}nm, HeadTailInterstitialAdlayers.tif', Mimg, format='tiff')

plt.figure(2)
plt.plot(m, CIRM_l3, 'r--o')
plt.savefig(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailInterstitialAdlayers.png')
np.savetxt(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailInterstitialAdlayers.dat', np.column_stack((m, CIRM_l3)), delimiter='\t')

Mimg = I_glass * np.ones((1024, 1024), dtype=np.uint16)
Mimg[212:812, 212:812] = I_glass * CIRM_l3[2]
Mimg[412:612, 412:612] = I_glass * CIRM_l3[3]
plt.imsave(f'nSol={n5:.2f},nIAd={n4:.2f},nTail={n3:.2f},nHead={n2:.2f},nBAd={n1:.2f}, wavelength={lambda3*1e9:.1f}nm, HeadTailInterstitialAdlayers.tif', Mimg, format='tiff')
