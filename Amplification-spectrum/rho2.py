import matplotlib.pyplot as plt
import numpy as np
import math

# constant values
H1 = 20
beta1 = 400
beta2 = 800
rho1 = 2.59
f = np.linspace(0.1, 50, num=500)

# variations of density in layer 2
rho2 = [2.1, 2.6, 3, 3.6, 4.1]

# amplification spectrum
AMP = []
for r2 in rho2:
    amp = []
    for i in range (0, len(f)):
        omega = 2 * np.pi * f[i]
        k1 = omega / beta1
        S1 = k1 * H1
        alpha1 = (rho1 * beta1) / (r2 * beta2)
        x2 = math.cos(S1)
        y2 = alpha1 * math.sin(S1)
        amp.append((x2**2 + y2**2)**(-0.5))
    AMP.append(amp)

# plots
i = 0
for amp in AMP:
    plt.title("$Amplifikacijski$ $spektar$ $za$ $\\rho_{2} = %.2f$ $g/cm^{3}$" %rho2[i])
    plt.xlabel("f [Hz]")
    plt.ylabel("AMP($\omega$)")
    plt.plot(f,np.array(amp))
    plt.show()
    i = i+1