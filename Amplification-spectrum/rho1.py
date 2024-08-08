import matplotlib.pyplot as plt
import numpy as np
import math

# constant values
H1 = 20
beta1 = 400
beta2 = 800
rho2 = 3
f = np.linspace(0.1, 50, num=500)

# variations of density in layer 1
rho1 = [1.6, 2.1, 2.59, 3.1, 3.6]

# amplification spectrum
AMP = []
for r1 in rho1:
    amp = []
    for i in range (0, len(f)):
        omega = 2 * np.pi * f[i]
        k1 = omega / beta1
        S1 = k1 * H1
        alpha1 = (r1 * beta1) / (rho2 * beta2)
        x2 = math.cos(S1)
        y2 = alpha1 * math.sin(S1)
        amp.append((x2**2 + y2**2)**(-0.5))
    AMP.append(amp)

# plots
i = 0
for amp in AMP:
    plt.title("$Amplifikacijski$ $spektar$ $za$ $\\rho_{1} = %.2f$ $g/cm^{3}$" %rho1[i])
    plt.xlabel("f [Hz]")
    plt.ylabel("AMP($\omega$)")
    plt.plot(f,np.array(amp))
    plt.show()
    i = i+1