import matplotlib.pyplot as plt
import numpy as np
import math

# constant values
H1 = 20
beta1 = 400
rho1 = 2.59
rho2 = 3
f = np.linspace(0.1, 50, num=500)

# variations of S-wave velocity in halfspace (layer 2)
beta2 = [400, 600, 800, 1000, 1200]

# amplification spectrum
AMP = []
for b2 in beta2:
    amp = []
    for i in range (0, len(f)):
        omega = 2 * np.pi * f[i]
        k1 = omega / beta1
        S1 = k1 * H1
        alpha1 = (rho1 * beta1) / (rho2 * b2)
        x2 = math.cos(S1)
        y2 = alpha1 * math.sin(S1)
        amp.append((x2**2 + y2**2)**(-0.5))
    AMP.append(amp)

# plots
i = 0
for amp in AMP:
    plt.title("$Amplifikacijski$ $spektar$ $za$ $\\beta_{2} = %i$ $m/s$" %beta2[i])
    plt.xlabel("f [Hz]")
    plt.ylabel("AMP($\omega$)")
    plt.plot(f,np.array(amp))
    plt.show()
    i = i+1