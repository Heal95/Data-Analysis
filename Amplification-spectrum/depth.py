import matplotlib.pyplot as plt
import numpy as np
import math
 
# constant values
beta1 = 400
beta2 = 800
rho1 = 2.59
rho2 = 3
alpha1 = (rho1 * beta1) / (rho2 * beta2)
f = np.linspace(0.1, 50, num=500)

# variations of layer 1 thickness
H1 = [5, 10, 20, 50, 100]

# amplification spectrum
AMP = []
for h1 in H1:
    amp = []
    for i in range (0, len(f)):
        omega = 2 * np.pi * f[i]
        k1 = omega / beta1
        S1 = k1 * h1
        x2 = math.cos(S1)
        y2 = alpha1 * math.sin(S1)
        amp.append((x2**2 + y2**2)**(-0.5))
    AMP.append(amp)

# plots
i = 0
for amp in AMP:
    plt.title("$Amplifikacijski$ $spektar$ $za$ $H_{1} = %i$ $m$" %H1[i])
    plt.xlabel("f [Hz]")
    plt.ylabel("AMP($\omega$)")
    plt.plot(f,np.array(amp))
    plt.show()
    i = i+1