import matplotlib.pyplot as plt
import numpy as np
import math

# constant values
H1 = 20
beta2 = 800
rho1 = 2.59
rho2 = 3
f = np.linspace(0.1, 50, num=500)

# variations of S-wave velocity in layer 1
beta1 = [200, 300, 400, 500, 600]

# amplification spectrum
AMP = []
for b1 in beta1:
    amp = []
    for i in range (0, len(f)):
        omega = 2 * np.pi * f[i]
        k1 = omega / b1
        S1 = k1 * H1
        alpha1 = (rho1 * b1) / (rho2 * beta2)
        x2 = math.cos(S1)
        y2 = alpha1 * math.sin(S1)
        amp.append((x2**2 + y2**2)**(-0.5))
    AMP.append(amp)

# plots
i = 0
for amp in AMP:
    plt.title("$Amplifikacijski$ $spektar$ $za$ $\\beta_{1} = %i$ $m/s$" %beta1[i])
    plt.xlabel("f [Hz]")
    plt.ylabel("AMP($\omega$)")
    plt.plot(f,np.array(amp))
    plt.show()
    i = i+1