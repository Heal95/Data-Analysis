import re
from datetime import date
import matplotlib.pyplot as plt
import numpy as np

# input earthquake hypocenter times
dat = open("CNSS_50.msh","r")
lista = []
for line in dat.readlines():
    _, year, month, day, hour, minute, sec, _, _, _, _, _, _, _, _, _, _, _, _, _ = re.split("\\s+", line)
    lista.append((date(int(year), int(month), int(day)), 3600*float(hour) + 60*float(minute) + float(sec)))

# time difference between earthquakes
raz = []
for i in range (0, len(lista)-1):
    dat1, sec1 = lista[i]
    dat2, sec2 = lista [i+1]
    a = (dat2 - dat1).days * 3600 * 24
    b = sec2 - sec1
    raz.append(a+b)

# return period
period = sum(raz) / float(len(raz))

# plot histograms
N = np.linspace(10, 100, num=3)
for n in N:
	plt.hist(raz,n)
	plt.title("Histogram vremenskih intervala - %i klasa" %n)
	plt.xlabel("dt [s]")
	plt.ylabel("ucestalosti")
	fig = plt.gcf()
	plt.show()

# print values
print ("min vremenski interval = %.2f s" %min(raz))
print ("max vremenski interval = %.2f s" %max(raz))
print ("povratni period = %f s" %period)