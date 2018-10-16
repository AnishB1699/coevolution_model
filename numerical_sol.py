import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import least_squares



print("Enter the parameters in the following order as per the equations: \na  --> Rate constant for logistic growth (Eg: 0.04) \nfr --> Fraction of bacteria that have gained resistance (Eg: 0.01) \nbr --> Fraction of bacteria that have lost resistance (Eg: 0.05)\nir --> Rate constant for phage-bacteria interaction (Eg: 0.009)\ndr --> Rate constant for Phage decay (Eg: 1)\nc  --> Maximum possible concentration of bacteria in medium (Eg: 20)\nbs --> Average number of phages relesed when bacteria is lysed (Eg: 150)\ng  --> Asymmetry factor (Eg: 0.5)")
para = [2.79292105e-03, 4.26925653e-01, 7.66729446e-03, 1.34160759e-01,
       5.86542753e-03, 5.83768743e-01, 2.26058998e+02, 9.00000000e+02]


a = para[0]
fr = para[1]
br = para[2]
ir = para[3]
dr = para[4]
c = para[5]
bs = para[6]
t0 = para[7]
dt = 0.1
i0 = int(t0/dt)



sb0 = 0.15
v0 = 0.001


data = [[sb0], [0], [v0]]
realdata=[sb0]
t=np.arange(0, 3000, dt)

for i,j in enumerate(t):
    sb = data[0][i]
    rb = data[1][i]
    v = data[2][i]

    data[0].append(sb+dt*(a*sb*(1-np.power((sb+rb)/c,1))-fr*sb+br*rb-ir*sb*v))
    data[1].append(rb+dt*(a*rb*(1-(np.power((sb+rb)/c,1)))+fr*sb-br*rb))
    if i>i0:
        data[2].append(v+dt*(-ir*sb*v-dr*v+bs*ir*data[0][i-i0-1]*data[2][i-i0-1]))
    else:
        data[2].append(v + dt * (-ir * sb * v - dr * v ))
    realdata.append(data[0][i+1]+data[1][i+1])

bin=data[0].pop()+data[1].pop()+data[2].pop()+realdata.pop()

csvfile="./out.csv"

with open(csvfile, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(data)

with open("file.txt", "w") as output:
    output.write(str(data[0]))


b =[[],[]]
i=0
with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',')
    for row in csv_reader:
        b[0].append(float(row[0]))
        b[1].append(int(i*30))
        i=i+1

a= np.asarray(b)

plt.plot(a[1], a[0])
plt.plot(t, realdata)
plt.show()



