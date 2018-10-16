import numpy as np
import csv
from scipy import optimize
from scipy.optimize import least_squares
import matplotlib.pyplot as plt


b =[[],[]]
i=0
with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',')
    for row in csv_reader:
        b[0].append(float(row[0]))
        b[1].append(int(i*30))
        i=i+1

a= np.asarray(b)


t = range(90)

print("Enter the parameters in the following order as per the equations: \na  --> Rate constant for logistic growth (Eg: 0.04) \nfr --> Fraction of bacteria that have gained resistance (Eg: 0.01) \nbr --> Fraction of bacteria that have lost resistance (Eg: 0.05)\nir --> Rate constant for phage-bacteria interaction (Eg: 0.009)\ndr --> Rate constant for Phage decay (Eg: 1)\nc  --> Maximum possible concentration of bacteria in medium (Eg: 20)\nbs --> Average number of phages relesed when bacteria is lysed (Eg: 150)\ng  --> Asymmetry factor (Eg: 0.5)")
para1 = [2.79301086e-03, 4.26925554e-01, 7.66731298e-03, 1.34153693e-01,
       5.86543261e-03, 5.83768687e-01, 2.26058997e+02, 9.99999989e-01, 900]




def model(para):
    a1 = para[0]
    fr = para[1]
    br = para[2]
    ir = para[3]
    dr = para[4]
    c = para[5]
    bs = para[6]
    g = para[7]
    t0 = para[8]
    dt = 0.1
    i0 = int(t0/dt)



    sb0 = 0.15
    v0 = 0.01


    data = [[sb0], [0], [v0]]
    realdata=[sb0]
    t=np.arange(0, 100000, dt)

    for i,j in enumerate(t):
        sb = data[0][i]
        rb = data[1][i]
        v = data[2][i]

        data[0].append(sb+dt*(a1*sb*(1-np.power((sb+rb)/c,g))-fr*sb+br*rb-ir*sb*v))
        data[1].append(rb+dt*(a1*rb*(1-(np.power((sb+rb)/c,g)))+fr*sb-br*rb))
        if i>i0:
            data[2].append(v+dt*(-ir*sb*v-dr*v+bs*ir*data[0][i-i0-1]*data[2][i-i0-1]))
        else:
            data[2].append(v + dt * (-ir * sb * v - dr * v ))
        realdata.append(data[0][i+1]+data[1][i+1])

    bin=data[0].pop()+data[1].pop()+data[2].pop()+realdata.pop()

    meas = []
    for i in range(90):
        meas.append(realdata[int(i/dt)*30])

    return np.asarray(meas)



def fun(para):

    return model(para)-a[0]


print("hi")



res = least_squares(fun, para1, bounds=([0,0,0,0,0,0,0,0,0], [10,1,1,5,5,1.5,300,1,1800]), max_nfev=50, verbose=2)

print(res)

plt.plot(a[1], a[0])

plt.show()

print("hi")