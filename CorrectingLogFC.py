import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv('e3fp/pS6_DE_1p_2heptanone.csv')
y = np.abs(data.loc[:,["logFC"]].values)
x = (data.loc[:,["No1"]].values + data.loc[:,["No2"]].values + data.loc[:,["No3"]].values)/3
testx = x.flatten()
testy = y.flatten()

params = np.array([1,1])

def funcinv(x, a, b):
    x+=.00000000001
    return b + a * np.log(x)

res = curve_fit(funcinv, testx, testy, params, maxfev = 1000000)
print(res)

plt.plot(testx,testy, 'rs')
plt.plot(testx, -0.04010694 * np.log(testx) + 0.66274938 - .00000000001, 'bs')