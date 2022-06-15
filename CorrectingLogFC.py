import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

data = pd.read_csv('e3fp/pS6_DE_1p_2heptanone.csv')
y = np.abs(data.loc[:,["logFC"]].values)
x = (data.loc[:,["No1"]].values + data.loc[:,["No2"]].values + data.loc[:,["No3"]].values)/3
testx = x.flatten()
testy = y.flatten()

model = np.polyfit(testx, testy, 3)

plt.plot(testx, testy, 'bs')
plt.plot(model[0] + model[1]*testx + model[2]*testx*testx + model[3]*testx*testx*testx, 'rs')
print(model)

rss = 0
for i in range(0, len(y), 1):
    rss += (testy[i] - (model[0] + model[1]*testx[i] + model[2]*testx[i]*testx[i] + model[3]*testx[i]*testx[i]*testx[i]))*(testy[i] - (model[0] + model[1]*testx[i] + model[2]*testx[i]*testx[i] + model[3]*testx[i]*testx[i]*testx[i]))
print(rss)