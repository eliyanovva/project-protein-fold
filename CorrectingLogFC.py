import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

data = pd.read_csv('Ligands_withSMILE/pS6_DE_1p_2heptanone.csv')
y = np.abs(np.array(data.loc[:,["logFC"]]))
x = (np.array(data.loc[:,["No1"]]) + np.array(data.loc[:,["No2"]]) + np.array(data.loc[:,["No3"]]))/3

model = np.polyfit(x, y)
r_sq = model.score(x, y)
print(f"coefficient of determination: {r_sq}")

plt.plot(x, y)