#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI

#Imports
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

#Number of trees in the forest
n_estimators = [int(x) for x in np.linspace(10, 100, 10)]
#Features considered at each split
max_features = ['auto', 'sqrt']
