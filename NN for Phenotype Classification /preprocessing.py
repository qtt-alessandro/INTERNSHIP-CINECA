import seaborn as sns
import pandas as pd
import numpy as np
import time
seed = 14
np.random.seed(seed)
from keras import backend as K, activations
from keras import layers
from keras import regularizers
from keras.models import Sequential, Model
from keras.models import load_model
from keras.layers import Dense, Activation, Input
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import recall_score, precision_score, f1_score
from matplotlib import pyplot as plt
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasClassifier
import tensorflow as tf
import matplotlib.pyplot as plt
# from tf_explain.callbacks.activations_visualization import ActivationsVisualizationCallback
import sklearn

myTrait = 'stimulation'  # IFNb, LPS, dNS1, unstim

features_table = pd.read_table('Input/features.txt')
features = features_table.set_index('geo_accession')
features = np.transpose(features)
GE_table = pd.read_table('Input/GE.txt')
GE = GE_table.set_index('ID_REF')
GE = np.transpose(GE)
normalize = True
# normalize = False
if normalize:
    GE = (GE - GE.mean())/(GE.std()) # normalizing GE by z-score for each gene
    print("Normalized GE")
else:
    print("Original GE")


dev = GE.sample(frac=0.30,replace=False, axis = 0, random_state = 1)
sample_patient_id = list(dev.index)
patient_stimualtion = []
for patient in sample_patient_id:
    patient_stimualtion.append(features.loc[patient].stimulation)

dev_arr = []
for index, row in dev.iterrows():
    dev_arr.append(np.array(row))
deviation = np.array(dev_arr)

num_genes = GE.shape[1]
num_samples = GE.shape[0]

print('Number of genes: %i' %num_genes)
print('Number of samples: %i' %num_samples)

X = GE.values.astype(float)

Y = features[myTrait]

classes, counts = np.unique(Y, return_counts=True)

num_classes = classes.size

print(dict(zip(classes, counts)))

encoder = LabelEncoder()
encoder.fit(Y)
encodedY = encoder.transform(Y)
dummyY = np_utils.to_categorical(encodedY)
