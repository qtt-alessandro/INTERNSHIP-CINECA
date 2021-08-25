import tensorflow as tf
from tensorflow.keras import backend as K
from tf_keras_vis.saliency import Saliency
from tf_keras_vis.utils import normalize
from tf_keras_vis.utils.model_modifiers import ReplaceToLinear
from tf_keras_vis.utils.scores import CategoricalScore
import numpy as np


def model_modifier_function(cloned_model):
    cloned_model.layers[-1].activation = tf.keras.activations.linear


def vanilla_gradient(model,cancer_class,cancer_class_idx,features_per_class):

    replace2linear = ReplaceToLinear()
    score = CategoricalScore([cancer_class_idx])

    saliency = Saliency(model,
                        model_modifier=replace2linear,
                        clone=True)

    A = np.zeros((71,100))
    for img in (features_per_class[cancer_class]):
        saliency_map = saliency(score, img)
        A+= normalize(saliency_map.squeeze())
    return normalize(A)


def call_gradient_given_class(model,cancer_class,cancer_class_idx,gradient_selection,features_per_class):
    if gradient_selection == "vanilla":
        grd = vanilla_gradient(model,cancer_class,cancer_class_idx,features_per_class)
    return grd
