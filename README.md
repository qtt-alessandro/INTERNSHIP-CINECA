# Biological Interpretability of Deep Neural Network for Cancer Type Prediction

<p align="center">
  <img src="https://github.com/qtt-alessandro/INTERNSHIP-CINECA/blob/main/README_img/logo.jpeg" />
</p>

<p align="center">
Faculty of Information Engineering, Informatics and Statistics <br />
UNIVERSITY OF "LA SAPIENZA" <br />
MSC IN DATA SCIENCE <br />
</p>

## Abstract 
This report presents some results obtained regarding the interpretability of using neural
networks for classification, particularly based on phenotype and cancer types from
gene expression values. The paper considers some of the best-known techniques for interpreting
deep neural networks and other techniques we have tested.

The datasets and part of the code are taken from two recently published articles
1. “Exploring Neural Networks and Related Visualization Techniques in Gene Expression Data" (Roni Wilentzik Müller and Irit Gat-Viks)
2. "Convolutional neural network models for cancer type prediction based on gene expression" (Milad Mostavi, Yu-Chiao Chiu, Yufei Huang and Yidong Chen) 

For the interpretation based on the first article two techniques are
compared: activation maximization and most weighted paths determination.
Thanks to this approach we can detect the two paths most frequently used for the classification given an input feature.
This method is compared with the one presented in the article and although those two methods agree in finding the most
activated genes for each class, maximizing activation can result in more significant scoring.
Regarding the second article, following the computer vision approach presented for the cancer type classification,
several gradient-based and activation maximization methods are developed. 
Finally, after identifying the genes that most determined the classification of a cancer type to a given category, the pathways are enriched and the statistical significance of the genes tested.  




<p align="right">
<b>Author:</b> Alessandro Quattrociocchi <br />
 <b>Supervisor:</b> Silvia Gioiosa <br />
 <b>Supervisor:</b> Bhaskar Agarwal<br />
</p>

