# Breast Cancer Polygenic Risk Score Analysis

## Overview
Reproducible R and Python pipeline analysing publicly available GWAS summary statistics from Michailidou et al. (2017, Nature), n = 122,977 cases, 105,974 controls. Calculating simulated PRS using logistic regression from the effect sizes obtained from Michailidou et al, 2017 achieving AUC = 0.542, consistent with modest predictive power expected for complex polygenic traits.

## Results
![ROC Curve](results/py_roc_curve.png)
![Manhattan Plot](results/manhattan_plot.png)

## R Pipeline
The R script preprocessed the GWAS data and obtained the descriptive statistics and summary plots such as Manhattan plot, QQ plot, effect size distribution plots to further export a clean SNP table

## Python Pipeline
The python notebook visualized the GWAS effect size landscape, simulated PRS using real GWAS effect sizes, performed logistic regression classification of cases vs controls and evaluated the accuracy of the model using AUC-ROC curve.

## Limitations
-obtain individual genetic data instead of simulation
-perform LD clumping and thresholding 
-an ethnically diverse population, not limited to European ancestry

## Data
Michailidou et al. (2017), Nature, 551:92-94.
Downloaded from: https://bcac.ccge.medschl.cam.ac.uk

## Skills demonstrated
R · Python · GWAS analysis · Quality control · Polygenic risk scores · Logistic regression · Scikit-learn · ggplot2 · Data visualisation
