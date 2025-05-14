# Diabetes-Prediction
The Pima Indian community near Phoenix, Arizona, has been studied since 1965 due to its high diabetes rates, with residents over age five taking regular health tests, including glucose tolerance tests. Diabetes is a long-term condition that raises blood sugar levels and can damage organs over time, so early diagnosis and access to treatment like insulin are vital. One way to help detect diabetes early is through statistical methods like discriminant analysis, which creates rules to classify new cases into groups. While standard discriminant analysis can struggle with outliers, robust versions use more reliable estimates, making them better suited for real-world data with varied group sizes and characteristics. Based on information above, diabetes needed an early detection as prevention to serious damage. At the same time, we have data from Pima Indian as sample for diabetes prediction. So, I would like to predict true positive rate of discriminant analysis and robust discriminant analysis.

# Data sampling and management
The dataset has been taken UCI Repository Machine Learning Databases but already converted to R format by Friedrich Leisch in mlbench package named PimaIndianDiabetes2.

# Data cleaning
![image](https://github.com/user-attachments/assets/080e8440-b6bc-4c18-a310-7559771cfb8c)

We can see from the original information of the dataset from [Smith et al. (1988)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2245318/), that the data has missing values. My decision to impute this data using mean or regression is based on [Buuren (2012)](https://stefvanbuuren.name/fimd/).

# EDA
![image](https://github.com/user-attachments/assets/19d6bb11-d4f1-46ba-9e72-a1cde74a1950)

We can see that several variables in the dataset exhibit outliers and skewed distributions. Meanwhile, the diabetes outcome variable is imbalanced, with approximately twice as many negative cases as positive ones. 
Discriminant analysis itself has assumptions that the features follow a normal multivariate distribution and between groups has equal matrices covariance. If normality assumption is violated, then it can be transform using Box-Cox. However, if the matrices covariance is unequal, then classification can be done with Quadratic discriminant analysis [Johnson and Wichern (2007)](https://books.google.co.id/books?id=gFWcQgAACAAJ).

# Data pre processing
In this section, I did compare mean vector (to know if discriminant analysis can be done or not) with [T^2 hotelling](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118391686), assess normal multivariate with QQ-plot and equal matrices covariance with BoxM, also balancing the dataset using [SMOTE](https://www.jair.org/index.php/jair/article/view/10302).
Result:
1. Mean vector is not the same, yes discriminant analysis can be done
2. Data seems not normal, so it being transform with box-cox
3. Matrices covariance not equal, so will be using quadratic discriminant analysis
4. SMOTE success make the data balance (1:1)

# Learning
After the data is processed, I'm going to compare 3 rules. 
First, quadratic discriminant rules for two groups (QDR) where I assume cost of misclassification is equal
![image](https://github.com/user-attachments/assets/3cd810ad-dde5-4dc0-b933-b4f5ee76b3cb)

second, quadratic discriminant score for several groups (can be used for two) (QDS)
![image](https://github.com/user-attachments/assets/a4bf985b-5f87-4398-a8b6-23be413b443a)

three, [robust quadratic discriminant score](https://www.researchgate.net/publication/245023580_Fast_and_Robust_Discriminant_Analysis) for several groups (can be used for two) using MCD estimator (RQDS)

![image](https://github.com/user-attachments/assets/b1bea4b5-b14f-4d4d-b56c-74fad51d82c2)

# Result and intepretation
| Model | Sensitivity |
|-------|-------------|
| QDR   | 86.87%      |
| QDS   | 86.87%      |
| RQDS  | 85.86%      |

Based on sensitivity results, QDR and QDS both achieved the highest sensitivity at 86.87%, while RQDS slightly lagged behind with 85.86%. Lagged is not an indication that robust methods are not needed. On the contrary this example shows that the performance of the robust approach is comparable with the classical one.

# Date
December 2024 - May 2025
