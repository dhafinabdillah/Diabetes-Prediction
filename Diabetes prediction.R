# Data sampling and management
library(mlbench)
data("PimaIndiansDiabetes2")
Diabetes <- PimaIndiansDiabetes2
head(Diabetes)

#Data cleaning
colSums(is.na(Diabetes))

glucose_fill_value <- mean(Diabetes$glucose[Diabetes$glucose >= 107.2 & Diabetes$glucose <= 123.1], na.rm = TRUE)
Diabetes$glucose[is.na(Diabetes$glucose)] <- glucose_fill_value

mass_fill_value <- mean(Diabetes$mass[Diabetes$mass >= 26.841 & Diabetes$mass <= 33.55], na.rm = TRUE)
Diabetes$mass[is.na(Diabetes$mass)] <- mass_fill_value

mean_pressure <- mean(Diabetes$pressure, na.rm = TRUE)
Diabetes$pressure[is.na(Diabetes$pressure)] <- mean_pressure

library(mice)
insulin <- Diabetes[, c("insulin", "glucose")]
imp <- mice(insulin, method = "norm.predict", seed = 1,
            m = 1, print = FALSE)
xyplot(imp, insulin ~ glucose)

triceps <- Diabetes[, c("triceps", "mass")]
imp1 <- mice(triceps, method = "norm.predict", seed = 1,
             m = 1, print = FALSE)
xyplot(imp1, triceps ~ mass)

completed_data <- complete(imp)
completed_data1 <- complete(imp1)

Diabetes$insulin <- completed_data$insulin
Diabetes$triceps <- completed_data1$triceps

# EDA
par(mfrow = c(3, 3))

boxplot(Diabetes$pregnant, main = "Pregnant")
boxplot(Diabetes$glucose, main = "Glucose")
boxplot(Diabetes$pressure, main = "Pressure")
boxplot(Diabetes$triceps, main = "Triceps")
boxplot(Diabetes$insulin, main = "Insulin")
boxplot(Diabetes$mass, main = "Mass")
boxplot(Diabetes$pedigree, main = "Pedigree")
boxplot(Diabetes$age, main = "Age")
barplot(table(Diabetes$diabetes), main = "Diabetes")

# Pre processing
Diabetes$diabetes <- as.character(Diabetes$diabetes)

Diabetes$diabetes[Diabetes$diabetes == "neg"] <- "1"
Diabetes$diabetes[Diabetes$diabetes == "pos"] <- "2"

x1 <- Diabetes[Diabetes$diabetes == "1", c("pregnant", "glucose","pressure","triceps","insulin","mass","pedigree","age")]
x2 <- Diabetes[Diabetes$diabetes == "2", c("pregnant", "glucose","pressure","triceps","insulin","mass","pedigree","age")]

x <- Diabetes[, 1:8]
y <- Diabetes$diabetes

n <- nrow(x)
p <- as.numeric(ncol(x))
g <- as.numeric(length(unique(y)))

mean_x1 <- colMeans(x1)
mean_x2 <- colMeans(x2)

cov_x1 <- cov(x1)
cov_x2 <- cov(x2)

n1 <- nrow(x1)
n2 <- nrow(x2)

spooled <- ((n1 - 1) * cov_x1 + (n2 - 1) * cov_x2) / (n1 + n2 - 2) # Calculate the pooled covariance matrix

mean_diff <- as.matrix(mean_x1 - mean_x2)

t2 <- (n1 * n2) / (n1 + n2) * t(mean_diff) %*% solve(spooled) %*% mean_diff
t2

F_critical <- qf(0.95, 8, 759)

T2_critical <- (8 * (768 - 2)) / (768 - 8 - 1) * F_critical
T2_critical # Comparing mean vectors

sd<-mahalanobis(x,colMeans(x),cov(x))

chikuadrat <- qchisq((n - seq_len(n) + 0.5) / n, 
                     df = p)

par(mfrow = c(1, 1))
qqplot(chikuadrat,sd,main="Q-Q plot",ylab="Squared Distance",xlab="Chi-square")
abline(a=0,b=1) # Normal multivariate

library(MASS)
x <- as.data.frame(lapply(x, function(col) as.numeric(as.character(col))))

if (any(x <= 0)) {
  shift <- abs(min(x)) + 1
  x <- x + shift
  message("Data was shifted by ", shift, " to make all values > 0.")
}

x_bc_list <- list()

for (i in seq_along(x)) {
  col_data <- x[[i]]
  model <- lm(col_data ~ 1)
  bc_result <- boxcox(model, lambda = seq(-2, 2, 0.1))
  lambda_opt <- bc_result$x[which.max(bc_result$y)]
  
  
  if (abs(lambda_opt) < 1e-6) {
    transformed <- log(col_data)
  } else {
    transformed <- (col_data^lambda_opt - 1) / lambda_opt
  }
  
  x_bc_list[[colnames(x)[i]]] <- transformed
} # Box cox transformation

x_bc <- as.data.frame(x_bc_list)

sd <- mahalanobis(x_bc, colMeans(x_bc), cov(x_bc))
chikuadrat <- qchisq((n - seq_len(n) + 0.5) / n, df = p)

qqplot(chikuadrat, sd, main = "Q-Q plot after Box-Cox")
abline(a=0, b=1) # Normal multivariate after transformation

library(dplyr)
normal_multivariate <- data.frame(sd)
normal_multivariate$chikuadrat <- chikuadrat
normal_multivariate$results<-ifelse(normal_multivariate$sd <= normal_multivariate$chikuadrat, 'True', 
                                    ifelse(normal_multivariate$chikuadrat > normal_multivariate$sd, 'No', 'None'))

sum_sd_less_than_chi_Diabetes<-length(which(normal_multivariate$results=="True"))
n_Diabetes<-nrow(normal_multivariate)
(sum_sd_less_than_chi_Diabetes/n_Diabetes)*100 # Normal multivariate after transformation

library(biotools)
boxM <- boxM(as.matrix(x_bc), y)
boxM

v <- 0.5 * p * (p + 1) * (g - 1)
chisq<-qchisq(c(0.05),df=v,lower.tail=FALSE)
chisq # Equality matrices covariance

x_bc$diabetes <- Diabetes$diabetes
x_bc$diabetes <- as.factor(x_bc$diabetes)

library(DMwR)
set.seed(123)
x_bc_smote <- SMOTE(diabetes ~ ., x_bc, perc.over = 87, perc.under = 215, k = 5)

table(x_bc_smote$diabetes)

# Classification for unequal two populations
s <- round(nrow(x_bc_smote)*0.8)
set.seed(12345)
sample <- sample(seq_len(nrow(x_bc_smote)), size=s)

train_diabetes <- x_bc_smote[sample, ]
test_diabetes <- x_bc_smote[-sample, ]

x1_train <- train_diabetes[train_diabetes$diabetes == "1", c("pregnant", "glucose","pressure","triceps","insulin","mass","pedigree","age")]
x2_train <- train_diabetes[train_diabetes$diabetes == "2", c("pregnant", "glucose","pressure","triceps","insulin","mass","pedigree","age")]

x_test <- as.matrix(test_diabetes[, c("pregnant", "glucose","pressure","triceps","insulin","mass","pedigree","age")])
y_test <- test_diabetes$diabetes

mean_x1_train <- colMeans(x1_train)
mean_x2_train <- colMeans(x2_train)

cov_x1_train <- cov(x1_train)
cov_x2_train <- cov(x2_train)

det_Sigma1_train <- det(cov_x1_train)
det_Sigma2_train <- det(cov_x2_train)

inv_Sigma1_train <- solve(cov_x1_train)
inv_Sigma2_train <- solve(cov_x2_train)

term1 <- 0.5 * log(det_Sigma1_train / det_Sigma2_train)
term2 <- 0.5 * (t(mean_x1_train) %*% inv_Sigma1_train %*% mean_x1_train - 
                  t(mean_x2_train) %*% inv_Sigma2_train %*% mean_x2_train)
k <- term1 + term2

p1 <- nrow(x1_train) / (nrow(x1_train) + nrow(x2_train))
p2 <- 1 - p1
c12 <- 1
c21 <- 1
threshold <- log((c12 / c21) * (p2 / p1))

compute_score <- function(x) {
  term_quad <- -0.5 * t(x) %*% (inv_Sigma1_train - inv_Sigma2_train) %*% x
  term_linear <- (t(mean_x1_train) %*% inv_Sigma1_train - t(mean_x2_train) %*% inv_Sigma2_train) %*% x
  score <- term_quad + term_linear - k
  return(score)
}

scores <- apply(x_test, 1, compute_score)
predicted_classes <- ifelse(scores >= threshold, "1", "2")

classification_df <- data.frame(
  score = as.numeric(scores),
  decision_boundary = threshold,
  predicted_class = predicted_classes,
  actual_class = y_test
)

head(classification_df)

library(caret)
confusionMatrix(as.factor(predicted_classes), as.factor(y_test), positive = "2")

# Classification for unequal several populations 
n1 <- nrow(x1_train)
n2 <- nrow(x2_train)
total_n <- n1 + n2

log_p1 <- log(n1 / total_n)
log_p2 <- log(n2 / total_n)

qda_score <- function(x, mean_x, inv_sigma, det_sigma, log_p) {
  d_quad <- -0.5 * log(det_sigma) -
    0.5 * t(x - mean_x) %*% inv_sigma %*% (x - mean_x) +
    log_p
  return(d_quad)
}

classification_results <- apply(x_test, 1, function(x) {
  score_x1 <- qda_score(as.numeric(x), mean_x1_train, inv_Sigma1_train, det_Sigma1_train, log_p1)
  score_x2 <- qda_score(as.numeric(x), mean_x2_train, inv_Sigma2_train, det_Sigma2_train, log_p2)
  predicted_x <- ifelse(score_x1 > score_x2, "1", "2")
  return(c(predicted_x, score_x1, score_x2))
})

classification__df <- data.frame(
  actual_class = y_test,
  predicted_class = classification_results[1, ],
  score_x1 = as.numeric(classification_results[2, ]),
  score_x2 = as.numeric(classification_results[3, ])
)

head(classification__df)

confusionMatrix(as.factor(classification__df$predicted_class), as.factor(y_test), positive = "2")

# Classification for unequal several populations with robust estimator 
library(rrcov)
mcd1 <- CovMcd(x1_train, alpha = 0.75, nsamp = 500)
mcd2 <- CovMcd(x2_train, alpha = 0.75, nsamp = 500)

mean_mcd_x1 <- mcd1@center
mean_mcd_x2 <- mcd2@center

cov_mcd_x1 <- mcd1@cov
cov_mcd_x2 <- mcd2@cov

det_mcd_Sigma1 <- det(cov_mcd_x1)
det_mcd_Sigma2 <- det(cov_mcd_x2)

inv_mcd_Sigma1 <- solve(cov_mcd_x1)
inv_mcd_Sigma2 <- solve(cov_mcd_x2)

n1_mcd <- sum(mcd1@raw.wt)
n2_mcd <- sum(mcd2@raw.wt)
total_n_mcd <- n1_mcd + n2_mcd

log_p1_mcd <- log(n1_mcd / total_n_mcd)
log_p2_mcd <- log(n2_mcd / total_n_mcd)

qda_score_mcd <- function(x, mean_mcd_x, inv_mcd_sigma, det_mcd_sigma, log_p) {
  d_quad_mcd <- -0.5 * log(det_mcd_sigma) -
    0.5 * t(x - mean_mcd_x) %*% inv_mcd_sigma %*% (x - mean_mcd_x) +
    log_p
  return(d_quad_mcd)
}

classification_results_mcd <- apply(x_test, 1, function(x) {
  score_mcd_x1 <- qda_score_mcd(as.numeric(x), mean_mcd_x1, inv_mcd_Sigma1, det_mcd_Sigma1, log_p1_mcd)
  score_mcd_x2 <- qda_score_mcd(as.numeric(x), mean_mcd_x2, inv_mcd_Sigma2, det_mcd_Sigma2, log_p2_mcd)
  predicted_mcd_x <- ifelse(score_mcd_x1 > score_mcd_x2, "1", "2")
  return(c(predicted_mcd_x, score_mcd_x1, score_mcd_x2))
})

classification_mcd_df <- data.frame(
  actual_class = y_test,
  predicted_class = classification_results_mcd[1, ],
  score_x1 = as.numeric(classification_results_mcd[2, ]),
  score_x2 = as.numeric(classification_results_mcd[3, ])
)

print(head(classification_mcd_df))

confusionMatrix(as.factor(classification_mcd_df$predicted_class), as.factor(y_test), positive = "2")
