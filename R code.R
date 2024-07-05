library(FactoMineR)
library(factoextra)
library(gridExtra)
library(survival)
library(ggrepel)
library(PCAmixdata)
df <- read.csv("C:/Users/SUJASH/OneDrive/Desktop/FAMD/pbc.csv")

# removed samples that are not randomized, and follow-up days:
df <- subset(df, drug != "not randomized") 
df <- subset(df, select = -fu.days)
df <- na.omit(df)

df <- subset(df, select = -stage)
df$spiders <- ifelse(df$spiders == 'present', 'spiders_y', 'spiders_n')
df$hepatom <- ifelse(df$hepatom == 'present', 'hepatom_y', 'hepatom_n')
df$ascites <- ifelse(df$ascites == 'present', 'ascites_y', 'ascites_n')
df$status <- ifelse(df$status == 1 , 'dead', 'alive')
df$drug <- ifelse(df$drug == 'D-penicillamine', 'd_pen', df$drug)
df$edema <- ifelse(df$edema == 'edema despite diuretic therapy', 'ed_wth',
                   ifelse(df$edema == 'edema, no diuretic therapy', 'ed_nth', 'n_ed'))

########################################################################################

res.afdm <- FAMD(df, graph = FALSE)
eig.val <- get_eigenvalue(res.afdm)
#plot(res.afdm,habillage = 4)
#fviz_screeplot(res.afdm, addlabels = TRUE, ylim = c(0, 25))

# Variable Contribution
a <- fviz_contrib(res.afdm, choice = "var", axes = 1)
b <- fviz_contrib(res.afdm, choice = "var", axes = 2)
#grid.arrange(a, b, nrow=1)

#Square Loading Plot
split <- splitmix(df)
res.pcamix <- PCAmix(X.quanti=split$X.quanti,  
                     X.quali=split$X.quali)
#plot(res.pcamix, choice="cor") 

###########################################################################################

#Kaplan-Meier survival curve for the entire dataset 
fit <- survfit(formula = Surv(df$fu.days, df$status == 1) ~ 1)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(fit, xlab = "Time", ylab = "Survival Probability", main = "Kaplan-Meier Estimator")

pc_fit <- survfit(formula = Surv(pc_df$fu.days, pc_df$status == 1) ~ 1)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(pc_fit, xlab = "Time", ylab = "Survival Probability", main = "Kaplan-Meier Estimator")

#mantel-Haenzel Test
survdiff(formula = Surv(df$fu.days, df$status == 1) ~ df$drug)

fit1 <- survfit(formula = Surv(df$fu.days, df$status == 1) ~ df$drug)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(fit1, col = c("blue", "red"), xlab = "Time (days)", ylab = "Survival Probability", main = "Difference in treatment variable")
legend("topright", legend = c("D-penicillamine", "Placebo"), col = c("blue", "red"), lty = 1)

############################################################################################
# K-fold cross validation 
set.seed(88)

k <- 4
folds <- cut(seq(1, nrow(df)), breaks = k, labels = FALSE)

train_concordance <- numeric(k)
test_concordance <- numeric(k)
test_se <- numeric(k)

for (i in 1:k) {
  test_indices <- which(folds == i, arr.ind = TRUE)
  train_indices <- setdiff(seq_len(nrow(df)), test_indices)
  
  train_df <- df[train_indices, ]
  test_df <- df[test_indices, ]
  
  m1 <- coxph(Surv(fu.days, status == 1) ~ ., data = train_df)
  
  train_con <- (concordance(m1))$concordance
  train_concordance[i] <- train_con
  
  test_con <- (concordance(m1, newdata = test_df))$concordance
  test_concordance[i] <- test_con
  
  predicted_probs <- predict(m1, newdata = test_df, type = "survival", se.fit = TRUE)
  
  test_se[i] <- mean(predicted_probs$se.fit)
}

mean_train_concordance <- mean(train_concordance)
mean_test_concordance <- mean(test_concordance)  
mean_test_se <- mean(test_se)   


zph_pc <- cox.zph(m1, transform = "identity") 
pvalue <- zph_pc$table['GLOBAL','p']  
par(mfrow = c(1,1)) 
#plot(zph_pc, main = 'Schoenfeld Individual Test', col = 'red', lines.col = 'black')

#########################################################################################

# stepwise with CV
My.stepwise.coxph <- function(Time, Status, variable.list, data, direction = "both", trace = FALSE) {
  surv_obj <- Surv(time = data[[Time]], event = data[[Status]])
  
  stepwise_model <- step(coxph(surv_obj ~ 1, data = data), 
                         scope = list(lower = as.formula("surv_obj ~ 1"), 
                                      upper = as.formula(paste("surv_obj ~", paste(variable.list, collapse = " + ")))),
                         direction = direction, trace = trace)
  
  return(stepwise_model)
}

my.variable.list <- setdiff(colnames(df), c("fu.days", "status"))

set.seed(98) 
num_folds <- 5 
concordance_vec <- numeric(num_folds)
concordance_se_vec <- numeric(num_folds)

for (i in 1:num_folds) {
  train_indices <- sample(seq_len(nrow(df)), size = 0.8 * nrow(df))
  train_df <- df[train_indices, ]
  test_df <- df[-train_indices, ]
  
  stepwise_model <- My.stepwise.coxph(Time = "fu.days", Status = "status", variable.list = my.variable.list, data = train_df)
  
  test_df$predicted_risk <- predict(stepwise_model, newdata = test_df, type = "risk")
  
  test_surv_obj <- Surv(time = test_df$fu.days, event = test_df$status)
  
  concordance <- survConcordance(test_surv_obj ~ test_df$predicted_risk)
  
  concordance_vec[i] <- concordance$concordance
  concordance_se_vec[i] <- concordance$std.err
}

mean_concordance <- mean(concordance_vec)
mean_concordance_se <- sqrt(sum(concordance_se_vec^2)) / num_folds

summary(stepwise_model)
model_aic <- AIC(stepwise_model)


zph_step_m <- cox.zph(stepwise_model, transform = "identity")
pvalue <- zph_step_m$table['GLOBAL', 'p']  
pvalue
test_df$predicted_risk <- predict(stepwise_model, newdata = test_df, type = "risk")

test_surv_obj <- Surv(time = test_df$fu.days, event = test_df$status)

concordance <- survConcordance(test_surv_obj ~ test_df$predicted_risk)

c_index <- concordance$concordance
c_index_se <- concordance$std.err

######################################################################################################

#Cox model with lasso
set.seed(1) 
train_indices <- sample(seq_len(nrow(df)), size = 0.8 * nrow(df))
train_df <- df[train_indices, ]
test_df <- df[-train_indices, ]

X_train <- model.matrix(~ . - fu.days - status, data = train_df)[, -1]
y_train <- Surv(train_df$fu.days, train_df$status)
X_test <- model.matrix(~ . - fu.days - status, data = test_df)[, -1]
y_test <- Surv(test_df$fu.days, test_df$status)

fit <- glmnet(X_train, y_train, family = "cox", alpha = 1)

cv.fit <- cv.glmnet(X_train, y_train, family = "cox", alpha = 1)
optimal_lambda <- cv.fit$lambda.min  

fit_lasso <- glmnet(X_train, y_train, family = "cox", alpha = 1, lambda = optimal_lambda)

coef <- fit_lasso$beta
nonzero_coef <- which(coef[, 1] != 0)
features <- rownames(coef)[nonzero_coef]

index_lasso <- which(names(train_df) %in% features)

selected_columns <- c("fu.days", "status", names(train_df)[index_lasso])
fit_lasso_final <- coxph(Surv(fu.days, status) ~ ., data = train_df[, selected_columns])

test_predictions <- predict(fit_lasso_final, newdata = test_df[, selected_columns], type = "risk")
concordance <- survConcordance(Surv(test_df$fu.days, test_df$status) ~ test_predictions)

c_index <- concordance$concordance  
c_index_se <- concordance$std.err  

predicted_probs <- predict(fit_lasso_final, newdata = test_df[, selected_columns], type = "survival", se.fit = TRUE)
mean_se_fit <- mean(predicted_probs$se.fit)

zph_lasso_final <- cox.zph(fit_lasso_final, transform = "identity")
pvalue <- zph_lasso_final$table['GLOBAL', 'p']  

################################################################################################

# No of PC's vs concordance
concordance_values <- numeric()

pc_range <- 1:30

for (ncp in pc_range) {
  res.famd <- FAMD(unl_df, ncp = ncp, graph = FALSE)
  
  pc_df <- as.data.frame(res.famd$ind$coord)
  
  pc_df$fu.days <- df$fu.days
  pc_df$status <- df$status
  
  crfit_pc <- coxph(Surv(fu.days, status == 1) ~ ., data = pc_df)
  
  concordance_values <- c(concordance_values, summary(crfit_pc)$concordance[1])
}

plot_data <- data.frame(
  PrincipalComponents = pc_range,
  Concordance = concordance_values
)

ggplot(plot_data, aes(x = PrincipalComponents, y = Concordance)) +
  geom_line() +
  geom_point() +
  labs(title = "Change in Concordance with Number of Principal Components",
       x = "Number of Principal Components",
       y = "Concordance Value") +
  theme_minimal()
##################################################################################################
# K fold CV on PC cox regression

unl_df <- subset(df, select = -status)
unl_df <- subset(df, select = -fu.days)
head(unl_df)

res.famd <- FAMD(unl_df,ncp = 14,graph = FALSE)
eig.val <- get_eigenvalue(res.famd)

loadings <- res.famd$var$coord

pc_df <- as.data.frame(res.famd$ind$coord)
pc_df$fu.days <- df$fu.days
pc_df$status <- df$status

set.seed(18)

folds <- cut(seq(1, nrow(pc_df)), breaks = k, labels = FALSE)

train_concordance <- numeric(k)
test_concordance <- numeric(k)
test_se <- numeric(k)

for (i in 1:k) {
  test_indices <- which(folds == i, arr.ind = TRUE)
  train_indices <- setdiff(seq_len(nrow(pc_df)), test_indices)
  
  train_df <- pc_df[train_indices, ]
  test_df <- pc_df[test_indices, ]
  
  m2 <- coxph(Surv(fu.days, status == 1) ~ ., data = train_df)
  
  train_con <- (concordance(m2))$concordance
  train_concordance[i] <- train_con
  
  test_con <- (concordance(m2, newdata = test_df))$concordance
  test_concordance[i] <- test_con
  
  predicted_probs <- predict(m2, newdata = test_df, type = "survival", se.fit = TRUE)
  
  test_se[i] <- mean(predicted_probs$se.fit)
}

mean_train_concordance <- mean(train_concordance)
mean_test_concordance <- mean(test_concordance)
mean_test_se <- mean(test_se)

zph_pc <- cox.zph(m2, transform = "identity") 
pvalue <- zph_pc$table['GLOBAL','p']