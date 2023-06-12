#install required packages

#install.packages("rio") 
#install.packages("dplyr") 
#install.packages("emmeans") 
#install.packages("ggplot2") 
#install.packages("caret")

#import required packages


library("rio")
library("dplyr")
library("emmeans")
library("ggplot2")
library("caret")

#function to enhence the data with the grouping (if it is even needed)
create_groups <- function(df){
  df <- df %>%
    mutate(
      TransferLearning = ifelse(model %in% c("B1", "B2", "B3"), 0, 1),
      TestDataType = case_when(
        TeD %in% c("TeD1", "TeD2", "TeD3", "TeD4") ~ "Classification",
        TeD == "TeD5" ~ "Reccomandation",
        TeD %in% c("TeD6", "TeD7") ~ "Regression",
        TRUE ~ NA_character_
      ),
      TrDCount = rowSums(select(., starts_with("TrD")) == 1)
      
      
    )%>%
    select(-score, everything(), score)
  
  return(df)
}

#data processing
df <- import("data.csv")
df <- create_groups(df)


ggplot(df, aes(model, score)) + geom_point(aes(color =TeD ), position = position_jitter(width = 0.25), size=1)
#plots
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TransferLearning)
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TestDataType)
ggplot(df, aes(score))  + geom_histogram()  + facet_grid(~TransferLearning)




#analyzing model types
#m1 <- lm(score ~ model * TeD, df)
# m2 <- lm(score ~ TransferLearning * TeD, df)
# m2 <- lm(score ~ TransferLearning * TeD, df)

#emmeans(m1, ~model)
# emmeans(m2, ~TransferLearning)
# emmeans(m2, ~TransferLearning)
#pairs(emmeans(m1, ~model))
# pairs(emmeans(m2, ~TransferLearning))
# pairs(emmeans(m2, ~TransferLearning))

#contrast_model <- contrast(emmeans_m1, list(StandardLearning = c("B1", "B2", "B3"), TransferLearning = c("M1", "M2", "M2", "MF", "MN", "MS")))
#contrast_m1 <- contrast(emmeans_m1, list(StandardLearning = c("B1", "B2", "B3"), TransferLearning = c("M1", "M2", "M2", "MF", "MN", "MS")))

#Models - analyze training data effect

#remove B1, B2, B3 models because they do not have training data.
# df_filtered <- subset(df, !(model %in% c("B1", "B2", "B3")))
# 
m3 <- lm(score ~ model * TeD * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7), df)
m3$coefficients <- na.omit(m3$coefficients)
summary(m3)
# as.data.frame(ref_grid(m3))
# pairs(emmeans(m3, ~TrD1))
# pairs(emmeans(m3, ~TrD2))
# pairs(emmeans(m3, ~TrD3))
# pairs(emmeans(m3, ~TrD4))
# pairs(emmeans(m3, ~TrD5))
# pairs(emmeans(m3, ~TrD6))
# pairs(emmeans(m3, ~TrD7))



#full model (16128 individual effects)
#m1 <-lm(score ~ model*TeD*TrD1*TrD2*TrD3*TrD4*TrD5*TrD6*TrD7*TrD8), df)

#simpler model (16128 individual effects)
m2 <-lm(score ~ model*TeD*(TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7+TrD8), df)
summary(m2)

# simplest model for RQ1
# m4 <- lm(score ~ model*TeD, df)

# take the number of datasets into account
# TrDCount: numerical
ggplot(df, aes(TrDCount, score)) + geom_smooth(se = FALSE, formula = y ~ s(x, bs = "cs", k = 4))


ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TeD)

library(gridExtra)

# Plot 1
p1 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD1, scales = "free")

# Plot 2
p2 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD2, scales = "free")

# Plot 3
p3 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD3, scales = "free")

# Plot 4
p4 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD4, scales = "free")

# Plot 5
p5 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD5, scales = "free")

# Plot 6
p6 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD6, scales = "free")

# Plot 7
p7 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD7, scales = "free")

# Plot 8
p8 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD8, scales = "free")

# Combine plots into one figure
combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3, ncol = 3)

# Display the combined plot
print(combined_plot)

# TrDCount: catergorical
df_subset <- subset(df, TrDCount %in% (2:7))
df_subset$TrDCount<- factor(df_subset$TrDCount, levels = c(2:7))
m3 <-lm(score ~ model * TrDCount * TeD, df_subset)
m3$coefficients <- na.omit(m3$coefficients)
summary(m3)

# find out the influence of the number of training datasets
pairs(emmeans(m3, ~TrDCount, rg.limit = 150000))
m2 <-lm(score ~ model*TeD*(TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7+TrD8), df)
summary(m2)

# simplest model for RQ1
# m4 <- lm(score ~ model*TeD, df)

# take the number of datasets into account
# TrDCount: numerical
ggplot(df, aes(TrDCount, score)) + geom_smooth(se = FALSE, formula = y ~ s(x, bs = "cs", k = 4))


ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TeD)

library(gridExtra)

# Plot 1
p1 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD1, scales = "free")

# Plot 2
p2 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD2, scales = "free")

# Plot 3
p3 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD3, scales = "free")

# Plot 4
p4 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD4, scales = "free")

# Plot 5
p5 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD5, scales = "free")

# Plot 6
p6 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD6, scales = "free")

# Plot 7
p7 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD7, scales = "free")

# Plot 8
p8 <- ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TrD8, scales = "free")

# Combine plots into one figure
combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3, ncol = 3)

# Display the combined plot
print(combined_plot)

# TrDCount: catergorical
df_subset <- subset(df, TrDCount %in% (2:7))
df_subset$TrDCount<- factor(df_subset$TrDCount, levels = c(2:7))
m3 <-lm(score ~ model * TrDCount * TeD, df_subset)
m3$coefficients <- na.omit(m3$coefficients)
summary(m3)

# find out the influence of the number of training datasets
pairs(emmeans(m3, ~TrDCount, rg.limit = 150000))

#load model
load("my_model2.rda")
# load("my_model1.rda")

#code to remove NA values
m2$coefficients <- na.omit(m2$coefficients)

#save models
# save(m2, file="my_model2.rda")
#save(m1, file="my_model1.rda")

#EEM
eem_orig <- emmeans(m2, ~model, rg.limit = 20000)
save(eem_orig, file="eem_orig_m1.rda")
load("eem_orig_m1.rda")
print(eem_orig)
#grouping results
eem_grouped <- add_grouping(eem_orig, "TransferLearning", "model", c("SL", "SL", "SL", "TL", "TL", "TL", "TL", "TL", "TL"))
eem_grouped <- emmeans(eem_grouped,  ~ TransferLearning)
eem_grouped_contrast <- pairs(eem_grouped)

print(summary(m2))
#p values for score
print(test(eem_grouped))
#CI for score
print(eem_grouped)
#p value for contrast
print(eem_grouped_contrast)
#CI for contrast
print(confint(eem_grouped_contrast))


# Q2
e <- rbind(
  confint(pairs(emmeans(m2, ~TrD1, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD2, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD3, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD4, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD5, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD6, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD7, rg.limit = 20000))),
  confint(pairs(emmeans(m2, ~TrD8, rg.limit = 20000)))
)

print(e)
#results plot

ggplot(summary(confint(contrast)),  aes(estimate, y = contrast) ) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL)) +
  labs(x = "score", y = "model")

ggplot(summary(confint(eem_grouped)),  aes(emmean, y = Group) ) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL)) +
  labs(x = "score", y = "TeD Group")

ggplot(summary(e),  aes(estimate, y = contrast) ) +
  geom_point() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL)) +
  labs(x = "score", y = "model")


pwpp(eem, type = "response")

plot(eem, comparisons = TRUE)

emmip(m2,  ~ model,  rg.limit = 20000)



