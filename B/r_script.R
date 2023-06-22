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


#function to extend the data with the grouping (if it is needed)
create_groups <- function(df){
  df <- df %>%
    mutate(
      TransferLearning = ifelse(model %in% c("B1", "B2", "B3"), 0, 1),
      TestDataType = case_when(
        TeD %in% c("TeD1", "TeD2", "TeD3", "TeD4") ~ "Classification",
        TeD == "TeD5" ~ "Recommendation",
        TeD %in% c("TeD6", "TeD7") ~ "Regression",
        TRUE ~ NA_character_
      ),
      TrDCount = rowSums(select(., starts_with("TrD")) == 1),
      RefineCount = case_when(
        model %in% c("B1", "B2", "B3", "MN") ~ 0,
        model %in% c("MS", "M1") ~ 1,
        model == "M2" ~ 2,
        model == "M3" ~ 3,
        model == "MF" ~ 4,
        TRUE ~ 0
      )
    )%>%
    select(-score, everything(), score)
  
  return(df)
}


#data loading and processing
df <- import("data.csv")
df <- create_groups(df)


#plots to explore the effects of model types and test data types
ggplot(df, aes(model, score)) + geom_point(aes(color =TeD ), position = position_jitter(width = 0.25), size=1)
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TransferLearning)
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TestDataType)
ggplot(df, aes(score))  + geom_histogram()  + facet_grid(~TransferLearning)


# analyzing model types
# m1 <- lm(score ~ model * TeD, df)
# m2 <- lm(score ~ TransferLearning * TeD, df)
# 
# emmeans(m1, ~model)
# emmeans(m2, ~TransferLearning)
# pairs(emmeans(m1, ~model))
# pairs(emmeans(m2, ~TransferLearning))


# Models - analyze training data effect
# 
# remove B1, B2, B3 models because they do not have training data.
# df_filtered <- subset(df, !(model %in% c("B1", "B2", "B3")))
# 
# m3 <- lm(score ~ model * TeD * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7), df_filtered)
# 
# summary(m3)
# as.data.frame(ref_grid(m3))
# pairs(emmeans(m3, ~TrD1))
# pairs(emmeans(m3, ~TrD2))
# pairs(emmeans(m3, ~TrD3))
# pairs(emmeans(m3, ~TrD4))
# pairs(emmeans(m3, ~TrD5))
# pairs(emmeans(m3, ~TrD6))
# pairs(emmeans(m3, ~TrD7))

# Q1

#complexer model (16128 individual effects) that considers the interactions among different TrDs
#m1 <-lm(score ~ model*TeD*TrD1*TrD2*TrD3*TrD4*TrD5*TrD6*TrD7*TrD8), df)

#simpler model (16128 individual effects) that doesn't consider these interactions
m2 <-lm(score ~ model*TeD*(TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7+TrD8), df)
summary(m2)

#code to remove NA values
m2$coefficients <- na.omit(m2$coefficients)

#save models
# save(m2, file="my_model2.rda")
#save(m1, file="my_model1.rda")

#load model
# load("my_model2.rda")
# load("my_model1.rda")

#EEM
eem_orig <- emmeans(m2, ~model, rg.limit = 20000)
# save(eem_orig, file="eem_orig_m1.rda")
# load("eem_orig_m1.rda")
print(eem_orig)
eem_orig_contrast <- pairs(eem_orig)
print(eem_orig_contrast)

#grouping results
eem_grouped <- add_grouping(eem_orig, "TransferLearning", "model", c("SL", "SL", "SL", "TL", "TL", "TL", "TL", "TL", "TL"))
eem_grouped <- emmeans(eem_grouped,  ~ TransferLearning)
eem_grouped_contrast <- pairs(eem_grouped)

#p values for score
print(test(eem_grouped))
#CI for score
print(eem_grouped)
#p value for contrast
print(eem_grouped_contrast)
#CI for contrast
print(confint(eem_grouped_contrast))

emmip(m2,  ~ model,  rg.limit = 20000)


# Q2

# individual effects of each TrD
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


# data exploration on the number of training datasets 

# TrDCount: numerical

# consider only TrDCount and score
# ggplot(df, aes(TrDCount, score)) + geom_smooth(se = FALSE, formula = y ~ s(x, bs = "cs", k = 4))

# see how scores are affected by TrDCount, taking into account TeD
ggplot(df, aes(TrDCount, score, color = TeD)) + geom_smooth(se = FALSE)

# see how scores are affected by TrDCount, taking into account TeD type
ggplot(df, aes(TrDCount, score, color = TestDataType)) + geom_smooth(se = FALSE, method = 'loess')

# see how scores are affected by TrDCount, taking into account Model
ggplot(df, aes(TrDCount, score, color = model)) + geom_smooth(se = FALSE)


# to see the TrDcount distributions over different TeD data sets
ggplot(df, aes(TrDCount)) + geom_histogram(binwidth = 5) + facet_grid(~TeD)


# TrDCount: catergorical

# df_subset <- subset(df, TrDCount %in% (2:7))

# the full model m3
df$TrDCount<- factor(df$TrDCount, levels = c(0:8))
df$RefineCount <- factor(df$RefineCount, levels = c(0:4))
# df2 <- df %>% group_by(TeD) %>% mutate_at(c('score'), ~(scale(.) %>% as.vector))
# find how RefineCount influences
ggplot(df, aes(TrDCount, score)) + geom_boxplot(notch = FALSE) + facet_wrap(~RefineCount)

m3 <-lm(score ~ model * TeD * TrDCount * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7) + RefineCount, df)
m3$coefficients <- na.omit(m3$coefficients)
summary(m3)
anova(m3)

# save(m3, file="full_model.rda")
# load("full_model.rda")


# find out the influence of the number of refined parts
m4 <-lm(score ~ model * TeD * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7) + RefineCount, df)
m4$coefficients <- na.omit(m4$coefficients)
summary(m4)
pairs(emmeans(m4, ~RefineCount, rg.limit = 362880))
emmeans(m4, ~RefineCount, rg.limit = 400000)

# find out the influence of the number of TrDs
m5 <-lm(score ~ model * TeD * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7) * TrDCount, df)
m5$coefficients <- na.omit(m5$coefficients)
summary(m5)
pairs(emmeans(m5, ~TrDCount, rg.limit = 362880))
emmeans(m5, ~TrDCount, rg.limit = 400000)

