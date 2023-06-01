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



#plots
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TransferLearning)
ggplot(df, aes( y=score)) + geom_boxplot(notch = FALSE) + facet_wrap(~TestDataType)
ggplot(df, aes(score))  + geom_histogram()  + facet_grid(~TransferLearning)



#analyzing model types
m1 <- lm(score ~ model * TeD, df)
m2 <- lm(score ~ TransferLearning * TeD, df)

summary(m1)
as.data.frame(ref_grid(m1))
emmeans(m1, ~model)
emmeans(m2, ~TransferLearning)
pairs(emmeans(m2, ~model))
pairs(emmeans(m2, ~TransferLearning))

#Models - analyze training data effect

#remove B1, B2, B3 models because they do not have training data.
df_filtered <- subset(df, !(model %in% c("B1", "B2", "B3")))

m3 <- lm(score ~ model * TeD * (TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7), df_filtered)



summary(m3)
as.data.frame(ref_grid(m3))
pairs(emmeans(m3, ~TrD1))
pairs(emmeans(m3, ~TrD2))
pairs(emmeans(m3, ~TrD3))
pairs(emmeans(m3, ~TrD4))
pairs(emmeans(m3, ~TrD5))
pairs(emmeans(m3, ~TrD6))
pairs(emmeans(m3, ~TrD7))

#full model (5376 individual effects)
m1 <-lm(score ~ model*TeD*TrD1*TrD2*TrD3*TrD4*TrD5*TrD6*TrD7, df_filtered)

summary(m1)
as.data.frame(ref_grid(m1))

pairs(emmeans(m1, ~model))
