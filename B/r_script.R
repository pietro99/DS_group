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



#data processing

df <- import("data.csv")

#dummy <- dummyVars(" ~ .", data=d)
#df <- data.frame(predict(dummy, newdata=d))


# Define the groups and their corresponding values
# group_mappings <- list(
#   TransferLearning = c("modelM1","modelM2", "modelM3", "modelMS", "modelMN", "modelMF"),
#   StandardLearning = c("modelB1", "modelB2", "modelB3"),
#   Classification = c("TeDTeD1", "TeDTeD2", "TeDTeD3", "TeDTeD4"),
#   Reccomandation = c("TeDTeD5"),
#   Regression = c("TeDTeD6", "TeDTeD7")
# )
# 
# 
# 
# for (group_name in names(group_mappings)) {
#   group_columns <- group_mappings[[group_name]]
# 
#   # Create a new column for the group based on if any of the group columns contain 1
#   df <- df %>%
#     mutate(!!group_name := as.integer(rowSums(select(., all_of(group_columns))) > 0))
# }
# 
# write.csv(df, "./data_processed.csv")


#Models

m <- lm(score ~ model*TeD*(TrD1+TrD2+TrD3+TrD4+TrD5+TrD6+TrD7), df)

summary(m)
as.data.frame(ref_grid(m))

pairs(emmeans(m, ~model))


