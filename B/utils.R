
create_one_hot <- function(df){
  group_mappings <- list(
    TransferLearning = c("modelM1","modelM2", "modelM3", "modelMS", "modelMN", "modelMF"),
    StandardLearning = c("modelB1", "modelB2", "modelB3"),
    Classification = c("TeDTeD1", "TeDTeD2", "TeDTeD3", "TeDTeD4"),
    Reccomandation = c("TeDTeD5"),
    Regression = c("TeDTeD6", "TeDTeD7")
  )
  
  
  
  for (group_name in names(group_mappings)) {
    group_columns <- group_mappings[[group_name]]
    
    # Create a new column for the group based on if any of the group columns contain 1
    df <- df %>%
      mutate(!!group_name := as.integer(rowSums(select(., all_of(group_columns))) > 0))
  }
  return(df)
}
