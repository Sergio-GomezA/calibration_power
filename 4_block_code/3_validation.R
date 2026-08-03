mc <- ifelse(local_run, 1, available_cores())

pow_threshold_label <- as.character(pow_threshold) %>% gsub("\\.", "_", .)


# 1. data preparation ####
cat("--------------------------------------------------------------------\n")
cat("Preparing data for low wind events analysis\n")
cat("--------------------------------------------------------------------\n")
