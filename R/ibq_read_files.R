#' read files form isobarquant workflow
#'
#' @param x the parent folder in which the the three result files are located
#' @return loads data into the global enironment of the R session alos prints the pooling files into the folders
#' @examples none yet
#' @export
#' @import tidyverse
ibq_read_results <- function(thefolder){

  alltxt_path <- list.files(path = thefolder, pattern =
                              "_proteins.txt|_summary.txt|_peptides.txt", recursive = TRUE)

  # generate path for reading ----------------------------------------

  results_files <-  data.frame(
    sample_name = alltxt_path %>%
      str_replace_all(pattern = "_proteins.txt|_summary.txt|_peptides.txt", "") %>%
      unique()
  )

  results_files <- results_files %>%
    mutate(
      prot_frame = map(sample_name, ~ paste(., "_proteins.txt", sep = "")),
      pep_frame = map(sample_name, ~ paste(., "_peptides.txt", sep = "")),
      summary_frame = map(sample_name, ~ paste(., "_summary.txt", sep = "")),
      bait_channel_frame = map(sample_name, ~ paste0(., "_pooling.csv"))
    )

  #  check if pooling.csv exists, write default ------------------------------
  default_pooling <-  tibble(
    "Bait_var_126",
    "Bait_var_127L",
    "Bait_var_127H",
    "Bait_var_128L",
    "Bait_var_128H",
    "Bait_var_129L",
    "Bait_var_129H",
    "Bait_var_130L",
    "Bait_var_130H",
    "Bait_var_131L"
  )


  default_pooling[2, ] <-  c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)



  colnames(default_pooling) <-  c(
    "signal_sum_126",
    "signal_sum_127L",
    "signal_sum_127H",
    "signal_sum_128L",
    "signal_sum_128H",
    "signal_sum_129L",
    "signal_sum_129H",
    "signal_sum_130L",
    "signal_sum_130H",
    "signal_sum_131L"
  )

  map(results_files[, "bait_channel_frame"],
      ~ ifelse(file.exists(.x) , "yes", write_csv(default_pooling, .x)))

  shortnamefilter <- ""
  results <- results_files %>% transmute(
    short_name = map_chr(sample_name, ~ gsub(shortnamefilter, "", .)),
    prot = map(
      prot_frame,
      ~ read_tsv(., col_types = "dcccdddddddddddddddddddddddddddddddddddddddddd")
    ),
    pep = map(
      pep_frame,
      ~ read_tsv(., col_types = "dcccddddddcddcdddddddddddddddddddddddd")
    ),
    bait_channel = map(
      bait_channel_frame,
      ~ read_csv(.)
    ),
    run_smry = map(summary_frame, ~ read_tsv(.))
  ) %>% as_tibble



  results <<- results %>%
    mutate(run = short_name %>%
             str_extract(., "^.*./") %>%
             str_remove_all(., "/")
    )
}
