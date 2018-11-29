#' this function collects summary information
#' @export
ibq_do_ms_identification_summary <- function(result){

# adding row_sums for protein and peptide table----------------------------
results$prot <-
  results$prot %>% map(~ mutate(., total_signal = rowSums(select(
    ., contains("signal_sum")
  ))))
results$pep <-
  results$pep %>% map(~ mutate(., total_signal = rowSums(select(., contains(
    "sig_"
  )))))

# add relative retention times --------------------------------------------
results$pep <- results$pep %>%
  modify_depth(.,1, ~mutate(., rel_RT = `retention time`/max(`retention time`)))



# add dataset summary overview list-tibble --------------------------------
results <- results %>%
  as_data_frame() %>%
  mutate(ID_ovw = map(run_smry, ~ filter(., sample_id == "all")))

# add protein info from summary table --------------------------------------
results <- results %>%
  mutate(
    ID_ovw = map2(
      ID_ovw,
      results$run_smry %>%
        map( ~ filter(., sample_id == "target") %>% .[, 2]) ,
      ~ add_column(.x, "protein_ID" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$run_smry %>%
        map( ~ filter(., sample_id == "decoy") %>% .[, 2]) ,
      ~ add_column(.x, "decoy_ID" = .y[[1]])
    )
  )
# add info from protein table and peptide table  ----------------------------
results <- results %>%
  mutate(
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ count(.)) ,
      ~ add_column(.x, "all_ID_table" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ filter(., .$qssm > 0) %>% count(.)) ,
      ~ add_column(.x, "quantified" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ filter(., .$qssm > 1) %>% count(.)) ,
      ~ add_column(.x, "quantified_min2" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ filter(., str_detect(protein_id, "###")) %>% count()),
      ~ add_column(.x, "decoy_ID_table" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot  %>%  map( ~ select(., contains("signal_sum")) %>% rowSums(na.rm = TRUE) %>% sum),
      ~ add_column(.x, "total_TMT_intSum" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ top_n(., 10, ssm) %>% select(gene_name) %>% str_c()),
      ~ add_column(.x, "TOP10_SPC" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map( ~ top_n(., 10, top3) %>% select(gene_name) %>% str_c(collapes = "")),
      ~ add_column(.x, "TOP10_intensity" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$prot %>% map(
        ~ top_n(., 10, total_signal) %>% select(gene_name) %>% str_c(collapes = "")
      ),
      ~ add_column(.x, "TOP10_TMT_intensity" = .y[[1]])
    ),
    ID_ovw = map2(
      ID_ovw,
      results$pep  %>%  map(~ max(.$`retention time`) %>% round()),
      ~ add_column(.x, "gradientLength" = .y[[1]])
    )
  )

results %>% colnames()
results$run
results$short_name


#### overview table
dataset_overview <<-
  results %>% select(run, ID_ovw) %>% unnest() %>%
  select(
    run,
    protein_ID,
    quantified,
    quantified_min2,
    all_ID_table,
    decoy_ID,
    decoy_ID_table,
    TOP10_intensity,
    TOP10_SPC,
    TOP10_TMT_intensity,
    acquired_spectra,
    mascot_matched_spectra,
    spectra_in_qc_proteins,
    quantified_spectra,
    mean_precursor_ion_accuracy,
    sd_precursor_ion_accuracy,
    mean_reporter_ion_accuracy,
    sd_reporter_ion_accuracy,
    total_TMT_intSum,
    gradientLength,
    run
  )

}
#remove optionally

#rm(clean)
# rm(results)
# rm(results_files)
