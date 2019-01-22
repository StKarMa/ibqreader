#' this function keeps the tibble structure and prepares the data for following steps
#' NA values are replaced by 0
#' it writes a file with all gene, protein, description information into the working directory
#' min_number_of_quantified_peps lets you choose how many peptides need to be quantified per run
#' to accept a protein as identidied and quantified
#' @import magrittr
#' @export
#'
#'

ibq_clean_data_TMTX <- function(results, min_number_of_quantified_peps, ... ){

# replacing missing values with 0 -----------------------------------------

results$prot <-
results$prot %>%
  modify_depth(., 1,
  ~ mutate_at(., vars(contains("signal_sum")), ~ replace(., is.na(.), 0)))

results$pep <-
  results$pep %>%
  modify_depth(., 1,
  ~ mutate_at(., vars(contains("sig_")), ~ replace(., is.na(.), 0)))

unique_protein_gene_ID <<-
  results$prot %>%
  modify_depth(., 1, ~ select(., gene_name, description, protein_id)) %>%
  bind_rows() %>% unique()

unique_protein_gene_ID %>% write_csv("Protein_Gene_Description.csv")

summarise_all(unique_protein_gene_ID, dplyr::n_distinct)

# "clean" as in no reverse hits and only peptides rank 1 and in quantification-----------------------------

clean <-
  results %>% select(run, short_name, bait_channel) %>% as_tibble

clean <- clean %>% mutate(
  protein = results$prot %>%
    map(~ filter(.,!str_detect(protein_id, "###"))) %>%
    map(~ filter(
      ., .$qupm >= min_number_of_quantified_peps
    )),

  peptide = results$pep %>%
    map(~ filter(.,!str_detect(protein_id, "###"))) %>%
    map(~ filter(., .$in_quantification_of_protein == 1))  %>%
    map(~ filter(., .$rank == 1))
)

clean <<- clean

# the bare data has only the TMT intensity information for the gen --------

data_bare <-
  clean %>%
  select(run, bait_channel) %>%
  mutate(

    protein = clean[["protein"]] %>%
      map(.,
          . %>% mutate(gene_name = protein_id) %>%
            select(., gene_name, contains("signal_sum")) ),



    peptide = clean[["peptide"]] %>%
      map(
        .,
        . %>% mutate(gene_name = protein_id) %>%
          select(., protein_id, sequence,  msms_id, source_file, contains("sig_")) %>%
          rowid_to_column(var = "rowID") %>%
          ##to accomodate for fractions in msms_id
          mutate(source_file =
                   str_trunc(source_file, 7, "left", ellipsis = "_") %>% str_remove(., ".raw")) %>%
          mutate(msms_id = paste0(msms_id, source_file)) %>% select(-source_file) %>%
         # unite(., "unique_pepID", rowID, protein_id, sequence, msms_id, sep = "_-_")
        unite(., "gene_name", protein_id, sequence, msms_id, sep = "_-_")
      )
  )

### quickfix to make msms_id even more unique
###(odd case here the same peptide is identified in 2 runs in the exact same MSMS scan number)
data_bare$peptide <- map2(data_bare$run, data_bare$peptide,
                          ~ .y %>% mutate(gene_name = paste(gene_name, .x, sep = "_"))
)

data_bare$protein <- map2(data_bare$bait_channel,
                          data_bare$protein,
                          ~ ibq_name_things(.x, .y))

data_bare$peptide <- map2(data_bare$bait_channel,
                          data_bare$peptide,
                          ~ ibq_name_things(.x, .y))



data_bare <<- data_bare

}


#' same as the other but attaches modification information
#' @import magrittr
#' @export
#'
#'

ibq_clean_data_keep_mods_TMTX <- function(results, min_number_of_quantified_peps, ... ){

  # replacing missing values with 0 -----------------------------------------

  results$prot <-
    results$prot %>% modify_depth(., 1,
                                  ~ mutate_at(., vars(contains("signal_sum")), ~ replace(., is.na(.), 0)))

  results$pep <- results$pep %>% modify_depth(
    .,
    1,
    ~ .x %>% mutate(
      sequence = modifications %>%
        str_replace_all(., "TMT6plex N-term;*", "") %>%
        str_replace_all(., "TMT6plex K..*;", "") %>%
        str_replace_all(., "Carbamidomethyl ", "Cam") %>%
        str_replace_all(., "Deamidated ", "Dea") %>%
        str_replace_all(., "Oxidation ", "Ox") %>%
        str_replace_all(., "Phospho ", "Ph") %>%
        str_replace_all(., " ", "") %>%
        paste0(sequence, ";", .)

    )
  )




  results$pep <-
    results$pep %>% modify_depth(., 1,
                                 ~ mutate_at(., vars(contains("sig_")), ~ replace(., is.na(.), 0)))

  unique_protein_gene_ID <<-
    results$prot %>%
    modify_depth(., 1, ~ select(., gene_name, description, protein_id)) %>%
    bind_rows() %>% unique()

  unique_protein_gene_ID %>% write_csv("Protein_Gene_Description.csv")

  summarise_all(unique_protein_gene_ID, dplyr::n_distinct)

  # "clean" as in no reverse hits and only peptides rank 1 and in quantification-----------------------------

  clean <-
    results %>% select(run, short_name, bait_channel) %>% as_tibble

  clean <- clean %>% mutate(
    protein = results$prot %>%
      map(~ filter(.,!str_detect(protein_id, "###"))) %>%
      map(~ filter(
        ., .$qupm >= min_number_of_quantified_peps
      )),

    peptide = results$pep %>%
      map(~ filter(.,!str_detect(protein_id, "###"))) %>%
      map(~ filter(., .$in_quantification_of_protein == 1))  %>%
      map(~ filter(., .$rank == 1))
  )

  clean <<- clean

  # the bare data has only the TMT intensity information for the gen --------

  data_bare <-
    clean %>%
    select(run, bait_channel) %>%
    mutate(

      protein = clean[["protein"]] %>%
        map(.,
            . %>% mutate(gene_name = protein_id) %>%
              select(., gene_name, contains("signal_sum"))),

      peptide = clean[["peptide"]] %>%
        map(
          .,
          . %>% mutate(gene_name = protein_id) %>%
            select(., protein_id, sequence,  msms_id, source_file, contains("sig_")) %>%
            rowid_to_column(var = "rowID") %>%
            ##to accomodate for fractions in msms_id
            mutate(source_file =
                     str_trunc(source_file, 7, "left", ellipsis = "_") %>% str_remove(., ".raw")) %>%
            mutate(msms_id = paste0(msms_id, source_file)) %>% select(-source_file) %>%
            # unite(., "unique_pepID", rowID, protein_id, sequence, msms_id, sep = "_-_")
            unite(., "gene_name", protein_id, sequence, msms_id, sep = "_-_")
        )
    )

  ### quickfix to make msms_id even more unique
  ###(odd case here the same peptide is identified in 2 runs in the exact same MSMS scan number)
  data_bare$peptide <- map2(data_bare$run, data_bare$peptide,
                            ~ .y %>% mutate(gene_name = paste(gene_name, .x, sep = "_"))
  )

  data_bare$protein <- map2(data_bare$bait_channel,
                            data_bare$protein,
                            ~ ibq_name_things(.x, .y))

  data_bare$peptide <- map2(data_bare$bait_channel,
                            data_bare$peptide,
                            ~ ibq_name_things(.x, .y))



  data_bare <<- data_bare

}
