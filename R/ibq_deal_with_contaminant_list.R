#' to filter out common contaminants via gene_name or protein id
#' a short list is delivered with the package and generated in the working folder
#' it will not be generated if a list already exists
#' so you can use your own
#' @export
ibq_filter_contam <- function(){

  # loading example file form package

cont <-
  system.file("extdata", "remove_contaminants.txt", package = "ibqreader") %>%
  read_tsv

  # generating example file in working directory if not existent
ifelse(file.exists("remove_the_contaminants.txt"), "yes", write_tsv(cont, "remove_the_contaminants.txt"))

 # laoding file form working dir and making it no matter if the list contains gene names or protein id
# veryveryelegant

contaminants <<-
  read_tsv("remove_the_contaminants.txt") %>%
  left_join(., unique_protein_gene_ID[-2],
            by = c("Contaminants" = "gene_name")) %>%
  left_join(., unique_protein_gene_ID[-2],
            by = c("Contaminants" = "protein_id")) %>% gather() %>% .[-1] %>%
  filter(!value == "NA") %>% unlist





data_bare$protein <- data_bare$protein %>%
  map(., ~ .x %>% filter(!gene_name %in% contaminants ))


### maybe just leave the rowid stuff out biut for now....
data_bare$peptide <- data_bare$peptide %>%
  map(., ~ .x %>%
        separate(
          .,
          gene_name,
          c("protein_id", "sequence", "msms_id"),
          remove = FALSE,
          sep = "_-_"
        ) %>%
        filter(!protein_id %in% contaminants ) %>%
        select( -"protein_id", -"sequence", -"msms_id")
      )

data_bare <<-data_bare

}

