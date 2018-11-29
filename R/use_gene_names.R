#' give the gene names to your tibbles
#' it also deals with duplicated stuff ... it will aggregate (sum for now) the single duplictes
#' @export
ibq_use_gene_names <- function() {

  data_bare$protein <-
    data_bare$protein %>%
    map(
      .,
      ~ left_join(., unique_protein_gene_ID, by = c("gene_name" = "protein_id")) %>%  #joining with the lookup
        ungroup %>%
        select(-"description", -"gene_name") %>%  # getting rid of unused
        select(ncol(.), 1:(ncol(.) - 1))  %>% # getting original order back
        rename(., gene_name = gene_name.y) %>%
        group_by(gene_name) %>% summarize_if(., is.numeric, sum) %>% ungroup
    )


  data_bare$peptide <-
    data_bare$peptide %>%
    map(
      .,
      ~ left_join(
        .x %>%
          separate(
            .,
            gene_name,
            c("protein_id", "sequence", "msms_id"),
            remove = TRUE,
            sep = "_-_"
          ),
        unique_protein_gene_ID
      ) %>%   #joining with the lookup
        ungroup %>%
        select(-"description", -"protein_id") %>% select(ncol(.), 1:(ncol(.) - 1)) %>%
        group_by(gene_name, sequence, msms_id) %>% summarize_if(., is.numeric, sum) %>% ungroup
    )


  data_bare <<- data_bare
}
