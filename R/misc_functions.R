#' this function names and selects TMT channels that are  specified in the pooling.csv
#' @import magrittr
#' @export

ibq_name_things <- function(bait_channel, data) {

  selecter <- bait_channel %>% colnames()
  identifier <- data %>% select_if(., ~ !is_numeric(.))
  data <- data %>% select_if(., ~ is_numeric(.))

# this is needed as the channels for peptides have different names.
# now we can use the same function for
  colnames(data) <-
    colnames(data) %>% str_replace_all(., "sig_", "signal_sum_")

    data <-
    data %>% select_if(., ~ is_numeric(.)) %>% select(selecter) %>%
    set_colnames(., c(paste0(bait_channel[1,])))

  bind_cols(identifier, data)
}

#' this function calculates the relative TMT intensities row so that the sum is one.
#' @export
ibq_relative_intensities <- function(bare_peptides_or_protein_data){
  identifier <- bare_peptides_or_protein_data %>% select_if(., ~!is_numeric(.))
  data <- bare_peptides_or_protein_data %>% select_if(., ~is_numeric(.))
  data <- as.matrix(data) %>% prop.table(., 1)
  cbind(identifier, data) %>% as_tibble(.)
}

#' get form fake combined unique peptide column to protein|sequence|msms_if
#' @export
ibq_split_pepident <- function(peptide_tibble){

  peptide_tibble %>%   separate(., gene_name,
                      c("gene_name", "sequence", "msms_id"),
                      sep = "_-_")
}
#' get
#' @export
ibq_combine_geneandsequence <- function(peptide_tibble){

  peptide_tibble %>%   mutate(gene_name = paste(gene_name, sequence, sep = "_"))
}


#' plotting a boxpot
#' @export
boxplot_10set_list_tibble <- function(x, y_intercept = 0.1){
    x %>%
    modify_depth(., 1,
                 ~ .x %>%
                   gather(., sample, value, -gene_name, -sequence, -msms_id)) %>%
    map(
      ~ .x %>% ggplot(aes(sample, (value))) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_size(range=c(0.5,0.5)) +
        theme(legend.position="none") +
        geom_hline(yintercept = y_intercept, color = "red")
    )
}
