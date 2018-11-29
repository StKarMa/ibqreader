#' if we want to apply general correction factors, specified int the pooling.csv here we are
#' @export
ibq_apply_coorection_factors <- function(){

  data_baret <-
    data_bare %>%
    mutate("korrection" =
             map(data_bare$bait_channel, ~.x %>%
                   slice(2) %>% unlist %>% as.vector %>% as.double)
    )

  map2(data_baret$peptide,
       data_baret$korrection,
       ~.x %>% .[, 2:ncol(.x)] %>% colnames(.) %>% bind_cols("sample" = ., "korrection" = .y) %>% t
  )



## this is still not too cool. i should be able to do it via select_if...
  ##then i alos do not need to care about if i have the peptide info separated
  # apply to peptide
  # data_baret$peptide <-
  #   map2(data_baret$korrection,
  #        data_baret$peptide,
  #        ~map2_dfc(.y[,4:ncol(.y)], .x, `*`) %>%
  #          bind_cols(.y[,1:3], .)   )

  data_baret$peptide <-
    map2(data_baret$korrection,
         data_baret$peptide,
         ~map2_dfc(.y[,2:ncol(.y)], .x, `*`) %>%
           bind_cols(.y[,1], .)   )


  # apply to protein level
  data_baret$protein <-
    map2(data_baret$korrection,
         data_baret$protein,
         ~map2_dfc(.y[,2:ncol(.y)], .x, `*`) %>%
           bind_cols(.y[,1], .)   )

data_bare <<- data_baret

}
