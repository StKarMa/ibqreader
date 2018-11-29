#'system.file("extdata/win", "mapDIA.exe", package = "ibqreader")
#' this does everything to do a mapDIA_analysis
#' @export
ibq_mapDIA_WIN <- function(combine_gene_sequence){

  #md_tib$peptide <- md_tib$peptide %>% modify_depth(.,1, ~ibq_split_pepident(.))


  if (combine_gene_sequence == TRUE){
    md_tib$peptide <- md_tib$peptide %>% modify_depth(.,1, ~ibq_combine_geneandsequence(.))
  }



  md_tib$peptide <-
    md_tib$peptide %>% modify_depth(., 1,
                                    ~select(., -contains("CONT")) %>%
                                      select(1,2,3,
                                             .[ ,4:ncol(.)] %>% names %>% order+3
                                      ))

  # function that does the tricks for the parameter file of mapDIA
  group_samples <- function(x) {
    x %>%
      names %>% .[-c(1:3)] %>%                           #column names without
      str_extract("^.*._") %>% str_replace("_$", "") %>% #remove the post string
      as_tibble() %>% group_by(value) %>%                #grouping
      tally %>% select(n)
  }


  md_tib$param_0 <- paste(parameter_mapDIA)

  ###adding the missing data dependant parameters from peptide table
  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(
      .x,
      "LABELS = ",
      .y %>% names %>% .[-c(1:3)] %>%  str_extract("^.*._") %>%
        str_replace("_$", "") %>% unique %>% str_c(collapse = " "), "\n"
    )
  )

  # setting up the sample name dependant paramteres for param file ----------


  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(
      .x,
      "SIZE = ",
      group_samples(.y) %>% t %>% str_c(collapse = " "), "\n"
    )
  )

  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(
      .x,
      "MIN_OBS = ",
      group_samples(.y) %>% transmute(n/n) %>%  as_vector() %>% str_c(collapse = " ") , "\n"
    )
  )

  md_tib$param_0 <-   map2(.x = md_tib$param_0,
                           .y = md_tib$peptide,

                           ~ paste0(
                             .x,
                             "CONTRAST = " ,

                             ifelse(
                               group_samples(.y) %>%
                                 count() == 2,
                               paste0(
                                 "
                                 - 0
                                 1 -")
                               ,

                               ifelse(
                                 group_samples(.y) %>% count() == 3,
                                 "
                                 - 0 1
                                 0 - 1
                                 0 0 -
                                 "
                                 ,
                                 ifelse(
                                   group_samples(.y) %>% count() == 4,
                                   "
                                   - 0 0 1
                                   0 - 0 1
                                   0 0 - 1
                                   0 0 0 -"
                                   ,
                                   ifelse(
                                     group_samples(.y) %>% count() == 5,
                                     "
                                     - 0 0 0 1
                                     0 - 0 0 1
                                     0 0 - 0 1
                                     0 0 0 - 1
                                     0 0 0 0 -"
                                     ,
                                     ifelse(
                                       group_samples(.y) %>% count() == 6,
                                       "
                                       - 0 0 0 0 1
                                       0 - 0 0 0 1
                                       0 0 - 0 0 1
                                       0 0 0 - 0 1
                                       0 0 0 0 - 1
                                       0 0 0 0 0 -"
                                       ,
                                       ifelse(
                                         group_samples(.y) %>% count() == 7,
                                         "
                                         - 0 0 0 0 0 1
                                         0 - 0 0 0 0 1
                                         0 0 - 0 0 0 1
                                         0 0 0 - 0 0 1
                                         0 0 0 0 - 0 1
                                         0 0 0 0 0 - 1
                                         0 0 0 0 0 0 -"
                                         ,

                                         ifelse(
                                           group_samples(.y) %>% count() == 8,
                                           "
                                           - 0 0 0 0 0 0 1
                                           0 - 0 0 0 0 0 1
                                           0 0 - 0 0 0 0 1
                                           0 0 0 - 0 0 0 1
                                           0 0 0 0 - 0 0 1
                                           0 0 0 0 0 - 0 1
                                           0 0 0 0 0 0 - 1
                                           0 0 0 0 0 0 0 - "
                                           ,
                                           ifelse(
                                             group_samples(.y) %>% count() == 9,
                                             "
                                             - 0 0 0 0 0 0 0 1
                                             0 - 0 0 0 0 0 0 1
                                             0 0 - 0 0 0 0 0 1
                                             0 0 0 - 0 0 0 0 1
                                             0 0 0 0 - 0 0 0 1
                                             0 0 0 0 0 - 0 0 1
                                             0 0 0 0 0 0 - 0 1
                                             0 0 0 0 0 0 0 - 1
                                             0 0 0 0 0 0 0 0 -"
                                             ,
                                             "
                                             - 0 0 0 0 0 0 0 0 0
                                             0 - 0 0 0 0 0 0 0 0
                                             0 0 - 0 0 0 0 0 0 0
                                             0 0 0 - 0 0 0 0 0 0
                                             0 0 0 0 - 0 0 0 0 0
                                             1 0 0 0 0 - 0 0 0 0
                                             0 1 0 0 0 0 - 0 0 0
                                             0 0 1 0 0 0 0 - 0 0
                                             0 0 0 1 0 0 0 0 - 0
                                             0 0 0 0 1 0 0 0 0 -"
                                           )))
                                     )))))))








  ### now generate a mapDIA folder and subfolders for the different runs
  dir.create("mapDIA")
  md_tib$run %>% map(., ~dir.create(paste0("mapDIA/", .)))


  walk2(md_tib$run,
        md_tib$peptide,
        ~write_tsv(.y, paste0("mapDIA/", .x, "/", "data.txt")))

  map2(md_tib$run,
       md_tib$param_0,
       ~write(.y , paste0("mapDIA/", .x, "/", "input")))




  mapdiaexe <- system.file("extdata/win", "mapDIA.exe", package = "ibqreader") %>% str_replace_all("/", "\\\\")

  call.mapDIA <- map(md_tib$run,
                     ~ (shell((
                       paste0(
                         "cd /D ",
                         normalizePath(getwd()),
                         "\\mapDIA\\",.x, " ",
                         "& ", mapdiaexe,  "  input"
                       )
                     ),
                     wait = TRUE)))



  # reading mapDIA result files ---------------------------------------------

  md_tib <- md_tib %>% mutate(mapDIA = map(md_tib$run,
                                           ~ read_tsv(
                                             paste0(
                                               normalizePath(getwd()),
                                               "\\mapDIA\\",
                                               .x,
                                               "\\analysis_output.txt"
                                             ), col_types = "ciiccddddd"
                                           )))

  md_tib <- md_tib %>% mutate(mapDIA_log2 = map(md_tib$run,
                                                ~ read_tsv(
                                                  paste0(
                                                    normalizePath(getwd()),
                                                    "\\mapDIA\\",
                                                    .x,
                                                    "\\log2_data.txt"
                                                  )#, col_types = "cccdddddddd"
                                                )))


  map2(md_tib$run,
       md_tib$mapDIA,
       ~write_tsv(.y %>% mutate(log_oddsDE = 1-FDR)
                  , paste0("mapDIA/", .x, "/", "ONEminusFDR.txt")))


md_tib <<-md_tib
}



