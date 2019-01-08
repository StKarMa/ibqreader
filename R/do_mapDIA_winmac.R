#'system.file("extdata/win", "mapDIA.exe", package = "ibqreader")
#' this does everything to do a mapDIA_analysis
#' it should recognize the operating system automatically
#' combine gene_sequence option to true will treat peptides as proteins
#' per definition the binary comparisons are made with the alphabetically last
#' sample name so you can influence it via the name.
#' system.file("extdata/win", "mapDIA.exe", package = "ibqreader")
#' @export
ibq_mapDIA <- function(combine_gene_sequence){

#md_tib$peptide <- md_tib$peptide %>% modify_depth(.,1, ~ibq_split_pepident(.))


  if (combine_gene_sequence == TRUE){
    md_tib$peptide <- md_tib$peptide %>%
      modify_depth(.,1, ~ibq_combine_geneandsequence(.))
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



# dealing with people that want freedom in choosing their min obse --------

  if (md_tib$param_0 %>% str_detect(.,  "MIN_OBS =")) {
    get_number_for_min_obs <-
      md_tib$param_0 %>% str_extract(.,  "MIN_OBS = \\d+") %>%
      str_remove(., "MIN_OBS = ") %>% as.double() %>% .[[1]]

    md_tib$param_0 <-
      md_tib$param_0 %>% str_remove(.,  "MIN_OBS = \\d+")
  }

  ###adding the missing data dependant parameters from peptide table
  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(
      .x,
      "LABELS = ",
      .y %>% names %>% .[-c(1:3)] %>%  str_extract("^.*._") %>%
        str_replace("_$", "") %>% unique %>% str_c(collapse = " "),
      "\n"
    )
  )

  # setting up the sample name dependant paramteres for param file ----------


  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(.x,
             "SIZE = ",
             group_samples(.y) %>% t %>% str_c(collapse = " "), "\n")
  )

  md_tib$param_0 <- map2(
    .x = md_tib$param_0,
    .y = md_tib$peptide,
    ~ paste0(
      .x,
      "MIN_OBS = ",
      group_samples(.y) %>% transmute(n / n) %>%  as_vector() %>%
        str_c(collapse = " ") %>%
        str_replace_all(., "1", paste0(get_number_for_min_obs)),
      "\n"
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
                               paste0("
                                      - 1
                                      0 -")
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
                                             - 0 0 0 0 0 0 0 0 1
                                             0 - 0 0 0 0 0 0 0 1
                                             0 0 - 0 0 0 0 0 0 1
                                             0 0 0 - 0 0 0 0 0 1
                                             0 0 0 0 - 0 0 0 0 1
                                             0 0 0 0 0 - 0 0 0 1
                                             0 0 0 0 0 0 - 0 0 1
                                             0 0 0 0 0 0 0 - 0 1
                                             0 0 0 0 0 0 0 0 - 1
                                             0 0 0 0 0 0 0 0 0 -"
                                           )
                                           )
                                           )
                                           )
                                           )
                                           )
                                           )
                                           )
                                           ))


  ### now generate a mapDIA folder and subfolders for the different runs
  dir.create("mapDIA")
  md_tib$run %>% map(., ~dir.create(paste0("mapDIA/", .)))


  walk2(md_tib$run,
        md_tib$peptide,
        ~write_tsv(.y, paste0("mapDIA/", .x, "/", "data.txt")))

  map2(md_tib$run,
       md_tib$param_0,
       ~write(.y , paste0("mapDIA/", .x, "/", "input")))



# doing this step dependent of operating system ---------------------------


  # function adopted from
  # https://www.r-bloggers.com/identifying-the-os-from-r/
  get_os <- function() {
    sysinf <- Sys.info()
    if (!is.null(sysinf)) {
      os <- sysinf['sysname']
      if (os == 'Darwin')
        os <- "osx"
    } else {
      ## mystery machine
      os <- .Platform$OS.type
      if (grepl("^darwin", R.version$os))
        os <- "osx"
      if (grepl("linux-gnu", R.version$os))
        os <- "linux"
    }
    tolower(os)
  }


  if (get_os() == "windows") {
    # finding the executable and running mapDIA  ------------------------------

    mapdiaexe <-
      system.file("extdata/win", "mapDIA.exe", package = "ibqreader") %>%
      str_replace_all("/", "\\\\")

    call.mapDIA <- map(md_tib$run,
                       ~ (shell((
                         paste0(
                           "cd /D ",
                           normalizePath(getwd()),
                           "\\mapDIA\\",
                           .x,
                           " ",
                           "& ",
                           mapdiaexe,
                           "  input"
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
                                               ),
                                               col_types = "ciiccddddd"
                                             )))

    md_tib <- md_tib %>% mutate(mapDIA_log2 = map(md_tib$run,
                                                  ~ read_tsv(
                                                    paste0(normalizePath(getwd()),
                                                           "\\mapDIA\\",
                                                           .x,
                                                           "\\log2_data.txt")
                                                  )))






  } else {

    ### this would be teh maccode now
    mapdiaexe <-
      system.file("extdata/mac", "mapDIA", package = "ibqreader")

    call.mapDIA <- map(md_tib$run,
                       ~ (system((
                         paste0(
                           "cd  ",
                           normalizePath(getwd()),
                           "/mapDIA/",
                           .x,
                           " ",
                           "&& ",
                           mapdiaexe,
                           "  input"
                         )
                       ),
                       wait = TRUE)))



    # reading mapDIA result files ---------------------------------------------
    md_tib <- md_tib %>% mutate(mapDIA = map(md_tib$run,
                                             ~ read_tsv(
                                               paste0(normalizePath(getwd()),
                                                      "/mapDIA/",
                                                      .x,
                                                      "/analysis_output.txt")
                                             )))

    md_tib <- md_tib %>% mutate(mapDIA_log2 = map(md_tib$run,
                                                  ~ read_tsv(
                                                    paste0(normalizePath(getwd()),
                                                           "/mapDIA/",
                                                           .x,
                                                           "/log2_data.txt")
                                                  )))

  }

# writeing a table for an alternative view in prohits viz
  map2(md_tib$run,
       md_tib$mapDIA,
       ~write_tsv(.y %>% mutate(log_oddsDE = 1-FDR)
                  , paste0("mapDIA/", .x, "/", "ONEminusFDR.txt")))


md_tib <<-md_tib
}



