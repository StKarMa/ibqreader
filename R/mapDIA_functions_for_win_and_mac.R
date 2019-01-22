#' this function actually calls mapDIA and takes care of windows pathnames
#' @export
call_mapDIA_windows <- function() {


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



}



#' this function actually calls mapDIA and takes care of windows pathnames
#' @export
call_mapDIA_mac <- function() {

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
  md_tib <<- md_tib %>% mutate(mapDIA = map(md_tib$run,
                                           ~ read_tsv(
                                             paste0(normalizePath(getwd()),
                                                    "/mapDIA/",
                                                    .x,
                                                    "/analysis_output.txt")
                                           )))

  md_tib <<- md_tib %>% mutate(mapDIA_log2 = map(md_tib$run,
                                                ~ read_tsv(
                                                  paste0(normalizePath(getwd()),
                                                         "/mapDIA/",
                                                         .x,
                                                         "/log2_data.txt")
                                                )))

  md_tib <<- md_tib

}
