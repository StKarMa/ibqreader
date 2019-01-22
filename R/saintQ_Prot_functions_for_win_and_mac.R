#' this function actually calls saintQprotein and takes care of windows pathnames
#' @export
call_saintQp_windows  <-  function() {

### here we tell the script where to find the executable
saintexe <-
  system.file("extdata/win", "saintq.exe", package = "ibqreader") %>%
  str_replace_all("/", "\\\\")

call.SAINTq <- map(saint_tib_p$run,
                   ~ (shell((
                     paste0(
                       "cd /D ",
                       normalizePath(getwd()),
                       "\\SAINTp\\",.x, " ",
                       "& ", saintexe,  "  param_prot_level_win"
                     )
                   ),
                   wait = TRUE)))


# reading saint result files ----------------------------------------------

saint_tib_p <-
  saint_tib_p %>%
  mutate(saint_scorelist =
           map(saint_tib_p$run,
               ~ read_tsv(
                 paste0(
                   normalizePath(getwd()),
                   "\\SAINTp\\",
                   .x,
                   "\\scores_list__data.txt__.tsv"
                 ),
                 col_types = "ccidd"
               )
           ))
}


#' this function actually calls saintQprotein and takes care of mac pathnames
#' @export
call_saintQp_mac  <-  function() {
saintexe <- system.file("extdata/mac", "saintq", package = "ibqreader")

call.SAINTq <- map(saint_tib_p$run,
                   ~ (system((
                     paste0(
                       "cd ",
                       normalizePath(getwd()),
                       "/SAINTp/",.x, " ",
                       "&& ", saintexe,  "  param_prot_level_win"
                     )
                   ),
                   wait = TRUE)))


# reading saint result files ----------------------------------------------

saint_tib_p <-
  saint_tib_p %>%
  mutate(
    saint_scorelist = map(saint_tib_p$run,
                          ~ read_tsv(
                            paste0(
                              normalizePath(getwd()),
                              "/SAINTp/",
                              .x,
                              "/scores_list__data.txt__.tsv"
                            ),
                            col_types = "ccidd")
))
}


