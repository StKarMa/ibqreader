#' system.file("extdata/win", "saintq.exe", package = "ibqreader")
#' this does everything to do a sanitQ analyis on a windows machine
#' @export
#'
call_saintQf_windows <- function() {

# locating executable -----------------------------------------------------



saintexe <-
  system.file("extdata/win", "saintq.exe", package = "ibqreader") %>%
  str_replace_all("/", "\\\\")

call.SAINTq <- map(saint_tib_f$run,
                   ~ (shell((
                     paste0(
                       "cd /D ",
                       normalizePath(getwd()),
                       "\\SAINTf\\",.x, " ",
                       "& ", saintexe, " param_prot_level_win"
                     )
                   ),
                   wait = TRUE)))


# reading saint result files ----------------------------------------------

saint_tib_f <-
  saint_tib_f %>%
  mutate(saint_scorelist =
           map(saint_tib_f$run,
               ~ read_tsv(
                 paste0(
                   normalizePath(getwd()),
                   "\\SAINTf\\",
                   .x,
                   "\\scores_list__data.txt__.tsv"
                 ),
                 col_types = "cciicdd")
           ))

}



#' system.file("extdata/win", "saintq.exe", package = "ibqreader")
#' this does everything to do a sanitQ analyis on a windows machine
#' @export
#'
call_saintQf_mac <- function() {



saintexe <-
  system.file("extdata/mac", "saintq", package = "ibqreader") %>%
  str_replace_all("/", "\\\\")

call.SAINTq <- map(saint_tib_f$run,
                   ~ (system((
                     paste0(
                       "cd  ",
                       normalizePath(getwd()),
                       "/SAINTf/",.x, " ",
                       "&& ", saintexe, " param_prot_level_win"
                     )
                   ),
                   wait = TRUE)))


# reading saint result files ----------------------------------------------

saint_tib_f <- saint_tib_f %>%
  mutate(saint_scorelist = map(saint_tib_f$run,
                               ~ read_tsv(paste0(
                                 normalizePath(getwd()),
                                 "/SAINTf/",
                                 .x,
                                 "/scores_list__data.txt__.tsv"
                               ),
                               col_types = "cciicdd")
))
}
