#' prepare data for saintQ SC
#'
#' @param saint_tib_SC

#'
#' @return
#' @export
#'
ibq_saintQ_SC <-
  function() {

saint_tib_SC$do_saint_form <-
saint_tib_SC$SC %>%

modify_depth(., 1,
             # from the protein tibbles we build a tibble with all "bait" names
             # and put them into one tibble where we can sepparate the 3-part name
             # and define the 2 first blocks of the name as exp = experiment
             ~ .x %>%
               colnames %>%
               enframe(value = "bait", name = NULL) %>%
               mutate(exp = bait) %>%
               mutate(iscontrol = .$bait %>%
                        str_detect(., "Controls") %>%
                        if_else(., "C", "T")) %>%
               arrange(iscontrol, bait) %>%
               t %>%
               as_tibble %>% slice(3:1) %>%
               set_colnames(., slice(., 2)) %>%
               select(-gene_name) %>%
               add_column(.,
                          gene_name = c("", "", "gene_name"),
                          .before = 1)
    )
saint_tib_SC$do_saint <-
    map2(
      saint_tib_SC$do_saint_form,
      saint_tib_SC$SC,
      ~bind_rows(.x,.y %>%
                   mutate_all(., ~as.character(.)))
      ### change format in order to allow merging
    )


dir.create("SAINTpSC")

saint_tib_SC$run %>%
   map(., ~dir.create(paste0("SAINTpSC/", .)))

walk2(saint_tib_SC$run,
        saint_tib_SC$do_saint,
        ~write_tsv(.y,
                   paste0("SAINTpSC/", .x, "/", "data.txt"), col_names = FALSE))

walk(saint_tib_SC$run,
       ~write(parameter.saintq.SC ,
              paste0("SAINTpSC/", .x, "/", "param_prot_level_win")))

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

##############################################################################
if (get_os() == "windows") {
  saint_tib_SC <- call_saintQP_SC_windows()
} else {
  saint_tib_SC <- call_saintQP_SC_mac()
}

# combine average channel intensity with saint results --------------------

# for a better description of this formula consult the saintq
#protein function of this package

average_channels <- function(a,b,...)
{
  map_dfc(
    a, ~b %>%
      select(contains(paste0(.x)))%>%
      transmute(rowMeans(.))
  ) %>% set_colnames(paste0(a)) %>% bind_cols(b[,1],.)
}

saint_tib_SC$prohitsviz <-
  map2(
    .x = saint_tib_SC$saint_scorelist %>%
      modify_depth(., 1,
                   ~ select(., Bait) %>%
                     distinct %>% as_vector) ,

    .y = saint_tib_SC$peptide,

    ~ average_channels(.x, .y)
  ) %>%

  modify_depth(.,
               1,

               ~ gather(., Bait, AvgSpec,-gene_name) %>%
                 unite(., combo, c(gene_name, Bait), sep = "-_-_"))




saint_tib_SC$prohitsviz <- map2(
  saint_tib_SC$prohitsviz,
  saint_tib_SC$saint_scorelist %>%
    modify_depth(., 1, ~ unite(., combo, c(Prey, Bait), sep = "-_-_")),
  left_join
) %>%
  modify_depth(., 1, ~ separate(., combo, c("Prey", "Bait"), sep = "-_-_"))


walk2(
  saint_tib_SC$run,
  saint_tib_SC$prohitsviz,
  ~write_tsv(.y %>%
               select(Bait, Prey, AvgP, BFDR, AvgSpec),
             paste0("SAINTpSC/", .x, "/", "for_prohits_viz.txt"),
             col_names = TRUE))


saint_tib_SC <<- saint_tib_SC





}




#' runs SaintQ on spectral counts for windows
#'
#' @return
#' @export
#'
call_saintQP_SC_windows <- function() {

    saintexe <-
      system.file("extdata/win", "saintq.exe", package = "ibqreader") %>%
      str_replace_all("/", "\\\\")

    call.SAINTq <- map(saint_tib_SC$run,
                       ~ (shell((
                         paste0(
                           "cd /D ",
                           normalizePath(getwd()),
                           "\\SAINTpSC\\",.x, " ",
                           "& ", saintexe, " param_prot_level_win"
                         )
                       ),
                       wait = TRUE)))


    # reading saint result files ----------------------------------------------

    saint_tib_SC$saint_scorelist <-

               map(saint_tib_SC$run,
                   ~ read_tsv(
                     paste0(
                       normalizePath(getwd()),
                       "\\SAINTpSC\\",
                       .x,
                       "\\scores_list__data.txt__.tsv"
                     ))
               )
  }

#' runs SaintQ on spectral counts for windows
#'
#' @return
#' @export
#'
call_saintQP_SC_mac <- function() {



  saintexe <-
    system.file("extdata/mac", "saintq", package = "ibqreader") %>%
    str_replace_all("/", "\\\\")

  call.SAINTq <- map(saint_tib_SC$run,
                     ~ (system((
                       paste0(
                         "cd  ",
                         normalizePath(getwd()),
                         "/SAINTpSC/",.x, " ",
                         "&& ", saintexe, " param_prot_level_win"
                       )
                     ),
                     wait = TRUE)))


  # reading saint result files ----------------------------------------------

  saint_tib_SC <- saint_tib_SC %>%
    mutate(saint_scorelist = map(saint_tib_SC$run,
                                 ~ read_tsv(paste0(
                                   normalizePath(getwd()),
                                   "/SAINTpSC/",
                                   .x,
                                   "/scores_list__data.txt__.tsv"
                                 ),
                                 col_types = "cciicdd")
    ))
}


