#' system.file("extdata/win", "saintq.exe", package = "ibqreader")
#' this does everything to do a mapDIA_analysis
#' @export
ibq_saintQprot_new <- function(combine_gene_sequence = FALSE){

## so when we have external controls then combine them with the experiment data
# frames each protein data frame is left_joined with the control
if (nrow(control_tib_p) > 0) {

    saint_tib_p$protein <-

      map(saint_tib_p$protein,
          ~ left_join(.,
            control_tib_p["protein"] %>%
              unnest %>%
              gather(.,key,value,na.rm = T,

                control_tib_p["protein"] %>%
                  unnest  %>%
                  select_if(is.numeric) %>%
                  colnames()
              ) %>%

              spread(key, value, fill = 0)
          )
          )
  }

  saint_tib_p$protein <-
  saint_tib_p$protein %>%
    modify_depth(., 1, ~.x %>% replace(., is.na(.), 0))



# preparing header lines for SaintQ ---------------------------------------
# it is rather complicated

saint_tib_p$do_saint_form <-

    saint_tib_p$protein %>%

    modify_depth(., 1,
    # from the protein tibbles we build a tibble with all "bait" names
    # and put them into one tibble where we can sepparate the 3-part name
    # and define the 2 first blocks of the name as exp = experiment
    ~ .x %>%
         select(-gene_name) %>%
         colnames %>%
         as_tibble() %>%
         set_colnames(value = "bait") %>%
         separate(bait, c("a", "b", "c"), sep = "_", remove = FALSE) %>%
         unite("exp", a, b, remove = FALSE) %>%
    # with the help of the naming convention we can set samples as T (test)
    # or C (control) needed by SaintQ in a new line
         mutate(iscontrol = .$a %>%
                  str_detect(., "CONT") %>%
                            if_else(., "C", "T")
                ) %>%
                   select(iscontrol, exp, bait) %>%
    # now we order and group the headers and and flip/transform the tibble
                   group_by(iscontrol) %>%
                   arrange(bait) %>% t %>% as_tibble()  %>%
    # now we can decide if we want to use bait_names or channel names for the
    # further processing depending if we want to match channels or baitnames
    # 4 = baitnames
                   set_colnames(., slice(., 3)) %>%
                    # slice(2:4) %>%  # keeping only what we need and that should be good to match with

    # and last we need to put the final touches to make it a great saintQ header
      add_column(.,gene_name = c("", "", "gene_name"), .before = 1)

)

# now we can unite data and header by matching of rows by name
  saint_tib_p$do_saint <-
    map2(
      saint_tib_p$do_saint_form,
      saint_tib_p$protein,
      ~bind_rows(.x,.y %>%
        mutate_all(., ~as.character(.))) ### change format in order to allow merging
    )


# preparing and writing the data a folder per run   ----------------------------

dir.create("SAINTp")

saint_tib_p$run %>%
  map(., ~dir.create(paste0("SAINTp/", .)))

walk2(saint_tib_p$run,
      saint_tib_p$do_saint,
        ~write_tsv(.y,
                   paste0("SAINTp/", .x, "/", "data.txt"), col_names = FALSE))

walk(saint_tib_p$run,
       ~write(parameter.saintq.protein ,
              paste0("SAINTp/", .x, "/", "param_prot_level_win")))





# figuring out which os we are on and which function to call --------------

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
  saint_tib_p <- call_saintQp_windows()
} else {
  saint_tib_p <- call_saintQp_mac()
}




# combine average channel intensity with saint results --------------------

# a is a vector of baits (aka experiments/runs)
# b is the protein intensities table...
# so it, keeps only the matching columns an calculates the

average_channels <- function(a,b,...)
  {
    map_dfc(a,
          ~b %>%
            select(contains(paste0(.x)))%>%
            transmute(rowMeans(.))
    ) %>% set_colnames(paste0(a)) %>% bind_cols(b[,1],.)
  }


# here we actually calculate the average  ---------------------------------


saint_tib_p$prohitsviz <-
map2(

    .x = saint_tib_p$saint_scorelist %>%
          modify_depth(., 1,
                     ~ select(., Bait) %>%
                       distinct %>% as_vector) ,

    .y = saint_tib_p$protein,

    ~ average_channels(.x, .y)
    ) %>%

    modify_depth(.,
                 1,
    ~ gather(., Bait, AvgSpec,-gene_name) %>%
         unite(., combo, c(gene_name, Bait), sep = "-_-_"))




  saint_tib_p$prohitsviz <- map2(
    saint_tib_p$prohitsviz,
    saint_tib_p$saint_scorelist %>%
      modify_depth(., 1, ~ unite(., combo, c(Prey, Bait), sep = "-_-_")),
    left_join
  ) %>%
    modify_depth(., 1, ~ separate(., combo, c("Prey", "Bait"), sep = "-_-_"))


  walk2(
    saint_tib_p$run,
    saint_tib_p$prohitsviz,
    ~write_tsv(.y %>%
                 select(Bait, Prey, AvgP, BFDR, AvgSpec) %>%
                 mutate(
                   AvgSpec = 50 * 10 * AvgSpec/max(AvgSpec) ## 50 is the standard of prohits viz
                 ),
               paste0("SAINTp/", .x, "/", "for_prohits_viz.txt"),
               col_names = TRUE))


  saint_tib_p %>% select(run, saint_scorelist) %>% unnest %>% split(.$Bait) %>% map(~

                                                                                      ggplot(data = ., aes(BFDR, fill = Bait)) +
                                                                                      geom_histogram(bins = 20) +
                                                                                      facet_wrap(~ run)
  )

  saint_tib_p <<- saint_tib_p

}
