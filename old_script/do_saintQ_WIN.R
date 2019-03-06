#' system.file("extdata/win", "saintq.exe", package = "ibqreader")
#' this does everything to do a mapDIA_analysis
#' @export
ibq_saintQprot_WIN <- function(combine_gene_sequence){

## so when we have external controls then combine them with the experiment data frames
  if (nrow(control_tib_p) > 0) {
    saint_tib_p$protein <-
      saint_tib_p$protein %>%
      map(.,
          ~ left_join(
            .,
            control_tib_p["protein"] %>% unnest %>%
              gather(
                .,
                key,
                value,
                na.rm = T,
                control_tib_p["protein"] %>%
                  unnest  %>% select_if(is.numeric) %>% colnames()
              ) %>%
              spread(key, value, fill = 0)
          ))
  }

  saint_tib_p$protein <-
  saint_tib_p$protein %>%
    modify_depth(., 1, ~.x %>% replace(., is.na(.), 0))



  # preparing header lines for SaintQ ---------------------------------------



  saint_tib_p$do_saint_form <-
    saint_tib_p$protein %>%
    modify_depth(., 1,
                 ~ .x %>% select(-gene_name) %>% colnames %>% as_tibble() %>% set_colnames(value = "bait") %>% # gather(., channel, bait) # %>% ### transform table
                   separate(bait, c("a", "b", "c"), sep = "_", remove = FALSE) %>% ### get naked identifiers a_b_c -> a | b | no c
                   unite("exp", a, b, remove = FALSE) %>%  ## combine them
                   mutate(iscontrol = ### make new line for saint header indicating if it is control or test sample
                            .$a %>% str_detect(., "CONT") %>%
                            if_else(., "C", "T")) %>%
                   select(iscontrol, exp, bait) %>% #then order it according
                   group_by(iscontrol) %>%
                   arrange(bait) %>% t %>% as_tibble()  %>%
                   set_colnames(., slice(., 3)) %>%  # depending if we want to match channels or baitnames 4 = baitnames
                   #slice(2:4) %>%  # keeping only what we need and that should be good to match with
                   add_column(.,
                              gene_name = c("", "", "gene_name"),
                              .before = 1) # and bringiing it into the right format for SAINTq
    )

  saint_tib_p$do_saint <-
    map2(
      saint_tib_p$do_saint_form,
      saint_tib_p$protein,
      ~bind_rows(.x,.y %>% mutate_all(., ~as.character(.))) ### needed to change format to character in order to allow merging
    )


  # preparing and writing the data ---------------------------------------------

  dir.create("SAINTp")
  saint_tib_p$run %>% map(., ~dir.create(paste0("SAINTp/", .)))

  walk2(saint_tib_p$run,
        saint_tib_p$do_saint,
        ~write_tsv(.y, paste0("SAINTp/", .x, "/", "data.txt"), col_names = FALSE))





  walk(saint_tib_p$run,
       ~write(parameter.saintq.protein , paste0("SAINTp/", .x, "/", "param_prot_level_win")))



  saintexe <- system.file("extdata/win", "saintq.exe", package = "ibqreader") %>% str_replace_all("/", "\\\\")

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

  saint_tib_p <- saint_tib_p %>% mutate(saint_scorelist = map(saint_tib_p$run,
                                                              ~ read_tsv(
                                                                paste0(
                                                                  normalizePath(getwd()),
                                                                  "\\SAINTp\\",
                                                                  .x,
                                                                  "\\scores_list__data.txt__.tsv"
                                                                ), col_types = "ccidd"
                                                              )
  ))




  # combine average channel intensity with saint results --------------------


  ####_ function

  average_channels <- function(a,b,...)
  {
    map_dfc(
      a, ~b %>%
        select(contains(paste0(.x)))%>%
        transmute(rowSums(.)) / length(a)
    ) %>% set_colnames(paste0(a)) %>% bind_cols(b[,1],.)
  }

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
