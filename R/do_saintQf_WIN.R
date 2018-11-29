#' system.file("extdata/win", "saintq.exe", package = "ibqreader")
#' this does everything to do a sanitQ analyis on a windows machine
#' @export
#'
ibq_saintQfrag_WIN <- function(combine_gene_sequence){




  ## in case we dont have external controls
  if(nrow(control_tib_f)>0){
    saint_tib_f$peptide <-
      map(saint_tib_f$peptide ,
          ~full_join(., control_tib_f["peptide"] %>%
                       unnest %>%
                       mutate(msms_id = paste0(msms_id,"c")) ))

  }



  ### replace NAs
  saint_tib_f$peptide <- saint_tib_f$peptide %>% modify_depth(., 1, ~.x %>% replace(., is.na(.), 0))

  ### bunch colnames


  # preparing header lines for SAINTf ---------------------------------------



  saint_tib_f$do_saint_form <-
    saint_tib_f$peptide %>%
    modify_depth(., 1,
                 ~ .x %>% select(-gene_name, -msms_id, -sequence) %>% colnames %>% as_tibble() %>% set_colnames(value = "bait") %>%
                   separate(bait, c("a", "b", "c"), sep = "_", remove = FALSE) %>% ### get naked identifiers a_b_c -> a | b | no c
                   unite("exp", a, b, remove = FALSE) %>%  ## combine them
                   mutate(iscontrol = ### make new line for saint header indicating if it is control or test sample
                            .$a %>% str_detect(., "CONT") %>%
                            if_else(., "C", "T")) %>%
                   select( iscontrol, exp, bait) %>% #then order it according
                   group_by(iscontrol) %>%
                   arrange(bait) %>% t %>% as_tibble()  %>%
                   set_colnames(., slice(., 3)) %>%  # depending if we want to match channels or baitnames 4 = baitnames
                   add_column(., msms_id = c("", "", "msms_id"),.before = 1) %>% # and bringiing it into the right format for SAINTf
                   add_column(., sequence = c("", "", "sequence"),.before = 1) %>%
                   add_column(., gene_name = c("", "", "gene_name"),.before = 1)
    )



  saint_tib_f$do_saint <-
    map2(
      saint_tib_f$do_saint_form,
      saint_tib_f$peptide,
      ~bind_rows(.x,.y %>% mutate_all(., ~as.character(.))) ### needed to change format to character in order to allow merging
    )


  # preparing and writing the data ---------------------------------------------

  dir.create("SAINTf")
  saint_tib_f$run %>% map(., ~dir.create(paste0("SAINTf/", .)))

  walk2(saint_tib_f$run,
        saint_tib_f$do_saint,
        ~write_tsv(.y, paste0("SAINTf/", .x, "/", "data.txt"), col_names = FALSE))

  # parameter.saintq.frag <-
  # "
  # normalize_control=true
  # input_filename=data.txt
  # input_level=fragment
  # protein_colname=gene_name
  # pep_colname=sequence
  # frag_colname=msms_id
  # compress_n_ctrl=100
  # compress_n_rep=100
  # min_n_pep=3
  # best_prop_pep=0.5
  # min_n_frag=3
  # best_prop_frag=0.5
  # "

  walk(saint_tib_f$run,
       ~write(parameter.saintq.frag , paste0("SAINTf/", .x, "/", "param_prot_level_win")))


  saintexe <- system.file("extdata/win", "saintq.exe", package = "ibqreader") %>% str_replace_all("/", "\\\\")
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

  saint_tib_f <- saint_tib_f %>% mutate(saint_scorelist = map(saint_tib_f$run,
                                                              ~ read_tsv(
                                                                paste0(
                                                                  normalizePath(getwd()),
                                                                  "\\SAINTf\\",
                                                                  .x,
                                                                  "\\scores_list__data.txt__.tsv"
                                                                ), col_types = "cciicdd"
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

  saint_tib_f$prohitsviz <-
    map2(
      .x = saint_tib_f$saint_scorelist %>%
        modify_depth(., 1,
                     ~ select(., Bait) %>%
                       distinct %>% as_vector) ,

      .y = saint_tib_f$peptide,

      ~ average_channels(.x, .y)
    ) %>%

    modify_depth(.,
                 1,

                 ~ gather(., Bait, AvgSpec,-gene_name) %>%
                   unite(., combo, c(gene_name, Bait), sep = "-_-_"))




  saint_tib_f$prohitsviz <- map2(
    saint_tib_f$prohitsviz,
    saint_tib_f$saint_scorelist %>%
      modify_depth(., 1, ~ unite(., combo, c(Prey, Bait), sep = "-_-_")),
    left_join
  ) %>%
    modify_depth(., 1, ~ separate(., combo, c("Prey", "Bait"), sep = "-_-_"))


  walk2(
    saint_tib_f$run,
    saint_tib_f$prohitsviz,
    ~write_tsv(.y %>%
                 select(Bait, Prey, AvgP, BFDR, AvgSpec) %>%
                 mutate(
                   AvgSpec = 50 * 10 * AvgSpec/max(AvgSpec) ## 50 is the standard of prohits viz
                 ),
               paste0("SAINTf/", .x, "/", "for_prohits_viz.txt"),
               col_names = TRUE))


  saint_tib_f <<- saint_tib_f




}
