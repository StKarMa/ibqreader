#' this little charm of a function will help us to group stuff bein in the same folder aka in the same analysis
#'
#'@export


ibq_group_by_folder <- function(){

  grouping_helper <- ## here we use the folders for grouping and combining single runs
    results %>% select(run, short_name) %>% as_tibble %>%
    separate(short_name, c("folder", "injection", "samplename"), sep = "/") %>%
    select(-samplename)

  ## to regain the the folde3r into the name
  data_bare_temp <-
    left_join(grouping_helper, data_bare)

  ## group files that are in the same folder
  data_bare_grouped <-
    grouping_helper %>% group_by(folder) %>% tally()

  ##

  ### this might me a little scetchy and error proen when i am doing it with the peptide data as well
  combine_single_run <-
    function(auswahl, datatype){

      tib <-  filter(data_bare_temp, folder == paste0(auswahl))
      tib <- tib %>% select(datatype)

      modify_depth(
        tib, 2, ~gather(.x, key, value,
                        .x %>% select_if(is.numeric) %>% colnames %>% as.vector )
      ) %>% unnest %>%  group_by( gene_name, key) %>% summarise(value = median(value)) %>%
        spread(key, value, fill = NA)

    }


  if (
    data_bare_grouped %>% filter(folder != "Controls") %>% summarise(mean(n)) > 1
  ) {

    data_bare_grouped$peptide <-
      map(data_bare_grouped$folder, ~combine_single_run(., "peptide"))



    data_bare_grouped$protein <-
      map(data_bare_grouped$folder, ~combine_single_run(., "protein" ) )





          ## by using the grouped data all other scrips schould work
          data_bare <<- data_bare_grouped %>% rename(run = folder)
  }



}



