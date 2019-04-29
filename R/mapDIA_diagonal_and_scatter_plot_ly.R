#' after mapDIA has been ran we can look at the results of the single pairwise
#' comparisons. so we need to do some data wrangling to generate the
#' dataframes that go into the draawings. also we set the FDR cut off.
#' we start form the md tibble
#' @export
#'
#'


ibq_format_mapDIA_for_drawing <-
  function(FDRcutoff = 0.01,  listofproteins = list("ACACA", "ADD3"))  {

    # generating a subset tibble for every comparion within
    # the already nested dataframe

    md_tib$mapDIA_nested_bait <-
      md_tib$mapDIA %>%
      map(.,
          ~ nest(.,-Label2) %>%
            separate(Label2,
                     into = c("compare", "base"),
                     sep = "/")
      )

    ### here we prepare the log2 data, removing fragment and peptide info,
    # merge replicates, calculate the mean value per protein and replicate

    md_tib$mapDIA_log2_aggegated_long <-
      md_tib$mapDIA_log2 %>%
      map(.,
          ~.x %>%
            select(-Peptide, -Fragment) %>%
            gather(key, value, -Protein) %>%
            separate(key, into = c("A","B"), ) %>%
            unite(col = "condition", A, B) %>%
            filter(!is.na(value)) %>%
            group_by(Protein, condition) %>%
            summarise( mean = mean(value)) %>%
            ungroup(.)
      )

    ## here we combine the data table with filterd and aggregated log2 intensities
    # of the respective baits

    md_tib$draw <-
      map2(md_tib$mapDIA_nested_bait,
           md_tib$mapDIA_log2_aggegated_long,
           function(.x, .y) {
             df <- .x
             df$log2 <-  list(.y)

             df$log2 <-
               pmap(
                 list(df$compare, df$base, df$log2),
                 ~ ..3  %>%
                   filter(condition == ..1 | condition == ..2)
               )

             df$log2 <-
               df$log2 %>%
               map(., ~ .x %>%
                     spread(condition, mean))

             df$data <-
               map2(df$data, df$log2,
                    ~ left_join(.x, .y))
             df

           })

    ## as all this nesting gets quite obsessive we can now spin off all the info
    # into a new tibble

    md_draw_tib <-
      md_tib %>% select(run, draw) %>% unnest(.)

    # and we mark our cut-offs


    md_draw_tib$data <-
      md_draw_tib$data %>%
      map(.,
          ~mutate(.,
                  color = if_else(FDR <= FDRcutoff, "SIG_COLOR", "DEFAULT_COLOR"),
                  group = if_else(FDR <= FDRcutoff, "Significant", "Not Significant"),
                  protein_listed = if_else(Protein %in% listofproteins, 6, 1.5)
          ))


    md_draw_tib <<- md_draw_tib
  }

#' mapDIA results that have been formatted with the function
#' ibq_format_mapDIA_for_drawing can be visualized as an interactiv plot_ly
#' graph with this function
#' @export
#' @import plotly

ibq_mapDIA_draw_diagonal <-
  function(bait, base, ddata,
           DEFAULT_COLOR = "#1976d2",
           SIG_COLOR = "#ff6600",
           protein_listed = list("YWHAE", "YWHAZ") ) {

    maxX = max(ddata[[base]])
    minX = min(ddata[[base]])

    interactivePlot = plot_ly(
      color = ~group,
      colors = c(DEFAULT_COLOR, SIG_COLOR),
      data = ddata,
      mode = "markers",
      text = ~Protein,
      type = "scatter",
      size = ~protein_listed,
      x = ~ ddata[[base]],
      y = ~ ddata[[bait]]  ) %>%
      layout(
        xaxis = list(title = paste(base), zeroline = FALSE),
        yaxis = list(title = paste(bait), zeroline = FALSE),
        shapes = list(
          type='line',
          x0 = minX,
          x1 = maxX,
          y0 = minX,
          y1 = maxX,
          line = list(color = '#757575', dash = 'dot', width = 1)
        )
      )
    interactivePlot
  }

#' mapDIA results that have been formatted with the function
#' ibq_format_mapDIA_for_drawing can be visualized as an interactiv plot_ly
#' graph with this function from https://www.biostars.org/p/214100/
#' @export
#' @import plotly
ibq_mapDIA_draw_volcano <-
  function(bait, base, ddata) {
  # make the Plot.ly plot
  p <- plot_ly(data = ddata,
               x = ddata$log2FC,
               y = -log10(ddata$FDR),
               text = ddata$Protein,
               mode = "markers",
               color = ddata$group) %>%
    layout(title ="Volcano Plot") #%>%
  #layout(annotations = a)
  p
}



