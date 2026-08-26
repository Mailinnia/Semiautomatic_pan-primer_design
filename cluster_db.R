options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
writeLines('Loading libraries')
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
writeLines('...')
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(plotly)))
suppressWarnings(suppressPackageStartupMessages(library(randomcoloR)))
suppressWarnings(suppressPackageStartupMessages(library(jsonlite)))
writeLines('Libraries loaded')

########################################## Colors ########################################## 
SSI_grey_100 <- '#E1E6EA'
SSI_charcoal_100 <- '#6B7585'
SSI_charcoal_80 <- '#89919D'
SSI_charcoal_60 <- '#A6ACB6'
SSI_charcoal_40 <- '#C4C8CE'
SSI_charcoal_20 <- '#E1E3E7'
SSI_red_100 <- '#B61616'
SSI_red_80 <- '#C54545'
SSI_red_60 <- '#D37373'
SSI_red_40 <- '#E2A2A2'
SSI_red_20 <- '#F0D0D0' 
petrol_100 <- '#006c9a'
petrol_80 <- '#3389AE'
petrol_60 <- '#66A7C2'
petrol_40 <- '#99C4D7'
green_100 <- '#96bf21'
green_80 <- '#ABCC4D'
green_60 <- '#C0D97A'
green_40 <- '#D5E5A6'

########################################## Shared functions ##########################################  
str_split_func <- function(x) {
  split_str <- strsplit(x, "\t", perl=TRUE)
  number <- as.numeric(data.frame(split_str)[3,])
  return(number)
}

modifyModebarSunburst <- function(fig, height=900, width=1682){
  js <- c(
    'function copyPlot(gd) {
       notifier = $(".plotly-notifier");
       if (notifier) {
          $("body").append("<div class=plotly-notifier></div>");
          notifier = $(".plotly-notifier");
        }
       notifier.append("<div style=padding:12px;width:180px;background:rgba(67,63,63,0.12);border-radius:5px;text-align:left;margin:3px;>Copying...</div>");
    
      var opts = {format: "png", width: WIDTHPX, height: HEIGHTPX, scale: 2};
      Plotly.toImage(gd, opts).then(async function(url) {
        try {
          const data = await fetch(url);
          const blob = await data.blob();
          await navigator.clipboard.write([ new ClipboardItem({ [blob.type]: blob }) ]);
          console.log("Image copied.");
        } catch (err) {
          console.error(err.name, err.message);
        }
        notifier.append("<br/>");
        notifier.append("<div style=padding:12px;width:180px;background:rgba(67,63,63,0.12);border-radius:5px;text-align:left;margin:3px;>Plot copied to clipboard</div>");
        notifier.fadeOut(2000, function() { $(this).remove(); });
      });
    }
    '
  )
  js <- gsub('HEIGHTPX', height, js)
  js <- gsub('WIDTHPX', width, js)
  Copy_SVGpath <- "M97.67,20.81L97.67,20.81l0.01,0.02c3.7,0.01,7.04,1.51,9.46,3.93c2.4,2.41,3.9,5.74,3.9,9.42h0.02v0.02v75.28 v0.01h-0.02c-0.01,3.68-1.51,7.03-3.93,9.46c-2.41,2.4-5.74,3.9-9.42,3.9v0.02h-0.02H38.48h-0.01v-0.02 c-3.69-0.01-7.04-1.5-9.46-3.93c-2.4-2.41-3.9-5.74-3.91-9.42H25.1c0-25.96,0-49.34,0-75.3v-0.01h0.02 c0.01-3.69,1.52-7.04,3.94-9.46c2.41-2.4,5.73-3.9,9.42-3.91v-0.02h0.02C58.22,20.81,77.95,20.81,97.67,20.81L97.67,20.81z M0.02,75.38L0,13.39v-0.01h0.02c0.01-3.69,1.52-7.04,3.93-9.46c2.41-2.4,5.74-3.9,9.42-3.91V0h0.02h59.19 c7.69,0,8.9,9.96,0.01,10.16H13.4h-0.02v-0.02c-0.88,0-1.68,0.37-2.27,0.97c-0.59,0.58-0.96,1.4-0.96,2.27h0.02v0.01v3.17 c0,19.61,0,39.21,0,58.81C10.17,83.63,0.02,84.09,0.02,75.38L0.02,75.38z M100.91,109.49V34.2v-0.02h0.02 c0-0.87-0.37-1.68-0.97-2.27c-0.59-0.58-1.4-0.96-2.28-0.96v0.02h-0.01H38.48h-0.02v-0.02c-0.88,0-1.68,0.38-2.27,0.97 c-0.59,0.58-0.96,1.4-0.96,2.27h0.02v0.01v75.28v0.02h-0.02c0,0.88,0.38,1.68,0.97,2.27c0.59,0.59,1.4,0.96,2.27,0.96v-0.02h0.01 h59.19h0.02v0.02c0.87,0,1.68-0.38,2.27-0.97c0.59-0.58,0.96-1.4,0.96-2.27L100.91,109.49L100.91,109.49L100.91,109.49 L100.91,109.49z"
  CopyImage <- list( name = "Copy", icon = list( path = Copy_SVGpath, width = 111, height = 123 ), click = htmlwidgets::JS(js))
  
  fig <- fig %>% config(displaylogo = FALSE) %>% 
    config(scrollZoom = TRUE) %>% 
    config(displayModeBar = TRUE ) %>% 
    config(modeBarButtonsToRemove = list('lasso2d','select2d','zoomIn2d','zoomOut2d', 'hoverClosestCartesian',
                                         'hoverCompareCartesian', 'autoScale2d')) %>% 
    config(toImageButtonOptions = list(format = "svg",width = 1600,height = 900)) %>%
    config(modeBarButtonsToAdd = list(CopyImage))
  
  return(fig)
}

modifyModebar <- function(fig, height=1009){
  js <- c(
    'function copyPlot(gd) {
       notifier = $(".plotly-notifier");
       if (notifier) {
          $("body").append("<div class=plotly-notifier></div>");
          notifier = $(".plotly-notifier");
        }
       notifier.append("<div style=padding:12px;width:180px;background:rgba(67,63,63,0.12);border-radius:5px;text-align:left;margin:3px;>Copying...</div>");
    
      var opts = {format: "png", width: 1972, height: HEIGHTPX, scale: 2};
      Plotly.toImage(gd, opts).then(async function(url) {
        try {
          const data = await fetch(url);
          const blob = await data.blob();
          await navigator.clipboard.write([ new ClipboardItem({ [blob.type]: blob }) ]);
          console.log("Image copied.");
        } catch (err) {
          console.error(err.name, err.message);
        }
        notifier.append("<br/>");
        notifier.append("<div style=padding:12px;width:180px;background:rgba(67,63,63,0.12);border-radius:5px;text-align:left;margin:3px;>Plot copied to clipboard</div>");
        notifier.fadeOut(2000, function() { $(this).remove(); });
      });
    }
    '
  )
  js <- gsub('HEIGHTPX', height, js)
  Copy_SVGpath <- "M97.67,20.81L97.67,20.81l0.01,0.02c3.7,0.01,7.04,1.51,9.46,3.93c2.4,2.41,3.9,5.74,3.9,9.42h0.02v0.02v75.28 v0.01h-0.02c-0.01,3.68-1.51,7.03-3.93,9.46c-2.41,2.4-5.74,3.9-9.42,3.9v0.02h-0.02H38.48h-0.01v-0.02 c-3.69-0.01-7.04-1.5-9.46-3.93c-2.4-2.41-3.9-5.74-3.91-9.42H25.1c0-25.96,0-49.34,0-75.3v-0.01h0.02 c0.01-3.69,1.52-7.04,3.94-9.46c2.41-2.4,5.73-3.9,9.42-3.91v-0.02h0.02C58.22,20.81,77.95,20.81,97.67,20.81L97.67,20.81z M0.02,75.38L0,13.39v-0.01h0.02c0.01-3.69,1.52-7.04,3.93-9.46c2.41-2.4,5.74-3.9,9.42-3.91V0h0.02h59.19 c7.69,0,8.9,9.96,0.01,10.16H13.4h-0.02v-0.02c-0.88,0-1.68,0.37-2.27,0.97c-0.59,0.58-0.96,1.4-0.96,2.27h0.02v0.01v3.17 c0,19.61,0,39.21,0,58.81C10.17,83.63,0.02,84.09,0.02,75.38L0.02,75.38z M100.91,109.49V34.2v-0.02h0.02 c0-0.87-0.37-1.68-0.97-2.27c-0.59-0.58-1.4-0.96-2.28-0.96v0.02h-0.01H38.48h-0.02v-0.02c-0.88,0-1.68,0.38-2.27,0.97 c-0.59,0.58-0.96,1.4-0.96,2.27h0.02v0.01v75.28v0.02h-0.02c0,0.88,0.38,1.68,0.97,2.27c0.59,0.59,1.4,0.96,2.27,0.96v-0.02h0.01 h59.19h0.02v0.02c0.87,0,1.68-0.38,2.27-0.97c0.59-0.58,0.96-1.4,0.96-2.27L100.91,109.49L100.91,109.49L100.91,109.49 L100.91,109.49z"
  CopyImage <- list( name = "Copy", icon = list( path = Copy_SVGpath, width = 111, height = 123 ), click = htmlwidgets::JS(js))
  
  fig <- fig %>% config(displaylogo = FALSE) %>% 
    config(scrollZoom = TRUE) %>% 
    config(displayModeBar = TRUE ) %>% 
    config(modeBarButtonsToRemove = list('lasso2d','select2d','zoomIn2d','zoomOut2d', 'hoverClosestCartesian',
                                         'hoverCompareCartesian', 'autoScale2d')) %>% 
    config(modeBarButtonsToAdd = list(CopyImage))
  
  return(fig)
}



########################################## Parsers ##########################################  
parseArgs <- function(){
  parser <- ArgumentParser()
  parser$add_argument("-c","--config", required=TRUE, 
                      help="Config containing parameters for pipeline")
  
  args <- parser$parse_args()
  return(args)
}

getParameters <- function(args){
  df1 <- read.csv(args$config, header = T, sep='\t') 
  named_list <-as.list(df1 %>% pivot_wider(names_from = 'Argument', values_from = 'Value'))
  args <- c(args, named_list)
  return(args)
}

getClusterDf <- function(args){
  cluster_df <- read.csv(args$major_tsv, sep='\t', header=T)
  return(cluster_df)
}

getReducedDB<- function(file_dir){
  files_list <- list.files(file_dir, recursive = F, full.names =T)
  
  tsv_df <- data.frame()
  for (file in files_list){
    if (str_detect(file, '_reduced.tsv')){
      df1 <- read.csv(file, sep='\t', header=T)
      if (nrow(tsv_df) == 0){
        tsv_df <- df1
      } else {
        tsv_df <- rbind(tsv_df, df1)
      }
      
    } 
  }
  return(tsv_df)
}


########################################## DT table ########################################## 

headerCallback <- c(
  "function(thead, data, start, end, display){",
  "  $('th', thead).css('border-bottom', '2px solid #6B7585');",
  "  $('th', thead).css('border-top', '1px solid #6B7585');",
  "  $('th', thead).css('border-right', '1px solid #6B7585');",
  "  $(thead).closest('thead').find('th').eq(0).css('border-left', '1px solid #6B7585');",
  "$('table.dataTable.no-footer').css('border-bottom', '1px solid #6B7585');",
  
  "}"
)

# getDTtable <- function(df1, scrollX=F, width='500px', height='75vh', headerAlign = 'dt-head-center', search=F, nowrap=F, filter=F){
#   class_str <- 'cell-border stripe'
#   if (nowrap){
#     class_str <- paste(class_str, 'nowrap')
#   }
#   if (filter){
#     library(jsonlite)
#     numeric_cols <- which(sapply(df1, is.numeric)) - 1 
#     df2 <- datatable(df1, rownames = F,class = class_str, escape = FALSE, width=width,
#                      options=list(searching=search, info=F, lengthChange=F, paginate=F, scrollX=scrollX, headerCallback=JS(headerCallback),
#                                   autoWidth=T, scrollY = height, scrollCollapse = TRUE,
#                                   search=list(regex=T),
#                                   columnDefs = list(list(className = headerAlign, targets = '_all')),
#                                   initComplete = JS(
#                                     sprintf(
#                                       "function(settings, json) {
#                                         this.api().columns().every(function(index) {
#                                           var column = this;
#                                           if (%s.includes(index)) {
#                                             var select = $('<select><option value=\"\"></option></select>')
#                                               .appendTo($(column.header()))
#                                               .on('change', function() {
#                                                 var val = $.fn.dataTable.util.escapeRegex($(this).val());
#                                                 column.search('^' + val + '$', true, false).draw();  // Exact match for dropdown
#                                               });
#                                             column.data().unique().sort().each(function(d, j) {
#                                               select.append('<option value=\"' + d + '\">' + d + '</option>');
#                                             });
#                                           }
#                                         });
#                                       }
#                                       ",
#                                       toJSON(numeric_cols, auto_unbox = TRUE) # Inject numeric column indices as JSON
#                                     )
#                                   )
#                      ),
#                      # filter = "none"
#                      filter = 'top'
#     )
#   } else{  
#     df2 <- datatable(df1, rownames = F,class = class_str,escape = FALSE, width=width,
#                      options=list(searching=search, info=F, lengthChange=F, paginate=F, scrollX=scrollX, headerCallback=JS(headerCallback), 
#                                   autoWidth=T, scrollY = height, scrollCollapse = TRUE, 
#                                   columnDefs = list(list(className = headerAlign, targets = '_all'))
#                      )
#     )
#   }
#   return(df2)
# }
# 

getDTtable <- function(df1, scrollX=F, width='500px', height='75vh', headerAlign = 'dt-head-center', search=F, nowrap=F, filter=F){
  class_str <- 'cell-border stripe'
  if (nowrap){
    class_str <- paste(class_str, 'nowrap')
  }
  if (filter){
    numeric_cols <- which(sapply(df1, is.numeric)) - 1 
    # searching is True, because otherwise filtering breaks when columns are all numeric. The dom decides whether the search box actually appears
    df2 <- datatable(df1, rownames = F,class = class_str, escape = FALSE, width=width,
                     options=list(searching=TRUE, dom = if (search) "ftip" else "t", 
                                  info=F, lengthChange=F, paginate=F, scrollX=scrollX, headerCallback=JS(headerCallback),
                                  autoWidth=T, scrollY = height, scrollCollapse = TRUE,
                                  search=list(regex=T),
                                  columnDefs = list(
                                    list(orderSequence = c('desc', 'asc'), targets = '_all'),
                                    list(className = headerAlign, targets = '_all')
                                    ),
                                  initComplete = JS(
                                    sprintf(
                                      "
                                        function(settings, json) {
                                          var api = this.api();
                                          api.columns().every(function(index) {
                                            var column = this;
                                            var select = $('<select><option value=\"\"></option></select>')
                                              .appendTo($(column.header()))
                                              .on('change', function() {
                                                var val = $.fn.dataTable.util.escapeRegex($(this).val());
                                                if (val) {
                                                  column.search('^' + val + '$', true, false).draw();
                                                } else {
                                                  column.search('', true, false).draw();   // clear filter
                                                }
                                                  
                                              });
                                    
                                            // Collect unique data
                                            var data = column.data().unique().toArray();
                                    
                                            // Numeric vs text sorting
                                            if (%s.includes(index)) {
                                              data = data.sort(function(a, b) {
                                                return parseFloat(a) - parseFloat(b);
                                              });
                                            } else {
                                              data = data.sort(); // alphabetical
                                            }
                                    
                                            // Append sorted options
                                            data.forEach(function(d) {
                                              if (d !== null && d !== undefined && d !== '') {
                                                select.append('<option value=\"' + d + '\">' + d + '</option>');
                                              }
                                            });
                                          });
                                        }
                                      ",
                                      toJSON(numeric_cols, auto_unbox = TRUE) # Inject numeric column indices as JSON
                                    )
                                  )
                     ),
                     # filter = "none"
                     filter = 'top'
    )
  } else{  
    df2 <- datatable(df1, rownames = F,class = class_str,escape = FALSE, width=width,
                     options=list(searching=search, info=F, lengthChange=F, paginate=F, scrollX=scrollX, headerCallback=JS(headerCallback), 
                                  autoWidth=T, scrollY = height, scrollCollapse = TRUE, 
                                  columnDefs = list(list(className = headerAlign, targets = '_all'))
                     )
    )
  }
  return(df2)
}


########################################## Sunburst ########################################## 
blendColor <- Vectorize(function(fg_hex, opacity, bg_hex="#FFFFFF") {
  # Convert hex to RGB
  fg_rgb <- col2rgb(fg_hex) / 255
  bg_rgb <- col2rgb(bg_hex) / 255
  
  # Blend formula
  result_rgb <- opacity * fg_rgb + (1 - opacity) * bg_rgb
  
  # Convert back to hex
  rgb(result_rgb[1], result_rgb[2], result_rgb[3], maxColorValue = 1)
})

safeDistinctPalette <- function(n) {
 if (n <= 1) {
    # one cluster — return a single consistent color
    return(petrol_100)  # or any color you like
  } else {
    # multiple clusters
    return(distinctColorPalette(n))
  }
}

sunburstSpecies <- function(fig, df1, name=NA_character_, depth=3,major_only=F, domain_x=c(0,1/3), domain_y=c(0.5,1), title=NA){
  set.seed(31)

  unique_species <- sort(unique(df1$species))
  unique_clusters <- unique(df1$cluster)
  n_species <- length(unique_species)
  n_clusters <- length(unique_clusters)

  cluster_colors <- data.frame(labels = c(unique_species, unique_clusters ), 
                                color = c(safeDistinctPalette(n_species),safeDistinctPalette(n_clusters)))
  
  # Summarize the data for the sunburst chart
  if (major_only){
    df2 <- df1 %>% 
      count(species, cluster) %>% mutate(cluster=as.character(cluster)) 
    df3 <- bind_rows(
      # Total level
      df2 %>% reframe(labels = 'Total', n = sum(n)),
      # Species level
      df2 %>%
        group_by(species) %>%
        reframe(labels = species, n = sum(n), parents = 'Total'),
      # Cluster level
      df2 %>%
        group_by(species, cluster) %>%
        reframe(labels = cluster, n = sum(n), parents = paste('Total', species, sep = ' - ')) %>% arrange(as.numeric(cluster))
    )
  } else{
    df2 <- df1 %>%
      count(species, cluster, sub_cluster) %>% mutate(cluster=as.character(cluster)) 
    df3 <- bind_rows(
      # Total level
      df2 %>% reframe(labels = 'Total', n = sum(n)),
      # Species level
      df2 %>%
        group_by(species) %>%
        reframe(labels = species, n = sum(n), parents = 'Total'),
      # Cluster level
      df2 %>%
        group_by(species, cluster) %>%
        reframe(labels = cluster, n = sum(n), parents = paste('Total', species, sep = ' - ')) %>% arrange(as.numeric(cluster)),
      # Species level
      df2 %>%
        rename(labels = sub_cluster, parents = cluster) %>%
        mutate(parents = paste('Total', species, parents, sep = ' - '))
    )
  }
  
  sunburst_data <- df3 %>% left_join(cluster_colors, by = c("labels" = "labels")) %>%  # Join colors for species
    mutate(
      ids = if_else(is.na(parents), labels, paste(parents, labels, sep = ' - ')),
      color = if_else(is.na(color), NA_character_, color)  # Ensure color is NA where missing
    ) %>%
    distinct() %>% arrange(as.numeric(cluster))
  # 
  # if (major_only){
  #   sunburst_data <- sunburst_data %>% filter(!labels %in% unique(df2$sub_cluster))
  # }
  
  
  # Create the sunburst chart
  fig <- fig %>% add_trace(
    data = sunburst_data,
    ids = sunburst_data$ids,
    labels = sunburst_data$labels,
    parents = sunburst_data$parents,
    values = sunburst_data$n,
    type = 'sunburst',
    branchvalues = 'total',
    insidetextorientation='radial',
    marker = list(colors = sunburst_data$color, line=list(width=0)),
    maxdepth=depth,
    sort=F, 
    domain = list(x=domain_x, y=domain_y),
    name=name
  )
  
  if (!is.na(title)){
    fig <- fig %>%
      add_annotations(text = paste0('<b>',title,'</b>'), textangle = 0, x = 0.50, xshift = 0, y =1.1, yshift=0, yref = "paper", xref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE,
                      font = list(size = 14))
  }
  
  return(fig)
}

xyPlacement <- function(n, i){
  if (n %% 2 ==0){
    max_cols <- 2
  } else {
    max_cols <- 3
  }
  
  if (n <= 3) {
    y_shift <- 0.1
  } else {
    y_shift <- 0.05
  }
  
  n_rows <- ceiling(n / max_cols)
  row <- floor((i - 1) / max_cols)
  col <- (i - 1) %% max_cols
  figs_in_row <- if (row == n_rows - 1 && n %% max_cols != 0) n %% max_cols else max_cols
  x_start <- (1 - figs_in_row / max_cols) / 2 + col / max_cols
  x_end <- x_start + 1 / max_cols
  y_start <- 1 - (row + 1) / n_rows
  y_end <- (y_start + 1 / n_rows)-0.05
  dynamic_x <- c(x_start, x_end)
  dynamic_y <- c(y_start,y_end)
  title_x <- mean(dynamic_x)
  title_y <- y_end+y_shift## might need to adjust this
  
  return(list(dynamic_x=dynamic_x, dynamic_y=dynamic_y, title_x=title_x, title_y=title_y, row=row))
}

figHeight <- function(n, px=400){
  if (n %%2==0){
    max_rows <- ceiling(n/2)
  } else {
    max_rows <- ceiling(n/3)
  }
  height <- max_rows*px
  return(height)
}

majorClusters <- function(cluster_df, add_title=F, px =400){
  major_ids <- unique(cluster_df %>% pull(major_id))
  n <- length(major_ids)
  height <- figHeight(n, px=px)
  
  fig <- plot_ly(height=height)
  
  i <- 0
  annotations <- list()
  for (id in major_ids){
    i <- i+1
    major_df <- cluster_df %>% filter(major_id == id)
    writeLines(str_glue('\t{id}'))
    title <- paste0(str_glue('<b>{id}</b>'))
    placement <- xyPlacement(n, i)
    annotations <- append(annotations, list(list(x=placement$title_x, y=placement$title_y, text=title, showarrow=F, xref="paper", 
                                                 yref="paper", xanchor='center', yanchor='top', font = list(size = 18))))
    fig <- sunburstSpecies(fig, major_df, major_only=T, domain_x=placement$dynamic_x, domain_y=placement$dynamic_y, name=title)
  }
  
  sub_fig <- fig %>%
    layout(
      grid = list(columns =3, rows=placement$row),
      margin = list(l = 0, r = 0, b = 0, t = 50),
      annotations = annotations
    )
  
  if (add_title) {
    sub_fig <- sub_fig %>% layout(title='<b>Major clusters</b>')
  }
  
  sub_fig <- modifyModebarSunburst(sub_fig)
  return(sub_fig)
}

subClusters <- function(df1, majorID, add_title=T, depth=4){
  sub_ids <- df1 %>% pull(sub_id) %>% unique() 
  n <- length(sub_ids)
  fig_height <- figHeight(n)
  fig <- plot_ly(height=fig_height)
  
  i <- 0
  annotations <- list()
  for (id in sub_ids){
    i <- i+1
    sub_df <- df1 %>% filter(sub_id == id)
    title <- str_glue('<b>Sub-cluster {id}</b>')
    placement <- xyPlacement(n, i)
    annotations <- append(annotations, list(list(x=placement$title_x, y=placement$title_y, text=title, showarrow=F, xref="paper", 
                                                 yref="paper", xanchor='center', yanchor='top', font = list(size = 18))))
    fig <- sunburstSpecies(fig, sub_df, depth=depth, domain_x=placement$dynamic_x, domain_y=placement$dynamic_y, name=title)
    
  }
  
  if (add_title){
    title <- str_glue('<b>Major cluster {majorID}</b>')
  } else {
    title <- NULL
  }
  sub_fig <- fig %>%
    layout(
      grid = list(columns =3, rows=placement$row),
      margin = list(l = 0, r = 0, b = 0, t = 50),
      annotations = annotations,
      title=list(text=title, font=list(size=22))
    )
  
  sub_fig <- modifyModebarSunburst(sub_fig)
  return(sub_fig)
}


subclusterOverview <- function(args){
  sub_dfs <- read.csv(args$sub_tsv, sep='\t', header=T)
  major_ids <- sub_dfs %>% pull(major_id) %>% unique()
  reduced_dfs <- getReducedDB(args$reduce_dir)
  
  sunburst_list <- list()
  for (id in major_ids){
    name <- paste0(id)
    writeLines(str_glue('\t{id}'))
    writeLines(str_glue('\t\tsubclusters'))
    df1 <- sub_dfs %>% filter(major_id == id )
    df2 <- reduced_dfs %>% filter(major_id == id)
    fig_cluster <- subClusters(df1, id, add_title=F)
    
    writeLines(str_glue('\t\tdataframe'))
    subcluster_df <-getDTtable(subclusterDF(df1), filter=T, search=T, width='90%', height = '45vh')
    
    writeLines(str_glue('\t\treduced db'))
    fig_extraction <- subClusters(df2, id, add_title=F, depth=4)
    
    writeLines(str_glue('\t\treduced db dataframe'))
    reduced_df <- getDTtable(reducedDF(df2), filter=T, search=T, width='90%', height = '45vh')
    
    writeLines(str_glue('\t\t...'))
    fig_list <- list(name =name, fig_cluster = fig_cluster,  fig_extraction=fig_extraction, 
                     subcluster_df=subcluster_df, reduced_df = reduced_df)
    sunburst_list <- append(sunburst_list, list(fig_list))
  }
  return(sunburst_list)
}

########################################## Dataframes ########################################## 
majorDF <- function(cluster_df){
  major_ids <- unique(cluster_df$major_id)
  cluster_counts <- data.frame()
  for (id in major_ids){
    df1 <- cluster_df %>% filter(major_id == id)
    df2 <- df1 %>% select(species, cluster) %>%
      count(species, cluster) %>% mutate(cluster=as.numeric(cluster), major_id = as.numeric(id)) %>% rename(count=n)
    if (nrow(cluster_counts)==0){
      cluster_counts <- df2
    } else {
      cluster_counts <- rbind(cluster_counts, df2)
    }
  }
  cluster_counts <- cluster_counts %>% arrange(major_id, cluster)
  return(cluster_counts)
}

subclusterDF <- function(df1){
  sub_ids <- df1 %>% pull(sub_id) %>% unique()
  cluster_counts <- data.frame()
  for (id in sub_ids){
    sub_df <- df1 %>% filter(sub_id == id)
    df2 <- sub_df %>%
      count(species, cluster, sub_cluster) %>% mutate(cluster=as.numeric(cluster), sub_id = as.numeric(id)) %>% rename(count=n)
    if (nrow(cluster_counts)==0){
      cluster_counts <- df2
    } else {
      cluster_counts <- rbind(cluster_counts, df2)
    }
  }
  # print(head(cluster_counts, n=5))
  cluster_counts <- cluster_counts %>% arrange(sub_id, cluster, sub_cluster)
  return(cluster_counts)
  
}

reducedDF <- function(df1){
  sub_ids <- df1 %>% pull(sub_id) %>% unique()
  cluster_counts <- data.frame()
for (id in sub_ids){
      df2 <- df1 %>% filter(sub_id == id) %>% select(species, cluster) %>%
      count(species, cluster) %>% mutate(cluster=as.numeric(cluster), major_id = as.numeric(id)) %>% rename(count=n)
    if (nrow(cluster_counts)==0){
      cluster_counts <- df2
    } else {
      cluster_counts <- rbind(cluster_counts, df2)
    }
  }
  cluster_counts <- cluster_counts %>% arrange(major_id, cluster)
  return(cluster_counts)
}


sharedClusters <- function(df1) {
  shared_clusters <- df1 %>%
    group_by(major_id, cluster, species) %>%
    summarise(seq_count = n(), .groups = "drop_last") %>%
    summarise(
      shared_species = paste(species, collapse = ", "),
      species_count  = n(),
      sequence_count = paste(seq_count, collapse = "; "),
      .groups = "drop"
    ) %>%
    filter(species_count > 1) %>%
    arrange(desc(species_count), major_id, cluster)
  return(shared_clusters)
}

splitSpecies <- function(df1){
  split_species <- df1 %>% group_by(species, major_id, cluster) %>%
    summarise(seq_count = n(), .groups = "drop_last") %>%
    summarize(
      clusters = paste(sort(unique(cluster)), collapse=', '), 
      cluster_count = n_distinct(cluster), 
      sequence_count = paste(seq_count, collapse = "; "),
      .groups = "drop") %>% 
    filter(cluster_count > 1) %>%
    arrange(desc(cluster_count), major_id, species)
  
  return(split_species)
}

########################################## Scores ########################################## 

scoresDF <- function(args){
  major_scores <- read.csv(args$major_scores, sep='\t', header=T) 
  sub_scores <- read.csv(args$sub_scores, sep='\t', header=T)
  
  return(list(major=major_scores, sub=sub_scores))
}

plotScores <- function(scores){
  if (min(scores$major$score, na.rm = T) < 0){
    major_range <- c(-1,1)
  } else {
    major_range <- c(0,1)
  }
  
  if (max(scores$sub$w_dev, na.rm = T) == 0){
    sub_range <- c(0,1)
  } else {
    sub_range <- c(0,max(scores$sub$w_dev) * 1.2)
  }
  
  fig_major <- plot_ly(
    data = scores$major %>% mutate(major_id = as.factor(major_id)),
    x = ~major_id,
    y = ~score,
    type = "bar",
    text = ~round(score, 3),        # text labels (rounded)
    textposition = "outside"        # place outside bar
  ) %>%
    layout(
      yaxis = list(title = list(text="Silhouette score", font=list(size=22)), range = major_range, tickfont = list(size = 16)),
      xaxis = list(title = list(text="Major cluster", font=list(size=22)), tickfont = list(size = 16))
    )
  
  fig_sub <- plot_ly(
    data = scores$sub %>% mutate(major_id = as.factor(major_id), sub_id = as.factor(sub_id)),
    x = ~major_id,
    y = ~w_dev,
    color = ~sub_id,             # color by sub_id
    type = "bar",
    text = ~round(w_dev, 3),
    textposition = "outside"
  ) %>%
    layout(
      barmode = "group",          # group bars side by side
      yaxis = list(title = list(text="Weighted penalty score", font=list(size=22)), range = sub_range, tickfont = list(size = 16)),
      xaxis = list(title = list(text="Major cluster", font=list(size=22)), type = "category", tickfont = list(size = 16))
    )
  
  
  fig_major <- modifyModebar(fig_major)
  fig_sub <- modifyModebar(fig_sub)
  
  return(list(major=fig_major, sub=fig_sub))
  
}




########################################## Settings ########################################## 

getSettings <- function(config_file, args){
  df1 <- read.csv(config_file, header = T, sep='\t')
  df2 <- df1 %>% filter(Value != 'None', Value != 'False')
  
  arg_df <- stack(args) 
  colnames(arg_df) <- c('Value', 'Argument')
  arg_df <- arg_df %>% filter(Argument != 'rscript', Argument != 'rmd', Argument != 'cluster_dir') %>% select(Argument, Value)
  df3 <- rbind(df2,arg_df) %>% distinct()
  
  return(df3)
}

########################################## Report ########################################## 
makeReport <- function(args){
  wd <- getwd()
  writeLines('Creating figures for scores')
  scores <- scoresDF(args)
  scores_fig <- plotScores(scores)
  major_scores <- getDTtable(scores$major %>% rename('Major id' = major_id, 'Score' = score), width='500px', height='65vh')
  sub_scores <- getDTtable(scores$sub %>% rename('weighted_score' = w_dev, 'mean_score' = mean_dev), 
                           filter=T, search=F, width='700px', height='65vh')
  cluster_df <- getClusterDf(args)
  writeLines('Creating figures for major clusters')
  major_clusters <- majorClusters(cluster_df, px=450)
  writeLines('Wrangling major cluster dataframes')
  major_counts <- getDTtable(majorDF(cluster_df), filter=T, search=F, width='90%', height = '65vh')
  shared_clusters <- getDTtable(sharedClusters(cluster_df),filter=T, search=F, width='90%', height = '65vh')
  split_species <- getDTtable(splitSpecies(cluster_df), filter=T, search=F, width='90%', height = '65vh') 

  writeLines('Creating data for subclusters')
  subcluster_figs<- subclusterOverview(args)
  settings_df <- getDTtable(getSettings(args$config, args), width='700px')
  html_out <- str_glue('{args$out_dir}/{args$prefix}_cluster-db_{args$date}.html')
  
  report <- rmarkdown::render(args$rmd,
                              output_file = html_out,
                              params = list(
                                date = args$date,
                                title = args$report_title,
                                settings = settings_df,
                                major_clusters = major_clusters,
                                major_counts = major_counts,
                                shared_clusters = shared_clusters,
                                split_species = split_species,
                                subcluster_figs = subcluster_figs,
                                scores_fig = scores_fig,
                                major_scores = major_scores,
                                sub_scores = sub_scores,
                                args=args
                              ),
                              envir = new.env(),
                              #intermediates_dir = tempdir()
  )
}

########################################## Main ########################################## 

main <- function(){
  args <- parseArgs()
  args <- getParameters(args)
  makeReport(args)
  writeLines('\nHtml report generated')
}

main()






