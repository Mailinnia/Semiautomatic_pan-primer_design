options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
writeLines('Loading libraries')
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(purrr)))
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
writeLines('...')
# suppressWarnings(suppressPackageStartupMessages(library(ggplotify)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
# suppressWarnings(suppressPackageStartupMessages(library(plotly)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
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
  df1 <- read.csv(args$config, sep='\t', header = T) 
  named_list <-as.list(df1 %>% pivot_wider(names_from = 'Argument', values_from = 'Value'))
  args <- c(args, named_list)
  return(args)
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

getDTtable <- function(df1, scrollX=F, width='500px', height='75vh', headerAlign = 'dt-head-center', search=F, nowrap=F, filter=F, filter_cols=NULL){
  class_str <- 'cell-border stripe'
  if (nowrap){
    class_str <- paste(class_str, 'nowrap')
  }
  if (filter){
    suppressWarnings(suppressPackageStartupMessages(library(jsonlite)))
    if (is.null(filter_cols)){
      numeric_cols <- which(sapply(df1, is.numeric)) - 1
    } else {
      numeric_cols <- which(sapply(df1, is.numeric)) - 1
      filtered_indices <- which(colnames(df1) %in% filter_cols)-1
      numeric_cols <- intersect(numeric_cols, filtered_indices)
    }

    numeric_cols_js <- toJSON(as.list(numeric_cols), auto_unbox = TRUE)

     
    df2 <- datatable(df1, rownames = F,class = class_str, escape = FALSE, width=width,
                     options=list(searching=search, info=F, lengthChange=F, paginate=F, scrollX=scrollX, headerCallback=JS(headerCallback),
                                  autoWidth=T, scrollY = height, scrollCollapse = TRUE,
                                  search=list(regex=T),
                                  columnDefs = list(list(className = headerAlign, targets = '_all')),
                                  initComplete = JS(
                                    sprintf(
                                      "function(settings, json) {
                                        this.api().columns().every(function(index) {
                                          var column = this;
                                          if (%s.includes(index)) {
                                            var select = $('<select><option value=\"\"></option></select>')
                                              .appendTo($(column.header()))
                                              .on('change', function() {
                                                var val = $.fn.dataTable.util.escapeRegex($(this).val());
                                                column.search('^' + val + '$', true, false).draw();  // Exact match for dropdown
                                              });
                                            column.data().unique().sort().each(function(d, j) {
                                              select.append('<option value=\"' + d + '\">' + d + '</option>');
                                            });
                                          }
                                        });
                                      }",
                                      numeric_cols_js
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

########################################## varVamp figs ########################################## 
getFigfiles <- function(args){
  files_list <- list.files(args$post_dir, recursive = F, full.names =T)
  
  fig_list <- c()
  for (file in files_list){
    if (str_detect(file, '.svg')){
      fig_list <- c(fig_list, file)
      
    } 
  }
  return(rev(fig_list))
}

getPrimers <- function(args){
  files_list <- list.files(args$out_dir, recursive = F, full.names =T)
  
  for (file in files_list){
    if (str_detect(file, '_primers.tsv')){
        primer_file <- file
    } 
  }
  
  df1 <- read.csv(primer_file, sep='\t', header=T)
  df2 <- df1 %>% select(-c(gc_best, temp_best, mode)) %>%  mutate(on_target_fails = case_when(
    str_count(on_target_fails, ";") < 3 ~ on_target_fails,
    TRUE ~ map_chr(
      str_split(on_target_fails, ";"),
      ~ paste(paste(.[1:3], collapse = ";"), "...")
    )
  ), off_targets = case_when(
    str_count(off_targets, ";") < 3 ~ off_targets,
    TRUE ~ map_chr(
      str_split(off_targets, ";"),
      ~ paste(paste(.[1:3], collapse = ";"), "...")
    )
  )
  )
  return(df2)
}

getStatus <- function(args){
  files_list <- list.files(args$post_dir, recursive = F, full.names =T)
  
  for (file in files_list){
    if (str_detect(file, '_status.tsv')){
      return(file)
      
    } 
  }
}



########################################## Report ########################################## 
makeReport <- function(args){
  wd <- getwd()
  writeLines('Creating primer plots')
  varvamp_plots <- getFigfiles(args)
  primer_df <- getDTtable(getPrimers(args), width='95%', height = '65vh', scrollX = T, filter=T, search=T, filter_cols = c('threshold'))
  status_df <- getDTtable(read.csv(getStatus(args), sep='\t', header=T) %>% select(-dir), width='550px')
  settings_df <- getDTtable(getSettings(args$config, args), width='700px', scrollX = T)
  html_out <- str_glue('{args$out_dir}/{args$prefix}_primer-design_{args$date}.html')
  
  report <- rmarkdown::render(args$rmd,
                              output_file = html_out,
                              params = list(
                                date = args$date,
                                title = args$report_title,
                                settings = settings_df,
                                varvamp_plots = varvamp_plots,
                                primer_df = primer_df,
                                status_df = status_df
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
  writeLines('\nReport generated')
}

main()


###########################
# setwd("C:/Users/b279093/OneDrive - Sundhedsdatastyrelsen/Github/Semiautomated_primer_design/entropy/jev_complex_auto/post_processing")
# args <- list(prefix='jev_complex')
# 
# entropy_tsv <- 'jev_complex_entropy.tsv' 
# potential_tsv <- 'jev_complex_potential_regions.tsv'
# amplicon_tsv <- 'jev_complex_amplicon_regions.tsv'
# entropy_df <- read.csv(entropy_tsv, sep='\t', header=T)
# potential_df <- read.csv(potential_tsv, sep='\t', header=T)
# amplicon_df <- read.csv(amplicon_tsv, sep='\t', header=T)
# 
# thresh <- 0.80
# 
# 
# shannon_fig <- shannonEntropy(entropy_df, thresh)
# regions_fig <- regionsOverview(potential_df, amplicon_df, thresh)
# 
# sub_fig <- subplot(shannon_fig, regions_fig, nrows=2, shareX=T, shareY = F)
# sub_fig <- sub_fig %>% 
#             layout(
#               title = list(text=str_glue('<b>Consensus threshold {thresh}</b>'), font=list(size=22)),
#               margin = list(t=50)
#               )
# sub_fig
# 
# df1 <- entropy_df %>% filter(threshold==thresh)
# df2 <- potential_df %>% filter(threshold==thresh)
# df3 <- amplicon_df %>% filter(threshold==thresh)
# shannon_figGG <- shannonEntropyGG(df1) + scale_x_continuous(limits = c(0, max(c(df1$pos,df2$end, df3$end))), expand=c(0,0))
# regions_figGG <- regionsPlotGG(df2, df3) + scale_x_continuous( limits = c(0, max(c(df1$pos,df2$end, df3$end))), expand=c(0,0))
# 
# combine_fig <- wrap_plots(grobs=list(shannon_figGG + theme(legend.position = "none"),
#                                 regions_figGG + theme(legend.position = "none")), nrow=2)+ plot_layout(heights= c(4, 1))
# 
# legend1 <- cowplot::get_legend(shannon_figGG + theme(legend.justification = c(0,1)))
# legend2 <- cowplot::get_legend(regions_figGG + theme(legend.justification = c(0,1)))
# combine_legend <- plot_spacer() + legend1 + legend2 + plot_spacer() + plot_layout(ncol = 1, heights = c(0, 0.5,0.5,5))
# sub_fig <- wrap_elements(combine_fig) + combine_legend + plot_layout(widths= c(4, 1))+
#   plot_annotation(
#     title = str_glue('Consensus threshold {thresh}'),  # Add your title here
#     subtitle = str_to_title(args$prefix),  # Optional subtitle
#     # caption = "Data source or notes about the plot",  # Optional caption
#     theme = theme(
#       plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title styling
#       plot.subtitle = element_text(size = 12, hjust = 0.5) # Subtitle styling
#       # plot.caption = element_text(size = 10, hjust = 1)  # Caption styling
#     )
#   )
# sub_fig
