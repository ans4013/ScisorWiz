#!/usr/bin/env python

###############################################################################
# Follow up script to PlotIsoforms.r only accessed when interactive = "yes":
#   - Input: input file name, output file name
#   - Output: Plot of all isoforms of gene across all cell type choices. File is
#             an interactive html file
#
# Command: ./interactivePlot.R input file name, output file name
#
# by Hagen U. Tilgner (2012); modified and updated by Alexander N. Stein (2021)
################################################################################

if(!require(plotly)) install.packages('plotly')
if(!require(htmlwidgets)) install.packages('htmlwidgets')
if(!require(RCurl)) install.packages('RCurl')

args<-commandArgs(trailingOnly=TRUE);
image_file <- args[1]
txt <- RCurl::base64Encode(readBin(image_file, "raw", file.info(image_file)[1, "size"]), "txt")
p <- plot_ly(x = 0:4, y = 0:4,
             type = 'scatter', mode = 'markers', alpha = 0) %>%
    layout(
        images = list(
            list(
                source =  paste('data:image/jpg;base64', txt, sep=','),
                xref = "x",
                yref = "y",
                x = 0,
                y = 4,
                sizex = 2,
                sizey = 4,
                sizing = "stretch",
                opacity = 1,
                layer = "below"
            )
        )
    )
#m = list(r=0, l=0, b=0, t=0)
fig <- p %>% #layout(margin = m) %>%
    layout(plot_bgcolor='#ffff',
           xaxis = list(showgrid = FALSE,showticklabels = FALSE),
           yaxis = list(showgrid = FALSE,showticklabels = FALSE)
    )
fig
f<-args[2]
htmlwidgets::saveWidget(fig, file.path(normalizePath(dirname(f)),basename(f)))
