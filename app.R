library(shiny)
library(httr)
library(ggplot2)
library(d3heatmap)
library(RColorBrewer)
library(shinythemes)
library(Matrix)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)


options(repos = BiocManager::repositories())

#loading dataset
data <- readMM("E-GEOD-75367.aggregated_filtered_counts.mtx")
data_col <- read_delim("E-GEOD-75367.aggregated_filtered_counts.mtx_cols", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)
data_row <- read_delim("E-GEOD-75367.aggregated_filtered_counts.mtx_rows", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)

#reshaping data
data <- as.matrix(data)
data_row <- data_row[, -1]
dim(data)

rownames(data)<-data_row$X2
colnames(data)<-data_col$X1

# #assign actual gene names to IDs 
gene.df <- bitr(data_row$X2, fromType = "ENSEMBL",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)



#duplicates can remove one row
gene_symbol <- gene.df[,2]


#subsetting new data
intersection <- gene.df[,1] %in% rownames(data) 
data[intersection, ]
data_subset <- data[gene.df[,1], ]
dim(data_subset)
data_subset <- as.data.frame(data_subset)


rownames(data_subset) = make.names(gene.df$SYMBOL, unique=TRUE)
data_subset<-head(data_subset, n = 50L)

genes <-rownames(data_subset)
data_matrix <- (data_subset)
data_matrix<- log10(as.matrix(data_matrix) + 1)




ui<-fluidPage(
    titlePanel("Heatmap of gene expression from scRNA-seq data of circulating breast cancer cells from 19 women with ER+/HER2- primary tumors"), 
    theme=shinytheme("journal"),
    selectInput(inputId = "normselect", 
                label = "Gene Expression Normalization by: ",
                choices = c("row", "col", "none")),
    textOutput(outputId = "textout"),
    d3heatmapOutput("heatmap", 
                    height="800px", 
                    width="90%"))





server <- function(input, output, session) 
{ 
  output$heatmap <- renderD3heatmap({d3heatmap(x = data_matrix, 
                                               col=brewer.pal(9,"YlGnBu"), 
                                               scale= input$normselect, 
                                               cellnote=data, 
                                               labRow=genes, 
                                               xaxis_font_size=7, 
                                               yaxis_font_size=12, 
                                               show_grid = TRUE,
                                               brush_color = "#0000FF",
                                               height=900)})
                                    
                                    
                          
                  
  output$textout <- renderText({
    x <- "The top 50 highly expressing genes normalized by "
    y <- "The heatmap is interactive! Use the mouse cursor to display row/column/value, highlight 
    a certain gene, or drag a rectangle over the map to zoom in."
    return(paste(x,input$normselect,".",y))
})

}



shinyApp(ui = ui, server = server)

