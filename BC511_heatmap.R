#library(shiny)
library(httr)
library(ggplot2)
library(readr)
library(shinyHeatmaply)
library(d3heatmap)
library(RColorBrewer)
library(shinythemes)
library(readr)
library(Matrix)


data <- read_csv("/Users/ericabrophy/Documents/BMI511_translational/BreastCancerDataset/marker-genes.csv")


row.names(data)<-data$Category
genes <-data$Category
data_matrix <- data[, -1]
data_matrix<- log10(as.matrix(data_matrix) + 1)



ui<-fluidPage(
  titlePanel("Heatmap of gene expression from scRNA-seq data of circulating breast cancer cells from 19 women with ER+/HER2- primary tumors"), 
              theme=shinytheme("cerulean"),
              selectInput(inputId = "normselect", 
                          label = "Gene Expression Normalization by: ",
                          choices = c("row", "col", "none")),
              d3heatmapOutput("heatmap", 
                              height="800px", 
                              width="80%"))


 
server <- function(input, output, session) 
{ output$heatmap <- renderD3heatmap({d3heatmap(x = data_matrix, 
                                                 col=brewer.pal(9,"Reds"), 
                                               scale= input$normselect, 
                                               cellnote=data, 
                                               labRow=genes, 
                                               xaxis_font_size=10, 
                                               yaxis_font_size=10, 
                                               height=900)})}





shinyApp(ui = ui, server = server)


