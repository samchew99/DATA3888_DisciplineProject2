library(shiny)
library(tidyverse)
library(plotly)
library(shinythemes)
library(DT)
load("app.RData")


# Define UI for application 
ui <- navbarPage("Kidney KNN Explorer",
  
  # Change app theme
  theme = shinytheme("flatly"),
  
  tabPanel("Classifier Performance",
           sidebarPanel(
             sliderInput(inputId = "k_no",
                         label = "Value of K in KNN classifier:",
                         min = 1,
                         max = 25,
                         value = 10),
             sliderInput(inputId = "feature_no",
                         label = "Number of features:",
                         min = 1,
                         max = 100,
                         value = 10)
           ),
           mainPanel(
             fluidRow(column(12,
                             h1("Classifier Performance"),
                             p("In a K nearest neighbours model, changing the value of K can affect the model performance. 
                               This is because a small value for K will result in noise having a higher influence on the result. 
                               On the other hand, a larger K value makes model construction more computationally expensive. 
                               Furthermore, the number of predictors included in the model will also affect the model performance. 
                               The boxplot below demonstrates how the performance (F1 score) of a KNN classifier will change as the value for K and the number of predictors change. 
                               The KNN model is constructed using highly differential genes as explanatory variables for graft rejection outcome and evaluated using its F1 score obtained through repeated 5-fold cross validation. 
                               Detailed information of these genes can be found in the gene information tab. "),
                             br(),
                             h4("Instructions"),
                             p("Use the sliders on the left to select the number of features and the value of K used to construct the KNN classifier."))),
             fluidRow(align = "center",plotOutput(outputId = "performPlot"))
           )),
  tabPanel("Gene Information",
           mainPanel(
             fluidRow(column(12,
                             h1("Gene Information"),
                             p("The data set used in this app was obtained from the public database GEO, which is accessible at GSE131179 [1]. 
                               The data set was from the 'Landscape of innate immune system transcriptome and acute T cell-mediated rejection of human kidney allografts.' study done by Mueller et al. 
                               The data set consisted of 60466 genes and 34 samples, which were collected via RNA sequencing of kidney allograft biopsy specimens from 34 adult kidney transplant recipients. 
                               Among the 34 specimens, 16 were categorized as Banff acute T cell-mediated rejection (TCMR) and 18 were normal."),
                             br(),
                             h4("Instructions"),
                             p("Hover over the points in the volcano plot or search the table below for details of each differentially expressed gene."))),
             fluidRow(align = "center",plotlyOutput(outputId = "genePlot")),
             fluidRow(dataTableOutput("data")))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  output$genePlot = renderPlotly({
    p <- ggplot(gene_table, aes(logFC,-log10(P.Value), key =hgnc_symbol))+
      geom_point(aes(colour=-log10(P.Value)), alpha=1/3, size=1) +
      scale_colour_gradient(low="blue",high="red")+
      xlab("log2 fold change") + ylab("-log10 p-value")
    ggplotly(p, source = "select", tooltip = c("key","x","y"))
  })
  
  output$performPlot <- renderPlot({
    #Most differentially expressed genes
    highExpGenes <- rownames(gene_table[1:input$feature_no,])
    i <- which(rownames(gse) %in% highExpGenes)
    sig_gse <- gse[i,]
    
    #Cross validation 
    set.seed(2020)
    
    X = as.matrix(t(sig_gse))
    y = rejection_status
    
    cvK = 5  # number of CV folds
    cv_50acc5_knn =c()
    cv_acc_knn =c()
    n_sim = 25 ## number of repeats
    for (i in 1:n_sim) {
      cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
      cv_acc_knn = c()
      
      for (j in 1:cvK) {
        test_id = cvSets$subsets[cvSets$which == j]
        X_test = X[test_id, ]
        X_train = X[-test_id, ]
        y_test = y[test_id]
        y_train = y[-test_id]
        
        ## KNN 
        fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = input$k_no) %>% as.factor() %>% caret::confusionMatrix(., factor(y_test, levels = c("Acute Cellular Rejection (ACR)","Normal/Non Specific")))
        cv_acc_knn[j] = fit5$table %>% F1score()
        
      }
      cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
      
      
    }  
    boxplot(list(KNN = cv_50acc5_knn), main="F1 score of KNN classifier",ylim = c(0,1))
    
  })

  
  output$data = renderDataTable({gene_table[,-c(1,2)] <- round(gene_table[,-c(1,2)],2)
                                gene_table})
}

# Run the application 
shinyApp(ui = ui, server = server)