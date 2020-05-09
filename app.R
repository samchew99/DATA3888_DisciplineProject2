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
  
  tabPanel("Motivation",
           sidebarLayout(position = "right",
                         sidebarPanel(position = "right",
                                      h3("Kidney Transplant Facts"),
                                      uiOutput("myList")),
                         mainPanel(h1("Motivation"),
                                  p("Kidney transplantation is one of the most common transplantation procedures performed in Australia but graft resources for these procedures remain scarce. 
                                    For such reasons, it is important to improve graft allocation techniques by optimising models that are used to differentiate between patients who will potentially experience rejection after transplant and patients who will not."), 
                                  p("K nearest neighbours is a commonly used supervised classification algorithm that assigns new data points to the class most common among their K nearest data points. 
                                    The hyper-parameter K is a positive integer that represents the number of neighbours considered when classifying the data points. 
                                    For example, if K = 1, the new data point will be assigned into the class of its closest neighbour. 
                                    Changing the value of K can affect the model performance as a small value for K will result in noise having a higher influence on the result, while a larger K value makes model construction more computationally expensive."), 
                                  p("Besides that, the number of predictors included in the KNN model will also affect the model performance. 
                                    Having less predictors may worsen model performance as less information is used to classify the data, but it is computationally less rigorous. 
                                    Thus, this app aims to allow users to investigate the impact of the number of features and the difference in hyper-parameters on classification performance, particularly in a K-nearest neighbours model.")
                                  )
           )),
  
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
                             br(),
                             h4("Instructions"),
                             p("Use the sliders on the left to select the number of features and the value of K used to construct the KNN classifier."))),
             fluidRow(align = "center",plotOutput(outputId = "performPlot")),
             br(),
             fluidRow("Figure 1: The boxplot above illustrates the performance (F1 score and accuracy) of the KNN classifier obtained through 5-fold cross validation depending on the number of genes and value of K selected")
           )),
  tabPanel("Gene Information",
           mainPanel(
             fluidRow(column(12,
                             h1("Gene Information"),
                             p("The data set used in this app was obtained from the public database GEO, which is accessible at GSE131179. 
                               The data set was from the 'Landscape of innate immune system transcriptome and acute T cell-mediated rejection of human kidney allografts.' study done by Mueller et al. 
                               The data set consisted of 60466 genes and 34 samples, which were collected via RNA sequencing of kidney allograft biopsy specimens from 34 adult kidney transplant recipients. 
                               Among the 34 specimens, 16 were categorized as Banff acute T cell-mediated rejection (TCMR) and 18 were normal."),
                             br(),
                             h4("Instructions"),
                             p("Hover over the points in the volcano plot or search the table below for details of each differentially expressed gene."))),
             fluidRow(align = "center",plotlyOutput(outputId = "genePlot")),
             br(),
             fluidRow("Figure 2: The volcano plot above shows the differentially expressed genes found in the GSE131179 dataset"),
             br(),
             fluidRow(dataTableOutput("data")),
             br(),
             fluidRow("Figure 3: The table above shows the differentially expressed genes found in the GSE131179 dataset sorted by p-value"))
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
    cv_50acc5_knn = cv50_f1 = c()
    cv_acc_knn = cv_f1= c()
    n_sim = 25 ## number of repeats
    for (i in 1:n_sim) {
      cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
      cv_acc_knn = cv_f1 = c()
      
      for (j in 1:cvK) {
        test_id = cvSets$subsets[cvSets$which == j]
        X_test = X[test_id, ]
        X_train = X[-test_id, ]
        y_test = y[test_id]
        y_train = y[-test_id]
        
        ## KNN 
        fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = input$k_no) %>% as.factor() %>% caret::confusionMatrix(., factor(y_test, levels = c("Acute Cellular Rejection (ACR)","Normal/Non Specific")))
        cv_acc_knn[j] = fit5$table %>% Accuracy()
        cv_f1 [j] = fit5$table %>% F1score()
        
      }
      cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
      cv50_f1 <-  append(cv50_f1, mean(cv_f1))
      
    }  
    boxplot(list("F1 score" = cv50_f1, Accuracy = cv_50acc5_knn), main="Performance of KNN classifier",ylim = c(0,1), ylab = "Score")
    
    
  })

  
  output$data = renderDataTable({gene_table[,-c(1,2)] <- round(gene_table[,-c(1,2)],2)
                                gene_table})
  
  output$myList <- renderUI(HTML("<ul><li>1 in 1400 Australians require kidney transplantation in their lifetime.</li>
                                 <li>Average waiting time for transplant is 3 years and can take up to 7 years.</li>
                                 <li>Almost 1% of patients died while waiting for transplant in 2016.</li>
                                 <li>15% of patients experience acute rejection within 6 months post-transplant.</li></ul>"))
}



# Run the application 
shinyApp(ui = ui, server = server)