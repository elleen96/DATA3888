library(shiny)
library(shinythemes)
library(tidyverse)
library(forecast)
library(shinythemes)
library(GEOquery)  


# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("flatly"),
  # Application title
  titlePanel("Gene Expression Classifiers"),
  sidebarLayout(
    sidebarPanel(
      helpText("This app shows the accuracy of KNN, SVM and Random Forest."),
      
      h4("Model"),
      helpText("You may select which model you want."),
      selectInput(inputId = "model", label = "Please select a model", choices = list("KNN" = 1, "SVM" = 2, "Random Forest" = 3), selected = 1),
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Accuracy Boxplot", plotOutput("c")),
        shiny::plotOutput(outputId = "c")
      )
      
    )
  )
)





fetch_genes = reactive ({
  
  datadir = '/Users/elleentiong/Downloads/GSE120396_RAW/'
  fileNames <- list.files(datadir)
  gse = c()
  for(i in 1:length(fileNames)){
    temptable <- read.delim(file.path(datadir, fileNames[i]), header=TRUE)
    gse <- cbind(gse, temptable[,2])
    colnames(gse)[i] <- colnames(temptable)[2]
  }
  return(gse)
  
})

rejection_status = reactive({
  
  clinical_outcome <-getGEO("GSE120396")
  clinical_outcome<- clinical_outcome$GSE120396_series_matrix.txt.gz
  rejection_status  <- clinical_outcome$characteristics_ch1.1
  rejection_status <- unlist(lapply( strsplit(as.character(rejection_status), ": " ) , `[[` , 2)  )
  return(rejection_status)
  

})

server <- function(input, output) {
   output$c = renderPlot({
   #index <- (input$c)
   #if(index==1) chosen = cv_50acc5_knn
   #if(index==2) chosen = cv_50acc5_svm
  # if(index==3) chosen = cv_50acc5_rf
   #boxplot(chosen)
   

    MATRIX = fetch_genes()
    # 
    largevar = apply(MATRIX, 1, var)
    ind = which(largevar > quantile(largevar, 0.9)) 
    # 
    X = as.matrix(t(MATRIX[ind,]))
    y = rejection_status()
    
    cvK = 5  # number of CV folds
    cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
    cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
    cv_50f15_knn = cv_50f15_svm = cv_50f15_rf = c()
    cv_f1_knn = cv_f1_svm = cv_f1_rf = c()
    
    n_sim = 25
    for (i in 1:n_sim) {
      
      cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
      cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
      cv_f1_knn = cv_f1_svm = cv_f1_rf = c()
      
      for (j in 1:cvK) {
        test_id = cvSets$subsets[cvSets$which == j]
        X_test = X[test_id, ]
        X_train = X[-test_id, ]
        y_test = y[test_id]
        y_train = y[-test_id]
        ## KNN
        fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = 5)
        cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        cv_f1_knn[j] = F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))
        
      
        ## SVM
        svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
        fit <- predict(svm_res, X_test)
        cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        cv_f1_svm[j] = F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))
        
        ## RandomForest
        rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
        fit <- predict(rf_res, X_test)
        cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        cv_f1_rf[j] = F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))
        
      }
      cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
      cv_50acc5_svm <- append(cv_50acc5_svm, mean(cv_acc_svm))
      cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
      cv_50f15_knn <- append(cv_50f15_knn, mean(cv_f1_knn))
      cv_50f15_svm <- append(cv_50f15_svm, mean(cv_f1_svm))
      cv_50f15_rf <- append(cv_50f15_rf, mean(cv_f1_rf))
    } 
 
    
    var(cv_50acc5_knn)
    var(cv_50acc5_svm)
    var(cv_50acc5_rf)
    
    index <- (input$model)
    if(index==1) boxplot(list(Accuracy = cv_50acc5_knn , "F1 score"= cv_50f15_knn)) 
    if(index==2) boxplot(list(Accuracy = cv_50acc5_svm , "F1 score"= cv_50f15_svm)) 
    if(index==3) boxplot(list(Accuracy = cv_50acc5_rf , "F1 score" = cv_50f15_rf)) 
    

  
})  
}



# Run the application 
shinyApp(ui = ui, server = server, options = list(height = 1080))

