library(shiny)
library(shinydashboard)
library(ggplot2)
library(FactoMineR) #Used for standard implementation of PCA
library(dplyr)
library(RSpectra)
require(dplyr)

ui <- dashboardPage(
  dashboardHeader(title = "PCA-Eigen Faces"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("PCA Data", tabName = "data", icon = icon("table")),
      menuItem("PCA Output", tabName = "output", icon = icon("play")),
      menuItem("PCA Variance", tabName = "variance", icon = icon("exchange")),
      menuItem("Eigen Faces", tabName = "eigen", icon = icon("laugh"))
      
    )
  ),
  
  dashboardBody(
    
    tags$head(tags$style(
      HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}')
    )),
    
    tabItems(
      # First tab content
      tabItem(tabName = "data",
              box(width = NULL, status = "primary", title="PCA Input",
                  titlePanel('PCA Dataset'),
                  fileInput('file1','Choose CSV File',
                            accept = c("text/csv", "text/comma-separated-values,text/plain",
                                       ".csv")),
                  textInput(inputId = 'startno', label = 'Start Column Number', placeholder = 'Give the column number of the starting feature'),
                  textInput(inputId = 'endno', label = 'End Column Number', placeholder = 'Give the column number of the ending feature'),
                  textInput(inputId = 'class_label', label = 'Class Feature Column Number', placeholder = 'Give the column number of the feature that defines the classes'),
                  checkboxInput(inputId = 'head', label = 'Head'),
                  checkboxInput(inputId = 'tail', label = 'Tail'),
                  tableOutput('table')
              )
      ),
      
      # Second tab content
      tabItem(tabName = "output",
              titlePanel('Custom PCA Implementation'),
              titlePanel('PCA Data'),
              htmlOutput('textdata'),
              tableOutput('var1'),
              plotOutput('plot1'),
              titlePanel('Standard PCA Implementation'),
              plotOutput('plot2')
      ),
      
      # Third tab content
      tabItem(tabName = "variance",
              titlePanel('Effect of Variance on the selection of Principal Components'),
              sliderInput('var','Variance', min = 0.0, max = 1.0, value = 1.0, step = 0.01, width = '30%'),
              tableOutput('pc')
      ),
      
      # Fourth tab content
      tabItem(tabName = "eigen",
              textInput(inputId = 'imagerow', label = 'Row Number (Value between 1-400)', placeholder = 'Insert the row number of the dataset whose image you want to view'),
              textInput(inputId = 'startimagerow', label = 'Start Row Number for Average Faces', placeholder = 'Insert the starting row number of the dataset whose average image you want to view'),
              textInput(inputId = 'endimagerow', label = 'End Row Number for Average Faces', placeholder = 'Insert the ending row number of the dataset whose average image you want to view'),
              titlePanel('Image 1: Face of the mentioned row  Image 2: Rotated Face of the mentioned row Image 3: Average Faces of the entered Row Numbers'),
              plotOutput('faceplot1'),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              titlePanel('Plot to show the eigen values'),
              plotOutput('faceplot2'),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              titlePanel('First two Eigen Faces'),
              plotOutput('faceplot3'),
      )
    )
  ),
)

server <- function(input, output) { 
  
  dataframe <- reactive({
    data <- read.csv(input$file1$datapath)
    if (input$head == TRUE){
      head(data, n = 10)
    } else if(input$tail == TRUE){
      tail(data, n = 10)
    } else{
      data
    }
  })
  
  pca <- reactive({
    
    data_pca <- dataframe()
    
    pca_calculation <- function(pca_dataset, variance_explained)
    {
      
      #Create an empty list to store the mean, standard deviation, eigen vectors, cumulative variance of each columns
      pca_components = list()
      
      #Initially we calculate the mean of each column and store that to our empty list with the name mean
      pca_components[['mean']]=colMeans(pca_dataset)
      
      #Once we have the mean we compute the standard deviation of the dataset which along with mean is needed for cov matrix
      pca_components[['standard_deviation']]=apply(pca_dataset,2,sd)
      
      #The data is standardized and normalized to get the data values within the range of [0,1]
      stand_pca_data=sweep(pca_dataset,2,pca_components[['mean']])
      stand_pca_data=stand_pca_data%*%diag(1/pca_components[['standard_deviation']])
      
      #The format for calculation covariance matrix here is A*A^T
      #As the dataframe has m = 400 rows n = 40 columns so if we take A^T*A the resultant matrix will be 400*400 so we take A*A^T which gives us 40*40 matrix  
      eigen_cov=eigen(crossprod(stand_pca_data,stand_pca_data))
      
      pca_components[['eigen_values']] <- eigen_cov$values
      
      #Once we have the cov matrix then from the eigen values we can check the contribution of each component towards the variance and the cumulative variance 
      pca_components[['cumvar']] <- cumsum(eigen_cov[['values']])
      
      #The individual contribution of each feature to the variance
      last <- length(pca_components[['cumvar']])
      
      for (i in (1:last)){
        if(i == 1){
          pca_components[['variance_contribution']]<- append(pca_components[['variance_contribution']],(pca_components[['cumvar']][1]*100/pca_components[['cumvar']][last])) 
        }
        if(i > 1){
          pca_components[['variance_contribution']]<- append(pca_components[['variance_contribution']],((pca_components[['cumvar']][i]-pca_components[['cumvar']][i-1])*100/pca_components[['cumvar']][last]))
        }
        i = as.integer(i) + 1
      }
      
      for (variance in pca_components[['cumvar']]){
        pca_components[['cum_variance_contribution']]<- append(pca_components[['cum_variance_contribution']],(variance*100/pca_components[['cumvar']][last])) 
      }
      
      #After the covariance matrix is calculated we can identify the number of required components and the transformation matrix for the same.
      #The column is selected if and only if it satisfies the variance condition 
      pca_components[['reqcomponents']] =sum((pca_components[['cumvar']]/sum(eigen_cov[['values']]))<variance_explained)+1
      
      #print(pca_components[['reqcomponents']])
      
      #We here select the Principal components which are required for the dataset and add them to the final_principal_components group in the pca_components list 
      pca_components[['fin_transform']] =eigen_cov[['vectors']][,1:pca_components[['reqcomponents']]]
      attr(pca_components, "class") <- "final_principal_components"
      
      pca_components
    }
    
    pca1 <- pca_calculation(as.matrix(data_pca[,input$startno:input$endno]),1)
    
    pca1
    
  })
  
  output$textdata <-renderUI({
    HTML('<strong> The rows of the matrix below are: Mean, Standard Deviation, Eigen Values, Individual Variance Contribution, Cumulative Variance Contribution </strong>')
  })
  
  output$table <- renderTable({
    if (!is.null(input$file1)){
      dataframe()
    }
  })
  
  output$var1 <- renderTable({
    data_pca <- dataframe()
    
    pca1 <- pca()
    
    i <- 0
    
    for (var in pca1[['eigen_values']]){
      i <- i+1
    }
    
    colname <- list()
    
    for (itr in 1:i){
      colname <- append(colname,paste0('PC',itr))
    }
    
    rowname <- c('Mean', 'Standard deviation', 'Eigen Values', 'Proportion of Variance', 'Cumulative Proportion')
    
    matx <- matrix(rbind(pca1[['mean']], pca1[['standard_deviation']], pca1[['eigen_values']], pca1[['variance_contribution']], pca1[['cum_variance_contribution']]), nrow=5, dimnames = list(rowname, colname))
    
    matx
  })
  
  output$plot1 <- renderPlot({
    #Mention the path to the dataset
    data_pca <- dataframe()
    
    pca1 <- pca()
    
    predict.pca_calculation <- function(pca_components,pca_data,..)
    {
      #The data is standardized and normalized to get the data values within the range of [0,1]
      stand_pca_data=sweep(pca_data,2,pca_components[['mean']])
      stand_pca_data=stand_pca_data%*%diag(1/pca_components[['standard_deviation']])
      
      #Then we finally return the data multiplied by the transform obtained from pca which gives us the direction in which the data variance increases
      return(stand_pca_data%*%pca_components[['fin_transform']])
    }
    
    pca_projection <- predict.pca_calculation(pca1,as.matrix(data_pca[,input$startno:input$endno]))
    
    xlabel <- paste0('PC1 with variance contribution of : ',pca1[['variance_contribution']][1])
    
    ylabel <- paste0('PC2 with variance contribution of : ',pca1[['variance_contribution']][2])
    
    color_val <- as.integer(input$class_label)
    
    ggplot(data = data_pca)+geom_point(aes(x=-pca_projection[,1],y=pca_projection[,2],color = data_pca[,color_val]))+xlab(xlabel)+ylab(ylabel)+ggtitle('Dataset projected on the main Principal Components')
    
  })
  
  output$plot2 <- renderPlot({
    data_pca <- dataframe()
    
    pca1 <- pca()
    
    pca_data <- PCA(as.matrix(data_pca[,input$startno:input$endno]))
    
    pca_final <- predict(pca_data, as.matrix(data_pca[,input$startno:input$endno]))$coord[,1:2]
    
    xlabel <- paste0('PC1 with variance contribution of : ',pca1[['variance_contribution']][1])
    
    ylabel <- paste0('PC2 with variance contribution of : ',pca1[['variance_contribution']][2])
    
    color_val <- as.integer(input$class_label)
    
    ggplot(data = data_pca)+geom_point(aes(x=pca_final[,1],y=pca_final[,2],color= data_pca[,color_val]))+xlab(xlabel)+ylab(ylabel)+ggtitle('Dataset projected on the main Principal Components')
  })
  
  output$pc <- renderTable({
    data_pca <- dataframe()
    
    pca_calculation <- function(pca_dataset, variance_required)
    {
      
      #Create an empty list to store the mean, standard deviation, eigen vectors, cumulative variance of each columns
      pca_components = list()
      
      #Initially we calculate the mean of each column and store that to our empty list with the name mean
      pca_components[['mean']]=colMeans(pca_dataset)
      
      #Once we have the mean we compute the standard deviation of the dataset which along with mean is needed for cov matrix
      pca_components[['standard_deviation']]=apply(pca_dataset,2,sd)
      
      #The data is standardized and normalized to get the data values within the range of [0,1]
      stand_pca_data=sweep(pca_dataset,2,pca_components[['mean']])
      stand_pca_data=stand_pca_data%*%diag(1/pca_components[['standard_deviation']])
      
      #The format for calculation covariance matrix here is A*A^T
      #As the dataframe has m = 150 rows n = 4 columns so if we take A^T*A the resultant matrix will be 150*150 so we take A*A^T which gives us 4*4 matrix  
      eigen_cov=eigen(crossprod(stand_pca_data,stand_pca_data))
      
      pca_components[['eigen_values']] <- eigen_cov$values
      
      #Once we have the cov matrix then from the eigen values we can check the contribution of each component towards the variance and the cumulative variance 
      pca_components[['cumvar']] <- cumsum(eigen_cov[['values']])
      
      #The individual contribution of each feature to the variance
      last <- length(pca_components[['cumvar']])
      
      for (i in (1:last)){
        if(i == 1){
          pca_components[['variance_contribution']]<- append(pca_components[['variance_contribution']],(pca_components[['cumvar']][1]*100/pca_components[['cumvar']][last])) 
        }
        if(i > 1){
          pca_components[['variance_contribution']]<- append(pca_components[['variance_contribution']],((pca_components[['cumvar']][i]-pca_components[['cumvar']][i-1])*100/pca_components[['cumvar']][last]))
        }
        i = as.integer(i) + 1
      }
      
      for (variance in pca_components[['cumvar']]){
        pca_components[['cum_variance_contribution']]<- append(pca_components[['cum_variance_contribution']],(variance*100/pca_components[['cumvar']][last])) 
      }
      
      #After the covariance matrix is calculated we can identify the number of required components and the transformation matrix for the same.
      #The column is selected if and only if it satisfies the variance condition 
      pca_components[['reqcomponents']] =sum((pca_components[['cumvar']]/sum(eigen_cov[['values']]))< variance_required)+1
      
      #print(pca_components[['reqcomponents']])
      
      #We here select the Principal components which are required for the dataset and add them to the final_principal_components group in the pca_components list 
      pca_components[['fin_transform']] =eigen_cov[['vectors']][,1:pca_components[['reqcomponents']]]
      attr(pca_components, "class") <- "final_principal_components"
      
      pca_components
    }
    
    pca1 <- pca_calculation(as.matrix(data_pca[,input$startno:input$endno]),input$var)
    
    pca1[['fin_transform']]
  })
  
  image.dataframe <- reactive({
    image.data <- read.csv("C://Users//dhruv//Documents//PCA-Eigen//PCA-EigenFaces//face_data.csv",header = F) %>% as.matrix()
    
    image.data
  })
  
  output$faceplot1 <- renderPlot({
    df_new <- data.frame()
    # Read the pixel data for face from the csv file
    image_data <- image.dataframe()
    
    # Function to plot image data
      
    image_plot <- function(x){ image(x, col=grey(seq(0, 1, length=256)))}
    
    # Set the graphical parameters
    par(mfrow=c(2,2))
    par(mar=c(2,2,2,2))
    
    if (!is.null(input$imagerow)){
      # We will display the image(data from the row parsed) from the image dataset which has 400 images
      image_1 <- matrix(as.numeric(image_data[as.integer(input$imagerow), ]), nrow=64, byrow=T)
      image_plot(image_1)
    }
    
    if (!is.null(input$imagerow)){
      straight_image <- t(apply(matrix(as.numeric(image_data[as.integer(input$imagerow), ]), nrow=64, byrow=T), 2, rev))
      image_plot(straight_image)
    }
    
    if (!is.null(input$startimagerow)){
      if (!is.null(input$endimagerow)){
        # All the images in the dataset are then rotated to make them straight and are then saved into a new file
        for(i in 1:nrow(image_data))
        {
          # Each row is looped through to convert the merge the pixels into a matrix and then rotate the resultant image by 90 degree
          straight_image <- as.numeric((apply(matrix(as.numeric(image_data[i, ]), nrow=64, byrow=T), 2, rev)))
          
          # Bind the new data iteratively with the new dataset
          df_new <- rbind(df_new,straight_image)
        }
        
        # Store the newly obtained data in a dataframe
        df_straight_image =as.data.frame(df_new)
        
        # Calculate the average of the data which are to be taken into consideration and parse it to the image_plot function to construct the average image
        average_data <- colMeans(data.matrix(df_straight_image[input$startimagerow:input$endimagerow,]))
        image_plot(matrix(average_data, nrow=64, byrow=T))
      }
    }
    
  }, height = 500, width = 500)
  
  output$faceplot2 <- renderPlot({
    df_new <- data.frame()
    
    # Read the pixel data for face from the csv file
    image_data <- image.dataframe()
    
    # Function to plot image data
    image_plot <- function(x){ image(x, col=grey(seq(0, 1, length=256)))}
    
    # Set the graphical parameters
    par(mfrow=c(2,2))
    par(mar=c(2,2,2,2))
    
    # All the images in the dataset are then rotated to make them straight and are then saved into a new file
    for(i in 1:nrow(image_data))
    {
      # Each row is looped through to convert the merge the pixels into a matrix and then rotate the resultant image by 90 degree
      straight_image <- as.numeric((apply(matrix(as.numeric(image_data[i, ]), nrow=64, byrow=T), 2, rev)))
      
      # Bind the new data iteratively with the new dataset
      df_new <- rbind(df_new,straight_image)
    }
    
    # Store the newly obtained data in a dataframe
    df_straight_image =as.data.frame(df_new)
    
    data_matrix <- data.matrix(df_straight_image)
    
    # Scale the data to standardize and normalize it to values between (0,1)
    scaled_data_new <- scale(data_matrix)
    
    # Use the A^T*A format to calculate the covariance matrix
    cov_mat <- t(scaled_data_new) %*% scaled_data_new / (nrow(scaled_data_new)-1)
    
    #Calculate the eigen values and vectors of the resultant covariance matrix    
    eigs <- Rspectra::eigs(cov_mat, 40, which = "LM")
    
    # Eigenvalues gives us the weight associated with each of the eigen faces
    eigenvalues <- eigs$values
    
    # Eigenvectors gives us the eigen faces
    eigenvectors <- eigs$vectors
    
    par(mfrow=c(1,1))
    par(mar=c(2.5,2.5,2.5,2.5))
    
    #Plotting the first 40 eigen faces values as they are the ones with the maximum weight
    
    y_cord=eigenvalues[1:40]
    
    plot(1:40, y_cord, type="o", log = "y", main="Weight associated with each Eigen Value", xlab="Eigenvalue Number", ylab="Weight")
  }, height = 500, width = 500)
  
  output$faceplot3 <- renderPlot({
    df_new <- data.frame()
    
    # Read the pixel data for face from the csv file
    image_data <- image.dataframe()
    
    # Function to plot image data
    image_plot <- function(x){ image(x, col=grey(seq(0, 1, length=256)))}
    
    # Set the graphical parameters
    par(mfrow=c(2,2))
    par(mar=c(2,2,2,2))
    
    # All the images in the dataset are then rotated to make them straight and are then saved into a new file
    for(i in 1:nrow(image_data))
    {
      # Each row is looped through to convert the merge the pixels into a matrix and then rotate the resultant image by 90 degree
      straight_image <- as.numeric((apply(matrix(as.numeric(image_data[i, ]), nrow=64, byrow=T), 2, rev)))
      
      # Bind the new data iteratively with the new dataset
      df_new <- rbind(df_new,straight_image)
    }
    
    df_straight_image =as.data.frame(df_new)
    
    ## Let's look at the average face, and need to be substracted from all image data
    average_face=colMeans(df_straight_image)
    AVF=matrix(average_face,nrow=1,byrow=T)
    
    data_matrix <- data.matrix(df_straight_image)
    
    # Scale the data to standardize and normalize it to values between (0,1)
    scaled_data_new <- scale(data_matrix)
    
    # Use the A^T*A format to calculate the covariance matrix
    cov_mat <- t(scaled_data_new) %*% scaled_data_new / (nrow(scaled_data_new)-1)
    
    #Calculate the eigen values and vectors of the resultant covariance matrix    
    eigs <- Rspectra::eigs(cov_mat, 40, which = "LM")
    
    # Eigenvalues gives us the weight associated with each of the eigen faces
    eigenvalues <- eigs$values
    
    # Eigenvectors gives us the eigen faces
    eigenvectors <- eigs$vectors
    
    sum(eigenvalues)/sum(eigen(cov_mat)$values)
    
    # Transformation Matrix obtained by multiplying the scaled data with the Eigen Vectors i.e. the Principal Components
    D_new <- scaled_data_new %*% eigenvectors
    
    # Set the graphical parameters the matrix is is 1*2 1 rows and 2 columns and a standard margin of 2 everywhere
    par(mfrow=c(1,2))
    par(mar=c(2, 2, 2, 2))
    
    #Plot the two Eigen faces against each other and see the characteristics they demonstrate 
    for (i in 1:2){
      image_plot(matrix(as.numeric(eigenvectors[, i]),nrow=64,byrow=T))
    }
    
  }, height = 500, width = 500)
  
}

shinyApp(ui, server)