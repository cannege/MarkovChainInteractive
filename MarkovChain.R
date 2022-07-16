
library(expm)

OperII.1 <- function() {
  
  options(digits = 6)
  
  # Number of states (N) (integer between [1,10]) (number of rows and columns of the matrix)
  N_states <- as.integer(readline(prompt = "Enter the number of states between 1 and 10 (N x N matrix) : "))
  
  while((N_states<1) | (N_states>10)){
    
    N_states <- as.integer(readline(prompt = "Number of states must be an integer between 1 and 10.Please enter a valid value : "))
  }
  
  
  matrix.n <- matrix(NA,N_states,N_states)
  t_prob <- " "
  i.1 = 1
  print("For Quit enter Q.")
  
  # Construction of matrix by user input(Transition Probabilities)
  while(i.1 < N_states+1){
    if (t_prob == "Q"){break}
    i.2 = 1
    
    while(i.2 < N_states+1){ 
      print(paste("Enter for row ",i.1," column ",i.2," :" ))
      t_prob <- readline()
      if (t_prob == "Q"){break}
      
      while ((t_prob < 0) | (t_prob > 1)) {
        print("Transition probablity cannot be greater than 1 or less than 0.")
        t_prob<- readline(prompt = "Please enter a valid value : ")
      }
      t_prob <- as.double(t_prob)
      matrix.n[i.1,i.2] <- t_prob
      i.2 = i.2+1
    }
    
    if (rowSums(matrix.n)[i.1] != 1){
      print(paste("Sum of probabilities in row ",i.1," is not equal to 1 ,Construct row ",i.1," again"))
    }
    else {i.1 = i.1+1}
  } 
  
  print("Input matrix : ")
  print(matrix.n) 
  
  # n - step transition matrix (to calculate any n-step transition probability matrix (TPM) 
  print("Enter the transition step number to calculate n-step transition probability matrix between 1 and 100 : ")
  step.n <- as.integer(readline())
  
  while ((1>step.n) | (step.n >100)) {
    
    print("Please enter an integer between 1 and 100.")
    step.n <- as.integer(readline())
  }
  
  print(paste("Input matrix at step ",step.n))
  
  matrix.n1 <- matrix.n%^%step.n
  print(matrix.n1)
  
  # If is there any Absorbing States or not.
  abs_states <- c()
  matrix_N <- matrix.n
  
  for (i.3 in N_states:1) {
    if ( matrix.n[i.3,i.3] == 1){ abs_states <- append(abs_states,i.3) 
    matrix_N <- matrix_N[-i.3,]
    matrix_N <- matrix_N[,-i.3]

    }   
  }
  # If there is at least one absorbing states.
  if(length(abs_states) > 0) {
    
          matrix_A <- matrix(NA,N_states,length(abs_states))
          col.lenA <- ncol(matrix_A)
          
          for (i.4 in abs_states) {
            matrix_A[,col.lenA] <- matrix.n[,i.4]
            col.lenA = col.lenA-1
          }
          
          for (i.4 in abs_states) {
            matrix_A <- as.matrix(matrix_A[-i.4,])
          }
          
          I1 <- diag(N_states-length(abs_states))
          x_one <- matrix(c(rep.int(1,N_states-length(abs_states))),N_states-length(abs_states),1)
          
          b1 <- solve(I1-matrix_N) # (I-N)^-1
          c <- b1%*%x_one # (I-N)^(-1)*1
          d <- b1%*%matrix_A # (I-N)^(-1)*A
          

          cat(paste(" ","Absorbing state(s)  : ",sep="\n"))
          print(abs_states)
          
          cat(paste(" ","Matrix N : ",sep="\n"))
          print(matrix_N)
          cat(paste(" ","Matrix A : ",sep="\n"))
          print(matrix_A)
          cat(paste(" ","(I-N)^-1 : ",sep="\n"))
          print(b1)
          cat(paste(" ","(I-N)^(-1)*1 : ",sep="\n"))
          print(c)
          cat(paste(" ","(I-N)^(-1)*A : ",sep="\n"))
          print(d)
  }
  
  # If there is no absorbing states.
  
    else{
      cat(paste(" ","There is no absorbing state.",sep="\n"))
      
      # Steady-State Probabilities
      I2 <- diag(N_states)
      A1 <- rbind(t(matrix.n-I2),rep.int(1,N_states))
      b2 <- c(rep.int(0,N_states),1) 
      
      cat(paste(" ","Steady-State Probabilities  : ",sep="\n"))
      Steady.State <- as.vector(qr.solve(A1,b2))
      print(Steady.State)
      
      # Mean Recurrence Time
      cat(paste(" ","Mean Recurrence Time  : ",sep="\n"))
      MRT <- c()
      for (i.5 in Steady.State) { MRT <- append(MRT,(1/i.5)) }
      print(MRT)
      
    }
  }
  





  




























