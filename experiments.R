
rm(list = ls())
source("utils.R")
library(igraph)
library(Rcpp)
source("gfl2d_admm.R")


sigma_grid = c(0.5,1,1.5)
l_grid = c(7)
lprime = 3
NMC = 50
num_scenarios = 4

haussdorff_original = array(0,c(num_scenarios,NMC,length(l_grid),length(sigma_grid))) 
haussdorff_tv = array(0,c(num_scenarios,NMC,length(l_grid),length(sigma_grid))) 

error_original = array(0,c(num_scenarios,NMC,length(l_grid),length(sigma_grid))) 
error_tv = array(0,c(num_scenarios,NMC,length(l_grid),length(sigma_grid))) 


average_haussdorff_original = array(0,c(num_scenarios,length(l_grid),length(sigma_grid))) 
average_haussdorff_tv = array(0,c(num_scenarios,length(l_grid),length(sigma_grid))) 

average_error_original = array(0,c(num_scenarios,length(l_grid),length(sigma_grid))) 
average_error_tv = array(0,c(num_scenarios,length(l_grid),length(sigma_grid))) 









iter = 1
ind_sigma = 1
lambda_grid_tv = 10^seq(-1,3,length= 20)
lambda_grid = seq(5,30,length=15)
scenario = 1
ind_l = 1 
for(scenario in 1:num_scenarios)
{
  print("scenario")
  print(scenario)
  for(ind_sigma in 1:length(sigma_grid) )
  {
    sigma = sigma_grid[ind_sigma]
    print("sigma")
    print(sigma)
    for(ind_l in 1:length(l_grid))
    {
      l = l_grid[ind_l]
      for(iter in 1:NMC)
      {
        print("iter")
        temp  = generate_data(l,scenario,sigma)
        y = temp$y
        theta0 = temp$theta0
        Lambda0 = temp$Lambda0
        
        temp =   AIC_path(y,l,lambda_grid)
        theta_hat1 = temp$theta_hat
        lambda =  temp$lambda
        Lambda_hat = merge_Dcart(theta_hat1,lambda,2^lprime)  
        
        n = 2^l
        theta_hat_after = matrix(0,n,n)
        for(j in 1:length(Lambda_hat))
        {
          theta_hat_after[Lambda_hat[[j]]] = mean(y[Lambda_hat[[j]]])
        }
        haussdorff_original[scenario,iter,ind_l,ind_sigma] = one_sided_haussdorff(Lambda0,Lambda_hat)
        error_original[scenario,iter,ind_l,ind_sigma] = abs(length(Lambda0)-length(Lambda_hat))
        
      
        ###############################
        ###  tv
        temp =   AIC_path_tv(y,lambda_grid_tv)
        theta_hat1 = temp$theta_hat
        lambda =  temp$lambda
        Lambda_hat_tv = merge_tv(theta_hat1,0.15,8)
        
        
        haussdorff_tv[scenario,iter,ind_l,ind_sigma] = one_sided_haussdorff(Lambda0,Lambda_hat_tv)
        error_tv[scenario,iter,ind_l,ind_sigma] = abs(length(Lambda0)-length(Lambda_hat_tv))
        
      }
      ###############3
      average_haussdorff_original[scenario,ind_l,ind_sigma] = mean(haussdorff_original[scenario,,ind_l,ind_sigma])  
      average_error_original[scenario,ind_l,ind_sigma] =  mean(error_original[scenario,,ind_l,ind_sigma] )
      
      average_haussdorff_tv[scenario,ind_l,ind_sigma] = mean(haussdorff_tv[scenario,,ind_l,ind_sigma])  
      average_error_tv[scenario,ind_l,ind_sigma] =  mean(error_tv[scenario,,ind_l,ind_sigma] )
      
     
      print("original")
      print(average_haussdorff_original[scenario,ind_l,ind_sigma] )
      print(average_error_original[scenario,ind_l,ind_sigma])

      print("tv")
      print(average_haussdorff_tv[scenario,ind_l,ind_sigma] )
      print(average_error_tv[scenario,ind_l,ind_sigma])
      
    }
  }
}


# scenario =1
# ind_sigma = 1
# print("original")
# print(average_haussdorff_original[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_original[scenario,,ind_l,ind_sigma] ))
# print(average_error_original[scenario,ind_l,ind_sigma])
# print(sd(error_original[scenario,,ind_l,ind_sigma]))
# 
# 
# print("tv")
# print(average_haussdorff_tv[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_tv[scenario,,ind_l,ind_sigma] ))
# print(average_error_tv[scenario,ind_l,ind_sigma])
# print(sd(error_tv[scenario,,ind_l,ind_sigma]))
# 
# 
# ind_sigma = 2
# print("original")
# print(average_haussdorff_original[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_original[scenario,,ind_l,ind_sigma] ))
# print(average_error_original[scenario,ind_l,ind_sigma])
# print(sd(error_original[scenario,,ind_l,ind_sigma]))
# 
# 
# print("tv")
# print(average_haussdorff_tv[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_tv[scenario,,ind_l,ind_sigma] ))
# print(average_error_tv[scenario,ind_l,ind_sigma])
# print(sd(error_tv[scenario,,ind_l,ind_sigma]))
# 
# 
# ind_sigma = 3
# print("original")
# print(average_haussdorff_original[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_original[scenario,,ind_l,ind_sigma] ))
# print(average_error_original[scenario,ind_l,ind_sigma])
# print(sd(error_original[scenario,,ind_l,ind_sigma]))
# 
# 
# print("tv")
# print(average_haussdorff_tv[scenario,ind_l,ind_sigma] )
# print(sd(haussdorff_tv[scenario,,ind_l,ind_sigma] ))
# print(average_error_tv[scenario,ind_l,ind_sigma])
# print(sd(error_tv[scenario,,ind_l,ind_sigma]))



