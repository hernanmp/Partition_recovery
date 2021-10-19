worstcasematrix = function(n){
  B = matrix(0.25,n,n)
  for (i in 1:n){
    for (j in 1:n){
      if (i > j){B[i,j] = 0.75}
    }
  }
  return(B)
}

dc = function(l,y,lambda){
  n = 2^l
  #opt would be a list of list of matrices. 
  ##opt[[d1]][[d2]] would be a matrix of stored opt values for nodes
  #in depths d1 and d2. 
  opt = list()
  #split = list()
  sum = list()
  sumsq = list()
  size = list()
  var = list()
  partition = list()
  ##partition is a list of lists. Each vertex has a list which denotes the optimal
  ##partition of that vertex. A partition is denoted by a list of two tuples 
  ##representing the end points of the intervals.
  for (i1 in 1:(l + 1)){
    opt[[i1]] = list()
    #split[[i1]] = list()
    sum[[i1]] = list()
    sumsq[[i1]] = list()
    var[[i1]] = list()
    size[[i1]] = list()
    for (i2 in 1:(l + 1)){
      #opt[[i1]] = list()
      opt[[i1]][[i2]] = matrix(rep(0,2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
      #split[[i1]] = list()
      #split[[i1]][[i2]] = matrix(rep(0,2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
      #sum[[i1]] = list()
      sum[[i1]][[i2]] = matrix(rep(0,2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
      #sumsq[[i1]] = list()
      sumsq[[i1]][[i2]] = matrix(rep(0,2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
      #var[[i1]] = list()
      var[[i1]][[i2]] = matrix(rep(0,2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
      #size[[i1]] = list()
      size[[i1]][[i2]] = matrix(rep(2^{2*l-i1-i2+2},2^{i1 + i2 - 2}),2^{i1 - 1},2^{i2 -1})
    }
  }
  #opt[[l + 1]] = list()
  opt[[l + 1]][[l + 1]] = matrix(rep(lambda,2^{2*l}),2^l,2^l)
  #split[[l + 1]][[l + 1]] = matrix(rep(0,2^{2*l}),2^l,2^l)
  sum[[l + 1]][[l + 1]] = y
  sumsq[[l + 1]][[l + 1]] = y^2
  size[[l + 1]][[l + 1]] = matrix(rep(1,n^2),n,n)
  var[[l + 1]][[l + 1]] = matrix(rep(0,n^2),n,n)
  
  ##defining partition. partition[[d1]][[d2]][[i]][[j]] is a list.
  ##This list encodes the optimal partition of the rectangle which
  ##a product of the ith interval at depth d1 in the first coordinate
  ##and a product of the jth interval at depth d2 in the second coordinate. 
  ##the optimal partition is encoded by a list of rectangles. Each
  ##rectangle is encoded by the product of two intervals. Each interval
  ##is encoded by its depth and its number from left to right in the
  ##binary tree.
  partition = list()
  for (d1 in 1:(l + 1)){
    partition[[d1]] = list()
    for (d2 in 1:(l + 1)){
      partition[[d1]][[d2]] = list()
    }
  }
  #partition[[l+1]] = list()
  #partition[[l + 1]][[l + 1]] = list()
  for (i in 1:n){
    partition[[l + 1]][[l + 1]][[i]] = list()
    for (j in 1:n){
      partition[[l + 1]][[l + 1]][[i]][[j]] = list(c(l+1,l+1,i,j))
    }
  }
  ####Figure out how to remove for loop here. Or merge this with the bottom
  ##up pass.
  
  ##bottom up pass
  ##Traverse the grid zigzag 
  A = NULL
  for (add in 2:(2*(l + 1))){
    for (i in 1:min((add - 1),(l + 1))){
      if (add - i < l + 2){
        A = rbind(A,c(i,(add - i)))
      }
    }
  }
  #############################
  for (row in 2:nrow(A)){
    i1 = A[row,1]
    i2 = A[row,2]
    ##d1,d2 are depths
    d1 = l + 2 - i1
    d2 = l + 2 - i2
    for (i in 1:(2^{d1 - 1})){
      partition[[d1]][[d2]][[i]] = list()
      for (j in 1:(2^{d2 - 1})){
        if (d1 == l + 1){
          sum[[d1]][[d2]][i,j] = sum[[d1]][[d2 + 1]][i,2*(j - 1) + 1] + sum[[d1]][[d2 + 1]][i,2*(j - 1) + 2]
        }
        if (d1 < l + 1){
          sum[[d1]][[d2]][i,j] = sum[[d1 + 1]][[d2]][2*(i - 1) + 1,j] + sum[[d1 + 1]][[d2]][2*(i - 1) + 2,j]
        }
        if (d1 == l + 1){
          sumsq[[d1]][[d2]][i,j] = sumsq[[d1]][[d2 + 1]][i,2*(j - 1) + 1] + sumsq[[d1]][[d2 + 1]][i,2*(j - 1) + 2]
        }
        if (d1 < l + 1){
          sumsq[[d1]][[d2]][i,j] = sumsq[[d1 + 1]][[d2]][2*(i - 1) + 1,j] + sumsq[[d1 + 1]][[d2]][2*(i - 1) + 2,j]
        }
        if (d1 < l + 1){
          size[[d1]][[d2]][i,j] = size[[d1 + 1]][[d2]][2*(i - 1) + 1,j] + size[[d1 + 1]][[d2]][2*(i - 1) + 2,j]
        }
        if (d1 == l + 1){
          size[[d1]][[d2]][i,j] = size[[d1]][[d2 + 1]][i,2*(j - 1) + 1] + size[[d1]][[d2 + 1]][i,2*(j - 1) + 2]
        }
        var[[d1]][[d2]][i,j] = sumsq[[d1]][[d2]][i,j] - (sum[[d1]][[d2]][i,j])^2/size[[d1]][[d2]][i,j]
        temp1 = Inf
        if (d1 < l + 1){
          temp1 = opt[[d1 + 1]][[d2]][2*(i - 1) + 1,j] + opt[[d1 + 1]][[d2]][2*(i - 1) + 2,j]
        }
        temp2 = Inf
        if (d2 < l + 1){
          temp2 = opt[[d1]][[d2 + 1]][i,2*(j - 1) + 1] + opt[[d1]][[d2 + 1]][i,2*(j - 1) + 2]
        }
        opt[[d1]][[d2]][i,j] = min(var[[d1]][[d2]][i,j] + lambda,temp1,temp2)
        if (var[[d1]][[d2]][i,j] + lambda <= min(temp2,temp1)){
          #split[[d1]][[d2]][i,j] = 0
          partition[[d1]][[d2]][[i]][[j]] = list(c(d1,d2,i,j))
        }
        else if (min(var[[d1]][[d2]][i,j] + lambda,temp2) >= temp1){
          #split[[d1]][[d2]][i,j] = 1
          partition[[d1]][[d2]][[i]][[j]] = c(partition[[d1+1]][[d2]][[2*(i-1)+1]][[j]],partition[[d1+1]][[d2]][[2*i]][[j]])
        }
        else if (min(var[[d1]][[d2]][i,j] + lambda,temp1) >= temp2){
          #split[[d1]][[d2]][i,j] = 2
          partition[[d1]][[d2]][[i]][[j]] = c(partition[[d1]][[d2+1]][[i]][[2*(j-1)+1]],partition[[d1]][[d2+1]][[i]][[2*j]])
        }
      }
    }
  }
  ####Need to compute the fit#######
  finalpart = partition[[1]][[1]][[1]][[1]]
  output = matrix(rep(0,n,n),n,n)
  for (i in 1:length(finalpart)){
    xstart = (finalpart[[i]][3] - 1)*(n/(2^{finalpart[[i]][1] - 1})) + 1
    xend = (finalpart[[i]][3])*(n/(2^{finalpart[[i]][1] - 1}))
    ystart = (finalpart[[i]][4] - 1)*(n/(2^{finalpart[[i]][2] - 1})) + 1
    yend = (finalpart[[i]][4])*(n/(2^{finalpart[[i]][2] - 1}))
    output[xstart:xend,ystart:yend] = matrix(rep(mean(y[xstart:xend,ystart:yend])),xend - xstart + 1,yend - ystart + 1) 
  }
  return(list(finalpart,output))
}


twopiecemat = function(n){
  twopiece1 = matrix(0,n,floor(n/2))
  twopiece2 = matrix(1,n,floor(n/2))
  twopiece = cbind(twopiece1,twopiece2)
  return(twopiece)
}

worstcasematrix = function(n){
  B = matrix(0.25,n,n)
  for (i in 1:n){
    for (j in 1:n){
      if (i > j){B[i,j] = 0.75}
    }
  }
  return(B)
}


sinematrix = function(n){
  ans = matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      ans[i,j] = sin(pi*i/n)*sin(pi*j/n)
    }}
  return(ans)
}

##given a square matrix, computes its tv
isotvcompute = function(theta){
  n = nrow(theta)
  a = 0
  for (i in 1:(n - 1)){
    for (j in 1:(n - 1)){
      a = a + abs(theta[i + 1,j] - theta[i,j])
      a = a + abs(theta[i,j + 1] - theta[i,j])
    }
  }
  v = theta[n,]
  v1 = v[2:n]
  v2 = v[1:(n - 1)]
  a = a + sum(abs(v2 - v1))
  v = theta[,n]
  v1 = v[2:n]
  v2 = v[1:(n - 1)]
  a = a + sum(abs(v2 - v1))
  return(a)
}



isotvd = function(truth,iter){
  if (!require("Rmosek")) {
    stop ("Rmosek not installed.")
  }
  n = nrow(truth)
  v = isotvcompute(truth)
  mse = rep(0,iter)
  #n = nrow(y)
  #Convert y to a vector 
  #y = c(t(y))
  library(Matrix)
  #n = nrow(y)
  theta.index <- seq(1,n^2)
  alphapos.index <- seq(n^2 + 1, n^2 + 2*n*(n - 1))
  alphaneg.index <- seq(n^2 + 2*n*(n - 1) + 1, n^2 + 4*n*(n - 1))
  alphanum = 4*n*(n - 1)
  qo1 <- list()
  qo1$sense <- "min"
  #qo1$c <- c(-2*y,rep(0,alphanum))
  #vrow = seq(1:(2*m))
  #vrow = vrow/2
  #vrow = ceiling(vrow)
  qo1$qobj <- list(i = c(seq(1:(n^2))),
                   j = c(seq(1:(n^2))),
                   v = c(rep(2,n^2)))
  ###Creating the edge incidence matrix D
  rowind = seq(1:(n*(n - 1)))
  colind = seq(1:n^2)
  delind = n*seq(1:n)
  colind = colind[-delind]
  rowind = rep(rowind,2)
  colind = c(colind,colind + 1)
  Drow = sparseMatrix(i = rowind,j = colind,
                      x = c(rep(-1,n*(n - 1)),rep(1,n*(n - 1))))
  ###
  Dcol = sparseMatrix(i = rep(seq(1,(n - 1)*n),2),j = c(seq(1,(n - 1)*n),seq(from = n + 1, to = n^2)),
                      x = c(rep(-1,(n - 1)*n),rep(1,(n - 1)*n)))
  
  D = rBind(Drow,Dcol)
  id = sparseMatrix(i = seq(1:(nrow(D))),j = seq(1:(nrow(D))),
                    x = c(rep(1,nrow(D))))
  mat = cBind(D,-id,id)
  tvcon = c(rep(0,n^2),rep(1,4*n*(n - 1)))
  qo1$A <- rBind(mat,tvcon)
  qo1$bc <- rBind(blc = c(rep(0,nrow(mat)),0),
                  buc = c(rep(0,nrow(mat)),v))
  qo1$bx <- rBind(blx = c(rep(-Inf,n^2),rep(0,4*n*(n - 1))),
                  bux = c(rep(Inf,n^2),rep(Inf,4*n*(n - 1))))
  #thetastar = twopiecemat(n)
  #se = rep(0,iter)
  for (i in 1:iter){
    z = matrix(rnorm(n^2),n,n)
    y = truth + z
    y = c(t(y))
    qo1$c <- c(-2*y,rep(0,alphanum))
    invisible(capture.output(r <- mosek(qo1)))
    fit <- r$sol$itr$xx[theta.index]
    #alphapos = r$sol$itr$xx[alphapos.index]
    #alphaneg = r$sol$itr$xx[alphaneg.index]
    #alphapos = round(alphapos,3)
    #alphaneg = round(alphaneg,3)
    fit = matrix(fit,n,n,byrow = T)
    #thetahat = fit
    mse[i] = mean((fit - truth)^2)
  }
  #}
  return(mse)
}

smoothing = function(y,lprime,l)
{
  #lprime = 3
  n = 2^l
  eta= 2^lprime
  m =  n/eta
  ytilde =  matrix(0, m,m)
  
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      ytilde[i,j] =  mean(y[(eta*(i-1)+1):(eta*i),  (eta*(j-1)+1):(eta*j) ])
    }
  }
  return(ytilde)
}

recover_sol = function(theta,lprime,l)
{
  n = 2^l
  theta_hat =  matrix(0,n,n)
  eta= 2^lprime
  m = n/eta
  
  
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      theta_hat[(eta*(i-1)+1):(eta*i),  (eta*(j-1)+1):(eta*j) ] =  theta[i,j]
    }
  }
  return( theta_hat)
}

##############################
merge_Dcart =  function(theta_hat,lambda_tilde,gam)
{
  #lambda_tilde = lambda*2
  n = dim(theta_hat)[1]
  Lambda_tilde = c()
  vals = unique(as.vector(theta_hat))
  for(j in 1:length(vals))
  {
    Lambda_tilde[[j]] =  which(theta_hat==vals[j])
  }
  
  E = matrix(0, length(vals),length(vals))
  for(j in 1:(length(vals)-1))
  {
    for(i in (j+1):length(vals))
    {
      #################
      ai = min(Lambda_tilde[[i]])
      bi = max(Lambda_tilde[[i]])
      
      ai1 = ai - floor(ai/n)*n
      if(ai1 ==0) ai1 =n
      ai2 = min(floor(ai/n)+1,n)
      
      bi1 = bi - floor(bi/n)*n
      if(bi1 ==0) bi1 =n
      bi2 = min(floor(bi/n)+1,n)
      
      pi  = list()
      pi[[1]] = c(ai1,ai2)
      pi[[2]] = c(bi1,ai2)
      pi[[3]] = c(ai1,bi2)
      pi[[4]] = c(bi1,bi2)
      
      #######
      aj = min(Lambda_tilde[[j]])
      bj = max(Lambda_tilde[[j]])
      
      aj1 = aj - floor(aj/n)*n
      if(aj1 ==0) aj1 =n
      aj2 = min(floor(aj/n)+1,n)
      
      bj1 = bj - floor(bj/n)*n
      if(bj1 ==0) bj1 =n
      bj2 = min(floor(bj/n)+1,n)
      
      pj  = list()
      pj[[1]] = c(aj1,aj2)
      pj[[2]] = c(bj1,aj2)
      pj[[3]] = c(aj1,bj2)
      pj[[4]] = c(bj1,bj2)
      ###############   
      
      temp = Inf
      for(l in 1:4)
      {
        for(k in 1:4)
        {
          temp2 = sqrt(sum((pi[[l]]-pj[[k]])^2))
          if(temp2 < temp)
            temp = temp2
        }
      }
      if(temp < gam )
      {
        a = length(Lambda_tilde[[j]])*length(Lambda_tilde[[i]])
        a = a/(length(Lambda_tilde[[j]])+length(Lambda_tilde[[i]])) 
        a = a*(  mean(y[Lambda_tilde[[i]]])-mean(y[Lambda_tilde[[j]]])  )^2
        if(a<lambda_tilde )
        {
          E[i,j] = 1
          E[j,i] = 1 
        }
      }
    }
  }
  G = graph_from_adjacency_matrix(E,mode="undirected")
  clu = components(G)
  temp = groups(clu)
  
  Lambda_hat =  list()
  for(j in 1:length(temp))
  {
    for(l in 1:length(temp[[j]]))
    {
      if(l==1)
        Lambda_hat[[j]] =  Lambda_tilde[[temp[[j]][l]]]
      
      if(l>1)
        Lambda_hat[[j]] = c(Lambda_hat[[j]],Lambda_tilde[[temp[[j]][l]]])
    }
  }
  return(Lambda_hat)
}
  
one_sided_haussdorff = function(Lambda0,Lambda_hat)
{
  
  temp = -Inf
  for(i in 1:length(Lambda0))
  {
    temp2 = Inf
    for(j in 1:length(Lambda_hat))
    {
      temp3 =  max(length(setdiff(Lambda0[[i]],Lambda_hat[[j]])),length(setdiff(Lambda_hat[[j]],Lambda0[[i]]))   )
      if(temp3 < temp2)
        temp2 = temp3
      
      # print(temp3)
      # print("")
    }
    if(temp2>temp)
      temp = temp2
  }
  return(temp) 
  
  
}

AIC_path =  function(y,l,lambda_grid)
{
  n = 2^l
  score =  rep(0,length(lambda_grid))
  
 ##aux = rep(0,length(lambda_grid))
  
  for(j in 1:length(lambda_grid))
  {
    ans = dc(l,y,lambda_grid[j])
    theta = ans[[2]]
    #aux[j] = mean((theta-theta0)^2)
    #image(theta)
    sigmahat2 = mean(diff(y)^2)/2 
    
    score[j] = sum((y - theta)^2) + sigmahat2*length(unique(as.vector(theta)))*log(n^2)
      #length(unique(theta))*2
      #sigmahat2*length(unique(theta))*2
    #length(unique(as.vector(theta)))*2
      #sigmahat2*length(unique(as.vector(theta)))*log(n^2)
      #sigmahat2*length(unique(as.vector(theta)))*2
      #length(unique(theta))*2
    
  }
  #best_j = which.max(abs(diff(score)))
  best_j =  which.min(score)
  ans = dc(l,y,lambda_grid[best_j])
  theta = ans[[2]]
  return(list(theta_hat=theta,lambda = lambda_grid[best_j]))
}


generate_data =  function(l,scenario,sigma)
{
  n = 2^l
  
  if(scenario==1)
  {
    theta0 =  matrix(0,n,n)
    
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if(i > n/4 &&  i< 3*n/4 && j<3*n/4  && j>n/4  )
          theta0[i,j] =1
      }
    }
    vals = unique(as.vector(theta0))
    Lambda0 = list()
    for(k in 1:length(vals))
    {
      Lambda0[[k]] = which(theta0==vals[k])
    }
  } 
  if(scenario== 2)
  {
    theta0 =  matrix(0,n,n)
    
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if( (i-n/4)^2+(j-n/4)^2< (n/5)^2   )
          theta0[i,j] =1
        
        if( (i-3*n/4)^2+(j-3*n/4)^2< (n/5)^2   )
          theta0[i,j] =-1
      }
    }
    vals = unique(as.vector(theta0))
    Lambda0 = list()
    for(k in 1:length(vals))
    {
      Lambda0[[k]] = which(theta0==vals[k])
    }
  }
  if(scenario ==3)
  {
    theta0 =  matrix(0,n,n)
    
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if(n/4<i &&  i<3*n/4 &&  j>n/4 &&  j < n/4+n/8)
          theta0[i,j] = 1
        
        if(n/2 + n/8<i && i <3*n/4 && j>= n/4+n/8 && j<3*n/4 )
          theta0[i,j] = 1
        
        #if((i-7*n/8)^2+(j-7*n/8)^2< (n/8)^2 )
        if(i > 6*n/8 && j>6*n/8) 
         theta0[i,j] = -1
      }
    }
    vals = unique(as.vector(theta0))
    Lambda0 = list()
    for(k in 1:length(vals))
    {
      Lambda0[[k]] = which(theta0==vals[k])
    }
  }
  if(scenario ==4)
  {
    theta0 =  matrix(0,n,n)
    
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if(i < n/5  && j <n/5)
          theta0[i,j] = 1
        
        if(i < n/5  && j > 4*n/5)
          theta0[i,j] = 2 
        
        if(i > 4*n/5  && j <n/5)
          theta0[i,j] = 3
        
        if(i > 4*n/5  && j > 4*n/5)
          theta0[i,j] = 4
        
        if( n/2 -n/8 < i &&  i<n/2+n/8 && n/2 -n/8 < j &&  j<n/2+n/8)
          theta0[i,j] = 5
      }
    } 
    vals = unique(as.vector(theta0))
    Lambda0 = list()
    for(k in 1:length(vals))
    {
      Lambda0[[k]] = which(theta0==vals[k])
    }
  }
  if(scenario==5)
  {
    theta0 =  matrix(0,n,n)
    
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        #if(i < n/5  && j <n/5)
        #  theta0[i,j] = 1
        
        if(i < n/5  && j > 2*n/5)
          theta0[i,j] = 2 
        
        if(i > 4*n/5  && j <3*n/5)
          theta0[i,j] = 3
        
        if(abs(i-n/2) <n/4.5 && abs(j-n/2) <n/4.5)
          theta0[i,j] = 4
      }
    }
    vals = unique(as.vector(theta0))
    Lambda0 = list()
    for(k in 1:length(vals))
    {
      Lambda0[[k]] = which(theta0==vals[k])
    }
  }
  #############################
  y = theta0 +   sigma*matrix(rnorm(n*n),n,n)
  
  return(list(theta0 =theta0,y=y,Lambda0=Lambda0))
  
}




ort = function(n,y,lambda){
  #n = 2^m
  #opt would be a list of list of matrices. 
  ##opt[[d1]][[d2]] would be a matrix of stored opt values for nodes
  #in depths d1 and d2. 
  opt = list()
  #split = list()
  sum = list()
  sumsq = list()
  var = list()
  partition = list()
  ##partition is a list of lists. Each vertex has a list which denotes the optimal
  ##partition of that vertex. A partition is denoted by a list of two tuples 
  ##representing the end points of the intervals.
  for (i in 1:n){
    opt[[i]] = list()
    sum[[i]] = list()
    sumsq[[i]] = list()
    var[[i]] = list()
    partition[[i]] = list()
    for (j in 1:n){
      opt[[i]][[j]] = list()
      sum[[i]][[j]] = list()
      sumsq[[i]][[j]] = list()
      var[[i]][[j]] = list()
      partition[[i]][[j]] = list()
      for (k in 1:n){
        opt[[i]][[j]][[k]] = list()
        sum[[i]][[j]][[k]] = list()
        sumsq[[i]][[j]][[k]] = list()
        var[[i]][[j]][[k]] = list()
        partition[[i]][[j]][[k]] = list()
        for (l in 1:n){
          opt[[i]][[j]][[k]][[l]] = list()
          sum[[i]][[j]][[k]][[l]] = list()
          sumsq[[i]][[j]][[k]][[l]] = list()
          var[[i]][[j]][[k]][[l]] = list()
          partition[[i]][[j]][[k]][[l]] = list()
        }
      }
    }
  }
  #This part should be done more efficiently. 
  
  for (i in 1:n){
    for (j in 1:n){
      opt[[i]][[i]][[j]][[j]] = lambda
      sum[[i]][[i]][[j]][[j]] = y[i,j]
      sumsq[[i]][[i]][[j]][[j]] = y[i,j]^2
      var[[i]][[i]][[j]][[j]] = 0
      partition[[i]][[i]][[j]][[j]] = list(c(i,i,j,j))
    }
  }
  
  
  
  #It should be possible to do the above more efficiently. 
  ##Traverse the grid zigzag 
  A = NULL
  for (add in 2:(2*n)){
    for (i in 1:min((add - 1),n)){
      if (add - i < n + 1){
        A = rbind(A,c(i,(add - i)))
      }
    }
  }
  ##################################
  for (row in 2:nrow(A)){
    n1 = A[row,1]
    n2 = A[row,2]
    for (scan1 in 1:(n - n1 + 1)){
      for (scan2 in 1:(n - n2 + 1)){
        ##Rectangle(scan1,scan1 + n1 - 1,scan2,scan2 + n2 -1)
        array1 = rep(Inf,n1)
        array2 = rep(Inf,n2)
        ##the last elements of array woud always be infinity
        if (n1 > 1){
          sum[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sum[[scan1]][[scan1 + n1 - 2]][[scan2]][[scan2 + n2 - 1]] + sum[[scan1 + n1 - 1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]
          sumsq[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sumsq[[scan1]][[scan1 + n1 - 2]][[scan2]][[scan2 + n2 - 1]] + sumsq[[scan1 + n1 - 1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]
          var[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sumsq[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]/(n1*n2) - (sum[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]/(n1*n2))^2
        }
        if (n1 == 1){
          sum[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sum[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 2]] + sum[[scan1]][[scan1 + n1 - 1]][[scan2 + n2 - 1]][[scan2 + n2 - 1]]
          sumsq[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sumsq[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 2]] + sumsq[[scan1]][[scan1 + n1 - 1]][[scan2 + n2 - 1]][[scan2 + n2 - 1]]
          var[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = sumsq[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]/(n1*n2) - (sum[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]/(n1*n2))^2
        }
        ###Remaining code
        if (n1 > 1){
          for (split1 in scan1:(scan1 + n1 - 2)){
            array1[split1 - scan1 + 1] = opt[[scan1]][[split1]][[scan2]][[scan2 + n2 - 1]] + opt[[split1 + 1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]
          }
        }
        if (n2 > 1){
          for (split2 in scan2:(scan2 + n2 - 2)){
            array2[split2 - scan2 + 1] = opt[[scan1]][[scan1 + n1 - 1]][[scan2]][[split2]] + opt[[scan1]][[scan1 + n1 - 1]][[split2 + 1]][[scan2 + n2 -1]]
          }
        }
        opt[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = min(var[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] + lambda,min(array1),min(array2))
        w = opt[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]]
        ####now defining partition
        if (w == var[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] + lambda){
          #split[[d1]][[d2]][i,j] = 0
          partition[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = list(c(scan1,scan1 + n1 - 1,scan2, scan2 + n2 - 1))
        }
        if (w < var[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] + lambda){
          if (w == min(array1)){
            ind = which.min(array1)
            a1 = scan1
            a2 = ind - 1 + scan1
            b1 = ind + scan1
            b2 = scan1 + n1 - 1
            partition[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = c(partition[[a1]][[a2]][[scan2]][[scan2 + n2 - 1]],partition[[b1]][[b2]][[scan2]][[scan2 + n2 - 1]])
          }
          else{
            ind = which.min(array2)
            a1 = scan2
            a2 = ind - 1 + scan2
            b1 = ind + scan2
            b2 = scan2 + n2 - 1
            partition[[scan1]][[scan1 + n1 - 1]][[scan2]][[scan2 + n2 - 1]] = c(partition[[scan1]][[scan1 + n1 - 1]][[a1]][[a2]],partition[[scan1]][[scan1 + n1 - 1]][[b1]][[b2]])
          }
        }
        #####end of the updates for this rectangle.
      }
    }
    #end of updating each rectangle with these sizes
  }
  #end of looping through all possible sizes
  ####Need to compute the fit#######
  finalpart = partition[[1]][[n]][[1]][[n]]
  output = matrix(rep(0,n,n),n,n)
  for (i in 1:length(finalpart)){
    xstart = finalpart[[i]][1]
    xend = finalpart[[i]][2]
    ystart = finalpart[[i]][3]
    yend = finalpart[[i]][4]
    output[xstart:xend,ystart:yend] = matrix(rep(mean(y[xstart:xend,ystart:yend])),xend - xstart + 1,yend - ystart + 1) 
  }
  return(list(finalpart,output))
}
#######The code is working!


AIC_path_ort =  function(y,l,lambda_grid)
{
  score =  rep(0,length(lambda_grid))
  for(j in 1:length(lambda_grid))
  {
    ans = ort(2^l,y,lambda_grid[j])
    theta = ans[[2]]
    #image(theta)
    sigmahat2 = mean(diff(y)^2)/2 
    
    score[j] = sum((y - theta)^2) +  sigmahat2*length(unique(theta))*2
    
  }
  best_j =  which.min(score)
  ans = dc(l,y,lambda_grid[best_j])
  theta = ans[[2]]
  return(list(theta_hat=theta,lambda = lambda_grid[best_j]))
}

AIC_path_tv =  function(y,lambda_grid)
{
  sigmahat2 = mean(diff(y)^2)/2 
  n = dim(y)[1]
  score =  rep(0,length(lambda_grid))
  #aux = rep(0,length(lambda_grid))
  
  for(j in 1:length(lambda_grid))
  {
    out = fusedlasso_weightedl2_admm2d(y,lambda_grid[j])
    #,lambda_grid[j])
    #aux[j] = mean((out$x - theta0)^2)
      #ort(2^l,y,lambda_grid[j])
    theta = out$x
    #theta = round(theta,3)
    #image(theta)
    
    
    score[j] = sum((y - theta)^2) +  sigmahat2*length(unique(as.vector(  round(theta,3)  )))*log(n^2)
    
  }
  best_j =  which.min(score)
  out = fusedlasso_weightedl2_admm2d(y,lambda_grid[best_j])
  theta = out$x
  return(list(theta_hat=theta,lambda = lambda_grid[best_j]))
}


#######################


merge_tv =  function(theta_hat,lambda_tilde,gam)
{
  #lambda_tilde = lambda*2
  n = dim(theta_hat)[1]
  Lambda_tilde = c()
  theta_hat = round(theta_hat,3)
  vals = unique(as.vector(theta_hat))
  c =0
  new_vals = c()
  for(j in 1:length(vals))
  {
    if(length( which(theta_hat==vals[j]))>gam)
    {
      c = c+1
      Lambda_tilde[[c]] =  which(theta_hat==vals[j])    
      #print(length( which(theta_hat==vals[j])))
      #print(vals[j])
      new_vals = c(new_vals,vals[j])
    }
  }
  if(length(Lambda_tilde)==1)
    return(lambda_tilde)
  
  E = matrix(0, c,c)
  for(j in 1:(c-1))
  {
    for(i in (j+1):c)
    {
      #################
      
        #a = length(Lambda_tilde[[j]])*length(Lambda_tilde[[i]])
        #a = a/(length(Lambda_tilde[[j]])+length(Lambda_tilde[[i]])) 
        #a = a*
        a=  (new_vals[i]-new_vals[j])^2
          #a*(  mean(y[Lambda_tilde[[i]]])-mean(y[Lambda_tilde[[j]]])  )^2
        #print(a)
        if(a<lambda_tilde^2)
        {
          E[i,j] = 1
          E[j,i] = 1 
        }

    }
  }
  G = graph_from_adjacency_matrix(E,mode="undirected")
  clu = components(G)
  temp = groups(clu)
  
  Lambda_hat =  list()
  for(j in 1:length(temp))
  {
    for(l in 1:length(temp[[j]]))
    {
      if(l==1)
        Lambda_hat[[j]] =  Lambda_tilde[[temp[[j]][l]]]
      
      if(l>1)
        Lambda_hat[[j]] = c(Lambda_hat[[j]],Lambda_tilde[[temp[[j]][l]]])
    }
  }
  return(Lambda_hat)
}



BIC_path_tv =  function(y,lambda_grid)
{
  sigmahat2 = mean(diff(y)^2)/2 
  n = dim(y)[1]
  score =  rep(0,length(lambda_grid))
  #aux = rep(0,length(lambda_grid))
  
  for(j in 1:length(lambda_grid))
  {
    out = fusedlasso_weightedl2_admm2d(y,lambda_grid[j])
    #,lambda_grid[j])
    #aux[j] = mean((out$x - theta0)^2)
    #ort(2^l,y,lambda_grid[j])
    theta = out$x
    #theta = round(theta,3)
    #image(theta)
    
    
    score[j] = sum((y - theta)^2) +  sigmahat2*length(unique(as.vector(  round(theta,1)  )))*log(n^2)
    
  }
  best_j =  which.min(score)
  out = fusedlasso_weightedl2_admm2d(y,lambda_grid[best_j])
  theta = out$x
  return(list(theta_hat=theta,lambda = lambda_grid[best_j]))
}