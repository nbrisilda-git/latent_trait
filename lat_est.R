library(statnet)
library(ggplot2)
library(igraph)
library(lme4)
library(nlme)
library(lavaan)
library(plm)
library(latentnet)
library("node2vec")
library(dplyr)
library(purrr)
library(amen)

#########################################
#200 simulation for one scenario, network:N=40, 
#time point T=5, 
# dimension d=3
#######################################

#############################
#n1: N=40,T=5
#############################
# Create dynamic network and behavior


data<-list(NULL)#need for SDNE to read data 
for ( y in 1:3)
  
{
  data[[y]]<-list()
  
  for (x in 1:1)# 1 sim
  {
    
    
    ##### Create a dynamic network with 40 nodes and high homophily
    g<-network(40,directed=TRUE,density=0.2)
    U<-list(NULL)
    for (i in 1:20)
    {U[[i]]<-matrix(0,40,40)}
    U[[1]]<-as.sociomatrix(g)
    
    # Create exposure matrix in influence model 
   
    A<-list(NULL)
    for (i in 1:20)
    {A[[i]]<-matrix(0,40,40)}
    
    
    s<-matrix(0,40,40)
    for (t in 1:40)
    {
      if (sum(U[[1]][t,])!=0)
        s[t,t]<-(0.3*y-0.2)/(sum(U[[1]][t,])) #(0.3*y-0.2)=beta_2
      
    }
    
    ww<-s%*%U[[1]]
    ss<-matrix(0,40,40)
    diag(ss)<-(1 - 0.3*y)#=beta_1
    A[[1]]<-ww+ss
    
    
    Y<-matrix(0,40,20)
    z<-rnorm(40,0,1)
    c<-rnorm(40,0,1)
    Y[,1]<-rnorm(40,0,2)
    
    for (k in 2:20)
    {
      Y[,k]<-(A[[k-1]])%*%Y[,k-1]+ 0.1*z + 0.1 * c  + rnorm(40,0,0.2)
      
      for (i in 1:40)
      {
        for (j in 1:40)
        {
          if (i!=j)
            # selection model
            U[[k]][i,j]= -0.25 -0.3*(abs(z[i]-z[j])) - 0.3 *(abs(c[i]-c[j]))+rnorm(1,0,1)>0
        }
      }
      s<-matrix(0,40,40)
      for (t in 1:40)
      {
        if (sum(U[[k]][t,])!=0)
          s[t,t]<-(0.3*y-0.2)/(sum(U[[k]][t,]))
        
      }
      
      ww<-s%*%U[[k]]
      ss<-matrix(0,40,40)
      diag(ss)<-(1 - 0.3*y)
      A[[k]]<-ww+ss
      
    }
    
    #sum(U[[20]])
    
    # Create exposure term 
    E<-matrix(0,40,19)
    for (j in 1:19)
    {
      for (i in 1:40)
      {
        if (sum(U[[j]][i,])!=0)
          E[i,j]<-(U[[j]][i,]%*%Y[,j])/sum(U[[j]][i,])
      }
      
    }
    
    Dep<-as.vector(t(Y[,16:20]))
    Prior<-as.vector(t(Y[,15:19]))
    Expo<-as.vector(t(E[,15:19]))
    
    infl<-data.frame(cbind(Dep,Prior,Expo))
    
    t<-rep(c(1:5),40)
    id<-rep(c(1:40),each=5)
    
    infl$t <- t
    infl$id <- id
    infl$id <- as.factor(infl$id)
    
    unobs<-rep(c,each=5)
    infl$unobs<-unobs
    cov<-rep(z,each=5)
    infl$cov<-cov
    

    ############################## Estimation Methods #######
    
    #1 :#############Latent adjusted  d=3 ###############
    X<-list(NULL)
    zz<-data.frame(z)
    g<-list(NULL)
    for (i in 1:5) #5 time points
    {
      g[[i]]<-network(U[[i+14]],vertex.attr=zz,directed=TRUE)
      X[[i]]<-ergmm(g[[i]] ~ euclidean(d = 3)+absdiff("z"),control=ergmm.control(burnin=20000,sample.size= 500,interval=5))$mkl$Z  # dimension: d=3 Run latent space model, get latent space estimates
      
    }
    
    df <- do.call("cbind",X)
    
    df<-df[rep(seq_len(nrow(df)), each = 5), ]
    infl<-cbind(infl, df)
    

    
    #2: ##############Latent Factor  d=3 ###############
    
    Xd<-matrix(,nrow=length(z),ncol=length(z)) 
    for(j in 1:length(z)){
      for(i in 1:length(z)){
        if(i!=j){
          Xd[i,j] = abs(z[i]-z[j])#generate dyad covariates
        }
      }
    }
    
    ameU<-list(NULL)
    ameV<-list(NULL)
    fitAME<-list(NULL)
    for (i in 1:2){
     
      fitAME[[i]] = ame(Y<-data[[y]][[x]][[i]],
                        Xdyad=Xd, # incorp dyadic covariates
                        symmetric=FALSE, # tell AME trade is directed
                        intercept=TRUE, # add an intercept             
                        family='bin', # model type
                        rvar=FALSE, # sender random effects (a)
                        cvar=FALSE, # receiver random effects (b)
                        dcor=TRUE, # dyadic correlation
                        R=3, # 3 dimensional multiplicative effects
                        nscan=500, burn=2000, odens=5,
                        plot=FALSE, print=FALSE, gof=TRUE)
      
      
      #gets estimatin 
      ameU[[i]]<-fitAME[[i]]$U[,] #latent est
      
      ameV[[i]]<-fitAME[[i]]$V[,] #latent est
    }
    dU <- do.call("cbind",ameU)
    
    dU<-dU[rep(seq_len(nrow(dU)), each = 2), ]
    
    dV<-do.call("cbind",ameV)
    dV<-dV[rep(seq_len(nrow(dV)), each = 2), ]
    
    
    infl<-cbind(infl[,-c(4:6)], dU,dV)
    names(infl)[c(11:16)]<-letters[1:6]#change come col names
    
    
 
    exp<-coef(summary(lm(Dep~.,data=infl)))[3,1]
  
  
    #3: ######################### SDNE d=3 #########################
  
    
    #get embedding from python) 
    data[[y]][[x]]<-list()
    for (l in 1:5)
    {
      data[[y]][[x]][[l]] <- read.table(paste("C:/file_", y,'-',x,'-',l, ".csv", sep=""))
    }
    mDF<-reduce(data[[y]][[x]],full_join, by='V1')
    
    mDF[] <- lapply(mDF, gsub, pattern =",", replacement = "")
    mDF[] <- lapply(mDF, as.numeric)
    mDF<-mDF[order(as.numeric(mDF$V1)),]
    mDF<-mDF[rep(seq_len(nrow(mDF)), each = 5), ]
    
    infl<-cbind(infl[,-c(4:6)], mDF[,-c(1)])#use only column nedded
    
   
    exp<-coef(summary(lm(Dep~.,infl)))[3,1]
  
    
    #4: node2vec d=3###################################
    #read embedding=latent
    data[[y]][[x]]<-list()
    for (l in 1:5)
    {
      data[[y]][[x]][[l]] <- read.table(paste("C:/file_", y,'-',x,'-',l, ".txt", sep=""), skip=1 ) #path
    }
    DF<-reduce(data[[y]][[x]],full_join, by='V1') #merge all 5 itme points
    DF<-DF[order(DF$V1),]
    DF<-DF[rep(seq_len(nrow(DF)), each = 5), ]   
    infl<-cbind(infl[,-c(4:6)], DF[,-c(1)])
    
    
    exp<-coef(summary(lm(Dep~.,infl)))[3,1]
    
  }
}



###########################################################################################

