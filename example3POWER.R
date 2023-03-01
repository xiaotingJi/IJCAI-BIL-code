####################################################################
#没有bootstrap的随机梯度下降法    ##################### 
#################################################################
result<-function(x,y,m){
  #W<-rnorm(300,1,1)
  # Wb<-matrix(data = 0,300,30)
  # for (b in 1:30) {
  #   Wb[,b]<-sample(W,300,replace = TRUE)
  # }
  #生成bootstrap数据
  
  
  #随机梯度下降得到估计值
  btsgd<-function(x,y,error,maxiter,step=0.0001)
  { 
    m<-nrow(x)
    x<-cbind(matrix(1,m,1),x)
    n<-ncol(x) 
    theta<-matrix(rep(0,n),n,1)  #ktheta初始值都设置为0
    theta1<-matrix(rep(0,n),n,1)
    iter<-0   #迭代次数
    k<-0  #第k个样本
    newerror<-1 
    
    while(iter<maxiter|newerror>error){
      
      iter<-iter+1
      k<-sample(1:m, 1, replace = TRUE)
      xk<-x[k,,drop=FALSE]
      yk<-y[k,,drop=FALSE]
      des<-t((xk%*%theta-yk)%*%xk)
      new_theta<-theta-step*des
      #newerror<-t(new_theta-theta)%*%(new_theta-theta)
      theta<-new_theta #/(iter)+(theta*(iter-1))/(iter)   
      if(iter==(maxiter-2)){theta1<-theta}
    }
    costfunction<-(t(x%*%theta-y)%*%(x%*%theta-y))/m
    result1<-list(theta,theta1,iter,costfunction)
    #names(result)<-c('系数','系数','系数','系数','系数','迭代次数','误差')
    result1
  }
  #########################
  #生成数据  100个样本#####
  #########################
  
  
  
  ############################
  ##   500次bootsrap      ###
  ############################
  
  
  result2<-btsgd(x,y,error=1,maxiter=m,step=0.001)
  r1=result2[[1]]
  r2=result2[[2]]
  r=list(r1,r2)
  return(r)
}


x1<-array(runif(400,0,5),c(200,2))#协变量
thet<-c(1.5,1,8) #要估计的参数
y1<-cbind(1,x1)%*%thet+rnorm(200)#模型
#################################################################################
nn<-seq(1000,5000,by=1000)
Power<-array(0,dim=c(length(nn),1))
PVALUE<-array(0,dim=c(length(nn),1))
for(l in 1:length(nn)){
  
  
  library(nleqslv)  
  FindRoot = function(x,gamma1,gamma2,alpha) {
    return(pnorm(x-gamma2)-exp(2*gamma2*x)*pnorm(-x-gamma2)-1+alpha)
  }
  MAX<-300
  max<-100
  Qn1<-rep(0,max)
  Times1<-rep(0,max)
  den<-rep(0,max)
  #n1<-50
  #n2<-50
  #sigma1<-1
  #sigma2<-1
  mu<-1
  d0<-0.007
  rm<-matrix(data=0,max,MAX)
  Z_alpha<-nleqslv(1, FindRoot, control=list(btol=.000001),jacobian=TRUE,method="Newton",gamma1=gamma1,gamma2=-d0,alpha= 0.05)$x
  #hatsigma1<-sigma1
  #hatsigma2<-sigma2
  for(k in 1:max){
    Tn1<-rep(0,MAX)
    Tn2<-rep(0,MAX)
    Tn3<-rep(0,MAX)
    Tn4<-rep(0,MAX)
    print(k)
    hatsigma1<-0.03
    #r1=rnorm(30,10,0.04)
    for(j in 1:MAX){
      #r1=matrix(data=rnorm(30,10,0.04),5,60)
      #sd(r1)
      #x<-sample(x1, replace = TRUE)
      #y<-sample(y1, replace = TRUE)
      
      r12=result(x1,y1,m=nn[l])
      r2=r12[[1]]-r12[[2]]
      r1=abs(r2[2])
      rm[k,j]=r1
      #r1=as.numeric(r2)
      #sd(r1)
      if(j<5){sdk=0.07}else
      {sdk=sd(rm[k,1:j])}
      
      xx=0.01-r1
      den[j]=xx
      hatsigma1<-sdk
      if(((sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1))<=0){
        Tn1[j]<-xx-d0
        Tn2[j]<-xx
        
        #hatsigma1<-(sd(xx))
      }else{
        
        Tn3[j]<-d0-xx
        Tn4[j]<--xx
        #hatsigma1<-(sd(xxx)/sqrt(30))
        
      }}
    Qn1[k]<-(sum(Tn2[1:j]+Tn4[1:j]))/((MAX))+sum(Tn1[1:j]+Tn3[1:j])/(sqrt(MAX)*hatsigma1)
    Times1[k]<-(abs(Qn1[k])>Z_alpha) 
  }
  densityf<-function(y,alpha,c){
    exp(-(y^2-2*alpha*(abs(y-c)-abs(c))+alpha^2)/2)/sqrt(2*pi)-
      alpha*exp(2*alpha*abs(y-c))*(1-pnorm(abs(c)+abs(y-c)+alpha,0,1))
  }   
  YY<-seq(-5,5,length=100)
  Alpha<--mu
  res<-array(0,dim=c(length(YY),length(Alpha)))
  c<-0
  for(j in 1:length(YY)){
    for(k in 1:length(Alpha))
      res[j,k]<-densityf(YY[j],Alpha[k],0)
  }
  #plot(YY,res,type="l",ylab="density",xlab="",col="blue",panel.first=grid(10,10,col="gray70"))
  #lines(density(Qn1),ylab="density",col="red",xlab="")
  Power[l,1]<-mean(Times1)
  #Power[l,2]<-mean(Times2)
  print(l)
  
  ############################################################################################3
  #plot(density(den))
  
  Pvalue<-rep(0,max)
  for (i in 1:max) {
    Pvalue[i]=1-(pnorm(abs(Qn1[i])+d0)-exp(2*(-d0)*abs(Qn1[i]))*pnorm(-abs(Qn1[i])+d0))
  }
  PVALUE[l,1]<-mean(Pvalue)
  
}

plot(nn,Power[,1],type="l",ylab="power",xlab="",col="blue")
plot(nn,PVALUE[,1],type="l",ylab="Pvalue",xlab="",col="red")
