# calucluate the maximum likelihood estiomation
# on mu for cauchy distribution


n<-100;
set.seed(40);
x<-rnorm(n,1,2);
#dnorm(x,mean=0,sd=1,log=F) 正态分布概率密度函数 dnorm(z)表示标准正态分布密度函数f(x)在x=z处的函数值
#pnorm 正态分布的分布函数值 pnorm(z)==P(X<=z)
#qnorm 给定概率p后的下分位点
#rnorm n个正态分布随机数构成的向量




#-------------------------------------------------------
#最大似然估计：

# method 1 use optimization solver
mlogl<-function(para,x){
  mu = para[1];     #期望
  sigma = para[2];  #方差
  sum(-dnorm(x,mu,sigma,log=TRUE))#概率密度函数求和
  # 因为nlm是non-linear-minimization
  # 然而MLE求最大
  # 故sum这里加负号
}

# or equivalent
mlogl2<-function(para,x){
    N= length(x);#  n个元素
    mu <- param[1];
    sigma <- param[2];
    -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)#带入正态分布密度函数
}

para.start= c(0,1)
max_likelihood = nlm(mlogl,para.start,x)
print(max_likelihood$estimate);


#不用
# method 2 use package

loglik<-function(para){
  sum(dnorm(x,para[1],para[2],log=TRUE));
}
para.start= c(0,1)
res=maxLik(loglik,para.start);
summary( res )



#-------------------------------------------------------
#矩估计
# calculate the estimation using the method of moments
#library: moments
moment_first = moment(x,1,FALSE); #一阶矩
#x:vector of data
#1: 阶矩
#moment(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE)
#central: 中心阶矩
moment_second = moment(x,2,FALSE);#二阶矩

model<-function(x){
  c(F1=x[1]-moment_first,
    F2=x[1]+x[2]^2-moment_second)
}
#library: rootSolve
res=multiroot(model,start=c(0,1));
print(res$root)







