
# confidence interval using simulation
N=100;
conf_level =0.95;
x<-runif(N,0,5)
plot(ecdf(x))  # only x is known 
med.hat=median(x)
median(sample(x,N,replace=TRUE))
T_boot_dist<-replicate(
  1e4,median(sample(x,N,replace=TRUE)))
se1<-sd(T_boot_dist)
se2<-sd(replicate(
  1e4,median(runif(N,0,5))))

normal_interval<-c(med.hat-qnorm((1+conf_level)/2)*se1,
             med.hat+qnorm((1+conf_level) / 2)*se1) #approx normal
pivotal_interval<-2*med.hat-quantile(T_boot_dist,
                              c((1+conf_level)/2,(1-conf_level)/2),
                              names=FALSE)
percentile_interval<-quantile(T_boot_dist,c((1-conf_level)/2,(1+conf_level)/2))















