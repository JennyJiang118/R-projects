#by 张志强&姜春妮
#此模型为预测武汉市疫情数据故取2020年1月23日至4月27日共96天的疫情数据进行研究
#添加以下改进：
#1 无症状感染者
#2 康复人群复发


#参数设置
#模型初值设定
S1 = 14186500
#武汉市人口总数
# E1 = 1766
E1 = 1766
#1月23至1月29日新增确诊病例数
# I1 = 441
I1 = 241
#感染者
Sq1 = 2609
#尚在接受医学观察的人数
Eq1 = 376
#估计值，为正在被隔离的潜伏者
# H1 = 817
H1 = 217
#正在住院的患者，为感染者和被隔离的潜伏者之和
R1 = 23
#官方公布的出院人数


#模型参数设定
c = 2
#接触率
deltaI = 0.13
#感染者的隔离速度
deltaq = 0.13
#隔离潜伏者向隔离感染者的转化速率
gammaI = 0.007
#感染者的恢复率
gammaH = 0.014
#隔离感染者的恢复速率
# 1/23 - 2/11 beta_1 = 7.35 * 10 ^ (-9) 
q = 1 * 10 ^ (-6)
#隔离比例

rho = 1
#有效接触系数，参考取1
theta1 = 1
#潜伏者相对于感染者的传染能力比值
lambda = 1 / 14
#隔离接触速度，为14天的倒数
sigma = 1 / 7
#潜伏者向感染者的转化速度，平均潜伏期为7天，为7天的倒数

#-----------change parameters-------------#
#上一模型中，为消除无症状感染者的误差，减小了beta_1,因此在本模型中增大
beta_1=1.3*10^(-8)
#传染概率

#由于无症状感染者的存在，病死率更小
alpha = 4 * 10 ^ (-3)
#病死率


#------------new parameters---------------#

beta_2=0.01
#无症状感染者传染概率，beta_2<beta1
eita_1=0.00005
#康复人群中复发成感染者的概率（包括误诊出院）
eita_2=0.00001
#康复人群中抗体失灵从而变成易感者的概率
miu=0.008#check
#无症状感染者转变成感染者的比例，数据来源于真实数据
#p=0.04
p=0.02
#无症状感染者被检测排查出的概率，反映的是检测速度
k=0.83
#接触者转变成感染者的概率，数据来源于真实数据

IN=H1/k*(1-k)
#无症状感染者，根据H1,k计算所得



#----------------iterate equations------------------------#
#----------------phase 1----------------------------------#

#差分迭代方程
T = 1:20
for (idx in 1:(length(T) - 1))
{
  S1[idx+1]=S1[idx]-(rho*c*beta_1+rho*c*q*(1-beta_1))*S1[idx]*(I1[idx]+theta1*IN[idx])+lambda * Sq1[idx]+eita_2*R1[idx]+theta1*beta_2*(1-miu)*(1-p)*IN[idx]+(1-deltaI)*beta_1*I1[idx]
  #易感人数迭代
  
  IN[idx+1]=IN[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])*(1-k)-miu*IN[idx]
  #无症状感染者人数迭代
  
  I1[idx+1]=I1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])*k+miu*(1-p)*IN[idx]+eita_1*R1[idx]-(deltaI + alpha + gammaI) * I1[idx]
  #感染者人数迭代
  
  Sq1[idx+1]=Sq1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])-lambda*Sq1[idx]
  #隔离易感染着人数迭代
  
  Eq1[idx+1]=Eq1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])-deltaq * Eq1[idx]
  #隔离潜伏者人数迭代
  
  H1[idx + 1] = H1[idx] + deltaI * I1[idx] + deltaq * Eq1[idx]+p*IN[idx] - (alpha +
                                                                              gammaH) * H1[idx]
  #住院患者人数迭代
  
  R1[idx + 1] = R1[idx] + gammaI * I1[idx] + gammaH * H1[idx]-(eita_1+eita_2)*H1[idx]
  #康复人数迭代
}

#I1
#H1
#plot(T,H1+I1)
#H1+I1


# 第一阶段数据 1/23 - 2/11 20天
list1=H1+I1

#----------------phase 2-----------------------------#


# 2/12之后
E1=12326
I1=2581
Sq1 = 116500
Eq1 = 43760 
H1 = 27257
R1 = 2286
beta_1 = 2.05 * 10 ^ (-9)
gammaH = 0.08
gammaI = 0.07

T = 1:76
for (idx in 1:(length(T) - 1))
{
  S1[idx+1]=S1[idx]-(rho*c*beta_1+rho*c*q*(1-beta_1))*S1[idx]*(I1[idx]+theta1*IN[idx])+lambda * Sq1[idx]+eita_2*R1[idx]+theta1*beta_2*(1-miu)*(1-p)*IN[idx]+(1-deltaI)*beta_1*I1[idx]
  #易感人数迭代
  
  IN[idx+1]=IN[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])*(1-k)-miu*IN[idx]
  #无症状感染者人数迭代
  
  I1[idx+1]=I1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])*k+miu*(1-p)*IN[idx]+eita_1*R1[idx]-(deltaI + alpha + gammaI) * I1[idx]
  #感染者人数迭代
  
  Sq1[idx+1]=Sq1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])-lambda*Sq1[idx]
  #隔离易感染着人数迭代
  
  Eq1[idx+1]=Eq1[idx]+rho*c*beta_1*(1-q)*S1[idx]*(I1[idx]+theta1*IN[idx])-deltaq * Eq1[idx]
  #隔离潜伏者人数迭代
  
  H1[idx + 1] = H1[idx] + deltaI * I1[idx] + deltaq * Eq1[idx]+p*IN[idx] - (alpha +
                                                                              gammaH) * H1[idx]
  #住院患者人数迭代
  
  R1[idx + 1] = R1[idx] + gammaI * I1[idx] + gammaH * H1[idx]-(eita_1+eita_2)*H1[idx]
  #康复人数迭代
}

list2=H1+I1
list=c(list1,list2)

#对武汉实际数据进行处理
wuhandata=data1
wuhandata$existing = wuhandata$cum_confirm - wuhandata$cum_heal - wuhandata$cum_dead
#模型取的时间段 2020年1月23日至4月27日共96天的疫情数据
useData=wuhandata[54:149,]
#数据进一步处理
useData <- read.csv("/Users/osx/Desktop/useData.csv")
useData=useData[,-1]

useData$prediction=list
useData$prediction=round(useData$prediction)

#绘图
library(ggplot2)
p1=ggplot(useData)+
  geom_point(aes(x=time,y=existing,fill ="real"),size=2,shape=21,color="green")+
  geom_point(aes(x=time,y=prediction,fill ="model"),size=2,shape=22)+
  geom_line(aes(x=time,y=prediction),color="red",group=1)+
  theme(axis.text = element_text(size = 14),axis.title=element_text(size = 16),title = element_text(size = 20),legend.title = element_text(size = 16), legend.text  = element_text(size = 16))+
  #设置字体
  labs(x="date",y="number",title="2020 Wuhan COVID-19 model forecast",fill="")+
  scale_x_discrete(breaks=c("2020-01-23","2020-02-12","2020-02-18","2020-03-25", "2020-04-27"), labels=c("01/23","02/12", "02/18","03/25","04/27"))
#设置x轴刻度




#数据保存
library(readr)
write.csv(useData,"useData.csv")


