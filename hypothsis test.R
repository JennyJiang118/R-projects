# z test with known sigma, normal or large samples
library('BSDA')
x = c(159, 280, 101, 212, 224, 379, 179, 264, 222, 362, 168, 250,
      149, 260, 485, 170);
z.test(x, alternative = "greater", mu=225,sigma.x=98,
       conf.level = 0.95)
#标准正态分布
# one of "greater", "less" or "two.sided",单侧/双侧
# sigma.x must be known


# t test with unknown sigma, normal
# t分布检验
t.test(x, alternative = "greater", mu=225, conf.level = 0.95)

# two samples,z-test and t-test
sdx = sd(x);
sdy = sd(y);
z.test(x,y, mu=0,alternative = "less",sigma.x=sdx,sigma.y=sdy)
t.test(x,y,mu=0,alternative = "less",var.equal = TRUE)



# one and two samples, variance
library('DescTools')
VarTest(x,alternative="two.sided",sigma.squared=3,conf.level=0.95)
VarTest(x,y,alternative = "two.sided",ratio=3);


#non-normal
#prop.test can be used for testing the null that the proportions (probabilities of success) in several groups are the same, or that they equal certain given values
prop.test(445, 500, p = 0.85, alternative = "greater")


#非参数统计
# Nonparametric test, ks-test
x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)
# Does x come from another distribution
ks.test(x+2, "pgamma", 3, 2) # two-sided, exact
ks.test(x, "pnorm")
ks.test(x+2, "pgamma", 3, 2, alternative = "gr")

#Nonparametric,lillie.test
LillieTest(rnorm(100, mean = 5, sd = 3))
LillieTest(runif(100, min = 2, max = 4))

#Nonparametric 


