### Data input
nl <- 5                #number of leaves
data <- matrix(c(44,61,2,16,41,1,100,41,0,36,38,1,13,61,0),ncol=nl) #number of YFP-only, CFP-only, and Co-infected sites in nl leaves

### Main body of bottleneck size estimation
##preparation for a calculation table
K <- 30 #maximum number of founders included in the calculation 
kv <- NULL
lv <- NULL
for (i in 0:K){
kv <- c(kv,rep(i,i+1)) #vector of k
lv <- c(lv,0:i)        #vector of l
}
klv <- kv-lv           #vector of k-l
lklv <- lv*klv         #if single-infection (0) or not (>0)
ln <- (K+2)*(K+1)/2    #line number

## Function for log likelihood for r(n) and lambda (-> rl)
rlmLL <- function(rl){
lambda <- rl[nl+1]
LL <- 0
for (l in 1:nl){
table <- matrix(rep(0,ln*5),ncol=5)
table[,1] <- kv
table[,2] <- lv
table[,3] <- klv
table[,4] <- lklv
table[,5] <- dpois(table[,1],lambda)*dbinom(table[,2],table[,1],rl[l])
pni <- table[1,5] #no infection (k=0) 
py <- sum(table[which(table[,2]>0&table[,3]==0),5])/(1-pni) #YFP-only
pc <- sum(table[which(table[,2]==0&table[,3]>0),5])/(1-pni) #CFP-only
pm <- 1-py-pc #Mixed (co-infected)
logL <- dmultinom(c(data[1,l],data[2,l],data[3,l]),prob=c(py,pc,pm),log=TRUE)
LL <- LL+logL
}
-LL
}

## Maximization of log likelihood to find most likely r(n) and lambda
init <- c(rep(0.5,nl),1)
rlmLL.opt <- optim(init,rlmLL, NULL, method="L-BFGS-B",hessian = TRUE,lower=c(rep(0.01,nl),0.002),upper=c(rep(0.99,nl),10))
rlmLL.opt$par        #r(n) and lambda estimates
v <- solve(rlmLL.opt$hessian)
se <- sqrt(diag(v))
se                   #standard errors

