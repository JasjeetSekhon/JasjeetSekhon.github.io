#Generate some data
a = rnorm(100,160,20)
b = rnorm(100,69,3) + (a - 160)/20
c = a+6*b+rnorm(100,0,22)

data = cbind(a,b,c)
head(data)
covMat = cov(data)

#What is the mahalanobis distance between unit 3 and unit 5?
#Two ways:
sqrt((data[3,] - data[5,])%*%solve(covMat)%*%(data[3,] - data[5,]))

sqrt(mahalanobis(data[3,],data[5,], covMat))

#Let's say that the first 50 units are treated, 
#and next 50 units are control.

trt = c(rep(1,50),rep(0,50))
data = cbind(data,trt)

#Look at covariate balance by treatment
aggregate(data[,1], by = list(trt = data[,4]), summary)
aggregate(data[,2], by = list(trt = data[,4]), summary)
aggregate(data[,3], by = list(trt = data[,4]), summary)

#Find minimum matchings for each treated unit
matches = rep(0,50)
count = 1
for(i in 1:50){
		dists = apply(data[51:100,],1, function(x) 
	   {
	   	mahalanobis(data[i,1:3],x[1:3], covMat)^(1/2)  
       })
    matches[count] = which.min(dists) +50
    count = count + 1  
}

mean(data[1:50,1])
mean(data[matches,1])

    