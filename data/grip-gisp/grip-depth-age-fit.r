file="~/Downloads/grip-gisp/DATA/GRIP/DEPTHAGE/GRIPAGE.DAT"
data=read.table(file, skip=30)
depth=data$V1
age=data$V2
plot(depth, age)

# ls <- loess(age~depth)
# pr.loess <- predict(ls)
# plot(depth, age, "l", las=1, xlab="Depth", ylab="Age")
# lines(pr.loess~depth, col="red", lwd=2)

exp.model=lm(log(age)~ depth)
depthvalues <- seq(0, 3000, 1)
Counts.exponential2 <- exp(predict(exp.model,list(depth=depthvalues)))
plot(depth, age, pch=21,cex=0.5)
lines(depthvalues, Counts.exponential2,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")

c=exp.model$coefficients[[1]]
b=exp.model$coefficients[[2]]
plot(depth, log(age))
lines(depth, b*depth+c, col="red", lwd=3)


# depth vs temp
file="~/Downloads/grip-gisp/DATA/GRIP/PHYSICAL/GRIPTEMP.DAT"
data=read.table(file, skip=39)
depth.new=data$V1
temp=data$V2

# fit
age.new=exp(b*depth.new+c)

plot(age.new, temp)

