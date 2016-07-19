testx <- read.csv("C:\\Users\\Keshev\\Documents\\School\\STA414\\A3\\testx",
			 header=FALSE, sep = " ")

testy <- read.csv("C:\\Users\\Keshev\\Documents\\School\\STA414\\A3\\testy", 
			header=FALSE, sep = " ")

train1x <- read.csv("C:\\Users\\Keshev\\Documents\\School\\STA414\\A3\\train1x",
			 header=FALSE, sep = " ")

train1y <- read.csv("C:\\Users\\Keshev\\Documents\\School\\STA414\\A3\\train1y",
			 header=FALSE, sep = " ")

library(cvTools) 
#Linear Model
ptm <- proc.time()
trX = train1x
train1y = as.matrix(train1y) 
train1x = as.matrix(train1x) 

trX$newcol <- rep(1, nrow(trX)) 
trX = trX[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
trX = as.matrix(trX) 
sqX = t(trX)%*%trX
invsqX = solve(sqX)
invsqX
weights = invsqX%*%t(trX)%*%train1y
weights
#H = trX%*%invsqX%*%t(trX)
#H
summary(trainlm)
tX = testx
tX$newcol <- rep(1, nrow(tX))
tX = tX[,c(9, 1, 2, 3, 4, 5, 6, 7, 8)]
tX = as.matrix(tX)
tX
pred = tX%*%weights
pred
Errors = (testy - pred)^2 
MSE = colMeans(Errors[1])
MSE
proc.time() - ptm

#Gaussian Process 

#Return Covariance matrix C with two cases of covariance function default set 
#linear (i.e. gamma="a")
get_C <- function(M1, M2, rho=NA, gamma="a") {
	C = matrix(0, nrow(M1), nrow(M1))
	n = nrow(M1)
	if (gamma!="a") {
	for (i in 1:n) { 
		for (j in 1:n) { 
 			K = get_hypeK(M1, M2, i, j, rho, gamma)
			C[i,j] = as.numeric(K) }}
	C <- C + diag(nrow(C))}
	else {
	for (i in 1:n) { 
		for (j in 1:n) { 
 			K = get_k(M1, M2, i, j)
			C[i,j] = as.numeric(K) }}
	C <- C + diag(nrow(C))}
	return(C) }

#Linear Covariance Function 
get_k <- function(M1, M2, i, j) { 
	K = (100^2)*M1[i,]%*%t(M2)[,j]	
	K <- as.matrix(K)
	return(K) } 

get_weights <- function(M1, M2, targets, rho=NA, gamma="a") { 
	C <- get_C(M1, M2, rho, gamma)	
	invC <- solve(C)
	weights = invC%*%targets
	return(weights) } 

#Return means of predictive distribution for observations in test matrix using
#covariance function K specified by gamma (for convenience) creating an 1xM 
#vector which the jth component is the value of the covariance function for 
#a certain observation in the test set and the jth observation in the training 
#matrix TrainM, predictions are found by multiplying this vetor by the weights 
get_pred <- function(TestM, TrainM, weights, rho=NA, gamma="a") { 
	n = nrow(TrainM)
	K = matrix(0, 1, n)
	pred = matrix(0, nrow(TestM), 1)
	if (gamma=="a") {
	for (j in 1:nrow(TestM)) {
		for (i in 1:n) {
			k = get_k(TestM, TrainM, j, i) 	
			K[i] = k}
		pred[j] = K%*%weights}}
	else { 
		for (j in 1:nrow(TestM)) {
		for (i in 1:n) {
			k = get_hypeK(TestM, TrainM, j, i, rho, gamma) 	
			K[i] = k}
		pred[j] = K%*%weights}}
	return(pred)}

get_error <- function(predictions, Ys) {
	MSE = mean((predictions - Ys)^2)
	return(MSE)}

#Script for GP 
ptm <- proc.time()
testx = as.matrix(testx)
weights = get_weights(train1x, train1x, train1y)
predictions = get_pred(testx, train1x, weights)
get_error(predictions, testy)
proc.time() - ptm

#Gaussian Process with hyperparameters 

#Exponential covariance function
get_hypeK <- function(M1, M2, i, j, rho, gamma) {
	K = 100^2 + gamma^2*exp(-1*(rho^2)*sum((M1[i,]- M2[j,])^2))
	return(K)}

#Return a matrix M which is either the ith block of horizontal components (if 
#is_train is FALSE) and otherwise Return whatever remains of M1 after the 
#ith block is removed (if is_train is TRUE) 
#For separating data for cross-validation, ith block corresponds to the
#ad hoc test set and the remaining corresponds to the ad hoc training set
splitdata <- function(M1, i, is_train, s=10) { 
	n = nrow(M1)
	each = n/s
	matlist <-list() 
	if (i==0) { 
		X1 <- M1[1:each,]
		X2 <- M1[(each+1):n,]}
	else if ((i+1)==s) { 
		X2 = M1[1:(n-each),]
		X1 = M1[(n-each+1):n,] }
	else {
		Before <- M1[1:(each*i),]
		X1 <- M1[(each*i + 1):(each*(i+1)),]
		After <- M1[(each*(i+1)+1):n,]
		X2 = rbind(as.matrix(Before), as.matrix(After))} 
	#matlist[[1]] = matrix(unlist(X1), ncol=ncol(M1), byrow=TRUE)
	#matlist[[2]] = matrix(unlist(X2), ncol=ncol(M1), byrow=TRUE)
	if (is_train==TRUE) {
		M = X2}
	else {M = X1}
	#M = matrix(M, ncol=ncol(M1), byrow=TRUE)
	return(M)}

#Return a list containing the Least Squared Prediction Error of cross 
#validation over all values of hyperparameters in the specified sets 
#if divideSeven is set to true the first and seventh column of the 
#training X matrix Mx are divided by ten 
get_hype<- function(rho_values,gamma_values, Mx, My, divideSeven=FALSE) {
	Error = 0 
	LSE = "q"
	if (divideSeven) {
		Mx[,1] = Mx[,1]*(1/10)
		Mx[,7] = Mx[,7]*(1/10)
		}
	for (rho in rho_values) { 
		for (gamma in gamma_values) { 
			for (i in 0:9) { 
				Testy = splitdata(My, i, FALSE)
				Trainy = splitdata(My, i, TRUE)
				Testx = splitdata(Mx, i, FALSE)
				Trainx = splitdata(Mx, i, TRUE) 
				weights = get_weights(Trainx, Trainx, Trainy, rho, gamma)
				predictions = get_pred(Testx, Trainx, weights, rho, gamma) 
				Error = Error + get_error(predictions, Testy) }
			if (Error < LSE | LSE == "q") { 
				LSE = Error 
				hypes = c(gamma, rho)
				Error = 0 }
			else {
			Error = 0 }}}
	matlist <- list() 
	matlist[[1]] = LSE 
	matlist[[2]] = hypes
	return(matlist)}

#script for GP with hyperparameters
ptm = proc.time()
weights = get_weights(TRX, TRX, TRY, 0.01, 0.1)
predictions = get_pred(TX, TRX, weights, 0.01, 0.1)
error = get_error(predictions, TY)
error*25
gamma_values = seq(0.1, 10, 0.5) 
gamma_values = c(gamma_values, 10) 
rho_values = seq(0.01, 1, 0.05)
rho_values = c(rho_values, 1) 
error_params = get_hype(rho_values, gamma_values, train1x, train1y)

proc.time() - ptm

ptm = proc.time() 
error_params2 = get_hype(rho_values, gamma_values, train1x, train1y, TRUE)
proc.time() - ptm
