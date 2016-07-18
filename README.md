# GaussianProcessRegression
Regression with Gaussian Processes (with and without hyperparameters made from first principles

Only dependeny is cvTooks so just a simple install.packages("cvTools") should sort you out. 

(I think) A kind of cool Gaussian process regression I coded to try to get my head around how the whole thing worked. 

Allows an option for Gaussian Process with and without hyperparameters (for without just don't specify a fourth parameter). 

There are some much better libraries with this and more floating all over the place but when I was first getting into ML I always found 
it helped me a lot to program them myself. Plus ya never know how much you may want to customize your model. Tough to add in a 
generalized Mahalanobis distance into the prebuilt ones if you ever need to include a cool distance function for categorical variables. 


