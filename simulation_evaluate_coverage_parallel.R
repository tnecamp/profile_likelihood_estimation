## Double check carefulness about the seed
## Things are done in the setting where I know the value of the likelihood of the MLE but I don't observe the MLE directly
## I also only do estimation for one set of initial estimates, I could change how I initialize my quadratic, 
## and how many initial points I use

require(doParallel)
registerDoParallel(16);
library(MASS)
library(mvtnorm)

## These will be updated additively, could have stored as a vector and taken mean, but no need


### Define optimization function to get estimate of profile likelihood ###

## This is where I define the meta model optimization function, it finds the quadratice with the highest
## likelihood for a give meta-model
## a_init,b_init,c_reparam define the initial quadratic guess
## x_star and y_star are the data values
## y_star_max is the maximum of all y_star, helpful when parameterizing c, bu could remove
## mcmc_sample_size is the number of monte carlo points you want to use to approximate the distribution
## sample_size is how many data points I have to evaluate the likelihood at
## This function works for any set of x_star, y_star but is particular to the specified distributions
meta_model_optimization = function(a_init, b_init, c_reparam, x_star, y_star, x_star_sd, y_star_max, mcmc_sample_size, sample_size){

	## Define likelihood function to optimize
	full_log_likelihood_alpha = function(a, b, c_reparam){
		## Set seed every time to make the mcmc sample consistent, i.e. for the same parameters, I will get the same likelihood
		## estimate
		set.seed(100)
		mcmc_sample = rnorm(mcmc_sample_size, mean = -b/(2*a), sd = sqrt(-1/(2*a)))
		temp  = a*mcmc_sample^2 + b*mcmc_sample + c_reparam + b^2/(4*a) + y_star_max
		cur_sum = 0
		  
		## Could I vectorize this?
		for(j in 1:sample_size){
			like_y_star = dgamma(temp - y_star[j], shape = .5, rate = 1)
			like_x_star = dnorm(mcmc_sample-x_star[j], mean = -.067, sd = x_star_sd)
			cur_sum = cur_sum + log(1/mcmc_sample_size*sum(like_y_star*like_x_star))
		}
		return(cur_sum)
	}

	## Within this I reparameterize to allow the optimization to not have any constraints, also prevents 
	## optimization errors
	optim_fun2 = function(arg_vec){
		a = -exp(arg_vec[1])
		b = arg_vec[2]
		return(-full_log_likelihood_alpha(a, b, c_reparam))
	}
	
	optim_sol = optim(c(log(-a_init), b_init), optim_fun2, hessian = TRUE)
	return(list(optim_sol$par[1], optim_sol$par[2], optim_sol$hessian))

}
                                 	

### Begin Simulation ###

	
                                      
total_sim = function(k, mu, sigma, n, t_g, sample){
	
	
	### Create data set for this iteration and true parameter values ###
	
	## Simulation parameters

	sigma_inv = solve(sigma)

	## Create data
	data = mvrnorm(n, mu, sigma)
	data_mean = colMeans(data)
	
	## Get the true profile likelihood function, need to redefine this everytime I get a new data set
	pl_eval = function(theta){
		## mu1 = theta case
		mu_2_max = data_mean[2] + sigma_inv[1,2]/sigma_inv[2,2]*(data_mean[1]-theta)
		mu_2_max = min(theta, mu_2_max)
		like_1 = dmvnorm(c(theta,mu_2_max), data_mean, sigma/n, log  = TRUE)
		  
		## mu2 = theta case
		mu_1_max = data_mean[1] + sigma_inv[1,2]/sigma_inv[1,1]*(data_mean[2]-theta)
		mu_1_max = min(theta, mu_1_max)
		like_2 = dmvnorm(c(mu_1_max,theta), data_mean, sigma/n, log  = TRUE)
		return(max(like_1,like_2))
	}
	
	## Getting MLE of my underlying data set
	mle_sample = max(data_mean)
	true_MLE = mle_sample
	
	## Get the likelihood value of MLE
	mle_like = dmvnorm(data_mean, data_mean, sigma/n, log  = TRUE)
	pl_eval(mle_sample) == mle_like  # should be true since the profile likelihood equals the true likelihood for the MLE
	
	## Plotting true profile likelihood
	theta_set = max(data_mean) + 0:5000/1000 - 2.5
	theta_set_like = sapply(theta_set, FUN = pl_eval)
	pl_cutoff = mle_like - 1.920729
	
	# true confidence interval lower and upper bound, i.e. if I knew the true profile likelihood
	true_L_cur = theta_set[(theta_set_like > pl_cutoff)][1] -.0005	# lower bound
	true_U_cur = as.numeric(tail(theta_set[(theta_set_like > pl_cutoff)], n=1)) + .0005 	  # upper bound
	true_L = true_L_cur   # storage
	true_U = true_U_cur	 # storage
	
	
	### Generate points 'below' the profile likelihood that I will use to estimate true profile likelihood ###

	
	mu_hat_max = vector(length = sample)	  ## Allocate the memory for storing our noisy x_star values
	
	## Get my sample of points to use for estimating the profile likelihood i.e. points below the PL
	mu_sample = mvrnorm(sample, data_mean, sigma/n)		# Generate points from a normal, this distribution can be specified, for example, sigma does not have to be the same, and it doesn't have to be normal
	mu_sample_eval = apply(mu_sample, MARGIN = 1, FUN = max)  # Get the max, i.e. x-coordinate without horizontal noise
	likehood_sample = dmvnorm(mu_sample, data_mean, sigma/n, log  = TRUE)	# Get their likelihood i.e. y coordinate	
	
	## Add in horizontal noise
	 mu_hat_max = apply(mvrnorm(sample, c(0,0), sigma/t_g) + mu_sample, MARGIN = 1, FUN = max)
	
	## This is to get an understanding of the horizontal error we introduced and how it comes
	## out in the maximum, this is not necessary for the final result, but for exploration purposes
	epsilon_vec = mu_hat_max - mu_sample_eval
	
	y_star_max = max(likehood_sample)	# Get largest likelihood value
	x_star_sd = sd(epsilon_vec)	# Get estimate of horizontal error. Note this should change since it's based on unknown truth 
	bias_g = sum(epsilon_vec)/sample	# get estimate of horizontal bias
	SD_g =  x_star_sd 	# storage
	
	
	### Given my points, I get an estimate of the profile likelihood  ###
	
	## Get initial quadratic guess
	## I can alter this to get better initial estimates
	curvature = -5		# Inital estimate of curvature
	center =  mean(mu_hat_max)		# Initial estimate of center of my quadratic
	height = mle_like		# My height is based on the likelihood of true mle which is known
	
	## Get the corresponding values for a quadratic function
	a_init = curvature
	b_init = -2*curvature*center
	c_reparam = height - y_star_max
	
	## Find the optimized quadratice parameters, i.e. my PL estimate
	optimized_parameters = meta_model_optimization(a_init, b_init, c_reparam, mu_hat_max, 
	likehood_sample, x_star_sd, y_star_max, 1000, sample)
	
	a_reparam = optimized_parameters[[1]]
	a = -exp(a_reparam)
	b = optimized_parameters[[2]]
	information_est = optimized_parameters[[3]]		## keep this positive since I minimized the negative log likelihood
	
	
	### Storage of values and finding new cut offs for our profile likelihood ###
	
	MLE_hat = -b/(2*a)		# Store estimate of MLE based on PL
	grad_mle = c(-b/(2*exp(a_reparam)), 1/(2*exp(a_reparam)))
	cur_inv = solve(information_est)		# get error estimates of parameters based on hessian
	cur_MLE_var = grad_mle%*%cur_inv%*%grad_mle		# get estimate of MLE variance
	MLE_var = cur_MLE_var	# store MLE variance estimate
	a_store = a		# store curvaturue
	b_store = b		# store b value in quadratice
	a_reparam_store = a_reparam
	a_var_store = exp(-2*a_reparam)*cur_inv[1,1]		# variance in a estimate
	b_var_store =  cur_inv[2,2]		# variance in b estimate
	a_reparam_var_store = cur_inv[1,1]		# variance in reparameterized a
	
	## obtain new profile likelihood cutoff based on estimated PL
	new_cut_off = y_star_max - 1.92 - 3.84*(-a)*max(cur_MLE_var,0)  # Tim double check this should be y_star_max vs mle_like
	
	
	L_vec_noise_cur =  -sqrt((new_cut_off - (c_reparam+y_star_max)) / a ) - b/(2*a)		# New estimated lower bound
	U_vec_noise_cur =  sqrt((new_cut_off - (c_reparam+y_star_max)) / a ) - b/(2*a)		# New estimated upper bound
	
	L_vec_noise =  L_vec_noise_cur
	U_vec_noise =  U_vec_noise_cur 
	 
	
	
	### Create a plots to check out the results while it's running ###
	
	## Plot the true profie likelihood in blue
	#plot(theta_set, theta_set_like, cex = .02, xlim = c(mle_sample-2.3,mle_sample+2.3), ylim = c(mle_like-3, mle_like+.3), col=4, xlab = "Theta", ylab = "Likelihood", main = "Sim2-mu_2, small MCMC sample")
	#points(mle_sample,mle_like, col=4, cex = .3)	## plot the mle
	
	## plot cut off and lower and upper bounds for true PL
	#abline(h = pl_cutoff)
	#points(true_L_cur, pl_cutoff, col = 4)
	#points(true_U_cur, pl_cutoff, col = 4)
	
	## Plot the initial estimate of the profile likelihood in turqouise
	#quad_fun = function(x) {return(a_init*(x+b_init/(2*a_init))^2 + c_reparam + y_star_max)}
	#points(theta_set, quad_fun(theta_set),col=5, cex = .02)
	
	## Plot the optimized PL estimate in pink
	#quad_fun = function(x) {return(a*(x+b/(2*a))^2 + c_reparam+y_star_max)}
	#points(theta_set, quad_fun(theta_set),col=6, cex = .02)
	
	## Plot the new cut off and lower and bounds for estimated PL
	#abline(h = new_cut_off)
	#abline(v = max(mu) )
	#points(L_vec_noise_cur, new_cut_off, col = 4)
	#points(U_vec_noise_cur, new_cut_off, col = 4)
	
	## Print the iteration
	#print(k)
	return(list(true_MLE = true_MLE, true_L = true_L, true_U = true_U, 
	MLE_hat = MLE_hat, MLE_var = MLE_var, a_store = a_store, b_store = b_store, 
	a_reparam_store = a_reparam_store, a_var_store = a_var_store, 
	b_var_store = b_var_store, a_reparam_var_store = a_reparam_var_store, 
	L_vec_noise = L_vec_noise, U_vec_noise = U_vec_noise, SD_g = SD_g, 
	bias_g = bias_g))
}

dat <- foreach(i = 1:100) %dopar% {
	## Here I have to reset the seed because, if not, the seed reset in my optimization will mess up 
	## my data
	set.seed(i+1000)
	
   	total_sim(i, mu = c(5.2, 5.15), sigma = 2*matrix(c(1,.5,.5,1), 2, 2), 
   	n = 20 , t_g = 10, sample = 100)
  }
combined_dat = mapply(c,dat)
combined_dat = data.frame(apply(t(combined_dat), unlist, MARGIN = 2))


### Check out results ###

# Compare variance of MLE_hat to estimated variance of MLE_hat
var(combined_dat$MLE_hat-combined_dat$true_MLE)
mean(combined_dat$MLE_var)

## Bias of MLE_hat
mean(combined_dat$MLE_hat-combined_dat$true_MLE)

which(!(combined_dat$true_U > 5.2 & combined_dat$true_L < 5.2))

which(!(combined_dat$U_vec_noise > 5.2 & combined_dat$L_vec_noise < 5.2))

sum(combined_dat$true_U > 5.2 & combined_dat$true_L < 5.2)
sum(combined_dat$U_vec_noise > 5.2 & combined_dat$L_vec_noise < 5.2)


mean(combined_dat$bias_g)


which((combined_dat$true_U < 5.2 & combined_dat$true_L < 5.2))
which((combined_dat$true_U > 5.2 & combined_dat$true_L > 5.2))


which((combined_dat$U_vec_noise < 5.2 & combined_dat$L_vec_noise < 5.2))
which((combined_dat$U_vec_noise > 5.2 & combined_dat$L_vec_noise > 5.2))
