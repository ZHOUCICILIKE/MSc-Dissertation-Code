####################################################################
#	A circle of routines that implement the EM algorithm 
# for mixtures of Markov chains.
# 
# mrm: Man Uni, 16 July 2019
# mrm & Rajenki Das: Man Uni, 5 Feb. 2020
#
# N.B. The convention here is that, in a transition
# matrix P, P_{ij} is the probability of a 
# transition i -> j.
######################################################
library( gtools )  # for rdirichlet()

####################################################################
#	Given a trajectory - a sequence of states drawn from
# {1, .., n.states} - build a matrix C whose i,j entry counts
# the number of times that the trajectory includes a transition 
# from state i to state j.
####################################################################

count.transitions <- function( traj, n.states=NA )
{
	# Deal with the possibility that n.states wasn't given 
	if( is.na(n.states) ) { n.states <- max( traj ) }
	
	# Initialise the result as a matrix full of zeroes.
	count.mat <- matrix( rep(0, n.states*n.states), nrow=n.states )
	
	# Run through the trajectory counting transitions
	n.transitions <- length( traj ) - 1
	for( t in 1:n.transitions ) {
		from.state <- traj[t] ;
		to.state <- traj[t+1] ;
		count.mat[from.state, to.state] <- count.mat[from.state, to.state] + 1 ;
	}
		
	return( count.mat )
}

if( exists("do.tests") && do.tests ) {
	# Do a small test
	test.traj <- c( 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2, 1, 2 ) ;
	count.transitions( test.traj, 2 )
	count.transitions( test.traj, 3 )
	count.transitions( test.traj )
}

####################################################################
#	Given a count matrix of the sort produced by count.transitions()
# and a matrix of transition probabilities, compute a log-likelihood.
####################################################################

# It'll prove helpful to have a separate test for the 
# presence of impossible transitions.
traj.impossibleQ <- function( count.mat, trans.mat ) {
	count.vec <- as.vector( count.mat )
	prob.vec <- as.vector( trans.mat ) ;

	return( any( (prob.vec == 0) & (count.vec != 0) ) ) ;
}

traj.loglike <- function( count.mat, trans.mat ) {
	count.vec <- as.vector( count.mat )
	prob.vec <- as.vector( trans.mat ) ;
	log.prob.vec <- ifelse( prob.vec==0, 0, log(prob.vec) ) ;
	loglike <- sum(count.vec * log.prob.vec)
	
	return( loglike  ) ;
}

if( exists("do.tests") && do.tests ) {
	# Another test
	small.trans.mat <- matrix( rep( 0.5, 4), nrow=2 ) ;
	test.loglike <- traj.loglike( count.transitions( test.traj, 2 ), small.trans.mat )
	test.loglike == log(0.5) * (length(test.traj) - 1)

	small.trans.mat <- matrix( c(1, 0, 0, 1), nrow=2 ) ;
	traj.impossibleQ( count.transitions( test.traj, 2 ), small.trans.mat )
}

####################################################################
#	Given a matrix gamma whose i,j entry is the probability
# that the i-th trajectory belongs to the j-th class, compute a 
# vector pi whose j-th entry is the probability that a 
# randomly-selected trajectory belongs to class j.
####################################################################

estimate.mixture.weights <- function( gamma.mat )
{
	n.trajectories <- nrow( gamma.mat )
	pi.vec <- colSums( gamma.mat ) / n.trajectories
	return( pi.vec )
}

if( exists("do.tests") && do.tests ) {
	# As ever, do a test
	test.gamma.entries <- c( 1.0, 0.0, 1.0, 0.0, 0.0, 1.0 )
	test.gamma.mat <- matrix( test.gamma.entries, byrow=TRUE, ncol=2 )
	test.gamma.mat
	estimate.mixture.weights( test.gamma.mat )
}

####################################################################
#	Given a matrix gamma as described above and a list
# of count matrices of the sort produced by count.transitions(),
# estimate the transition matrices for all the classes.
####################################################################
# 文件：MarkovChainMixtureEM.R

estimate.transition.matrices <- function( gamma.mat, count.mats )
{
  # Examine gamma.mat and count.mat.list to get 
  # the numbers of trajectories, classes and states
  n.trajectories <- nrow( gamma.mat )
  n.classes <- ncol( gamma.mat )
  n.states <- nrow( count.mats[[1]] )
  
  # Initialise a list of empty pseudocount matrices, one per class.
  empty.mat <- matrix( rep( 0.0, n.states*n.states), nrow=n.states )
  pseudocount.mats <- lapply( 1:n.classes, function(j){ empty.mat }) ;
  
  # Accumulate the pseudocount matrices, adding a fraction gamma_{kj} 
  # of trajectory k's counts to the pseudocount matrix for the j-th class.
  for( k in 1:n.trajectories ) {
    for( j in 1:n.classes ) {
      pseudocount.mats[[j]] <- pseudocount.mats[[j]] + gamma.mat[k,j] * count.mats[[k]] ;
    }
  }
  
  # Normalize the pseudocount matrices to get transition matrices.
  # Given that the i,j entry counts transitions from i to j, we want
  # the *rows* of our transition matrices to sum to one.
  for( j in 1:n.classes ) {
    
    # --- 新增的修正逻辑 ---
    # 检查是否有任何行的和为零（即，从未作为起始状态的状态）
    crnt.row.sums <- rowSums( pseudocount.mats[[j]] )
    zero.sum.rows <- which( crnt.row.sums == 0 )
    
    if( length(zero.sum.rows) > 0 ) {
      # 对于这些行，我们假设它以概率1转移到自身
      for( row_idx in zero.sum.rows ) {
        pseudocount.mats[[j]][row_idx, row_idx] <- 1.0
      }
    }
    # --- 修正逻辑结束 ---
    
    # 现在，重新计算行和并进行归一化
    # 所有行的和现在都应该大于0
    crnt.row.sums <- rowSums( pseudocount.mats[[j]] )
    
    # 创建归一化矩阵 D，其中 D_{ii} = 1 / (第i行和)
    # 旧代码中对0值的处理现在不再是必需的，但保留它以增加稳健性
    diag.entries <- ifelse( crnt.row.sums==0, 1, 1/crnt.row.sums )
    
    # %*% 是矩阵乘法
    pseudocount.mats[[j]] <- diag( diag.entries ) %*% pseudocount.mats[[j]]
  }
  
  return( pseudocount.mats )
}

if( exists("do.tests") && do.tests ) {
	# As ever, do a test
	test.count.mat <- count.transitions( test.traj, 2 )
	test.count.mat
	test.count.mats = list( test.count.mat, test.count.mat, t(test.count.mat))
	estimated.trans.mats <- estimate.transition.matrices( test.gamma.mat, test.count.mats )
	estimated.trans.mats
}

####################################################################
#	Given a Markov-process mixture model and some count 
# matrices of the sort produced by count.transitions(),
# compute gamma, the matrix of class membership probs.
####################################################################

compute.gamma.mat <- function( count.mats, pi.vec, trans.mats )
{
	# Examine pi.vec and count.mats to get 
	# the numbers of trajectories and classes.
	n.trajectories <- length( count.mats )
	n.classes <- length( pi.vec )
	
	# Initialise gamma.mat with NA's. That way, if we somehow fail 
	# to fill some entry in, the mistake will be evident. 
	gamma.mat <- matrix( rep(NA, n.trajectories*n.classes), nrow=n.trajectories )
	
	# Arrange for the k,j entry to be the likelihood of the k-th count matrix 
	# being generated from the j-th transition matrix
	for( k in 1:n.trajectories ) {
		for( j in 1:n.classes ) {
			if( traj.impossibleQ( count.mats[[k]], trans.mats[[j]] ) ) {
				gamma.mat[k,j] = 0.0 ;
			} else {
				crnt.loglike <- traj.loglike( count.mats[[k]], trans.mats[[j]] ) ;
				gamma.mat[k,j] <- pi.vec[j] * exp( crnt.loglike ) ;
			}
		}
	}
	
	# Normalize the matrix we've just constructed to get probabilities of membership.
	gamma.row.sums <- rowSums( gamma.mat ) ;
	if( any(gamma.row.sums <= 0) ) {
		warning( "Row sum of zero?!?" ) ;
	}
	
	diag.entries <- 1/gamma.row.sums ;
	gamma.mat <- diag( diag.entries ) %*% gamma.mat ;
	return( gamma.mat ) ;
}

if( exists("do.tests") && do.tests ) {
	compute.gamma.mat( test.count.mats, c(2/3, 1/3), estimated.trans.mats )

	#### Do some speed tests ##########
	library( microbenchmark )

	microbenchmark(
		count.transitions( test.traj, 2 ),
		estimate.mixture.weights( test.gamma.mat ),
		estimate.transition.matrices( test.gamma.mat, test.count.mats ),
		compute.gamma.mat( test.count.mats, c(2/3, 1/3), estimated.trans.mats ) 
	)
}

#####################################################################
#	Generate a random gamma matrix to initialise the EM algorithm
#####################################################################

random.gamma.mat <- function( n.trajectories, n.classes, alpha.val=0.5 ) {
	# Draw the values from a Dirichlet distribution
	alpha <- rep( alpha.val, n.classes )
	
	return( rdirichlet( n.trajectories, alpha ) )
}

if( exists("do.tests") && do.tests ) {
	random.gamma.mat( 4, 3 )
}

#####################################################################
#	Given a collection of count matrices and a fitted mixture
# of Markov chains, compute the log-likelihood, AIC and BIC
#####################################################################

log.sum.exp <- function( x ) {
	# Use a standard trick to compute the log of a sum of exponentiated terms
	# See, e.g., https://en.wikipedia.org/wiki/LogSumExp
	x.star <- max( x )
	return( x.star + log(sum(exp(x - x.star))) )
}

goodness.of.fit <- function( mc.mix, count.mats ) {
	# Begin by computing the log-likelihoods
	loglike <- 0.0 
	for( j in 1:mc.mix$n.trajectories ) {
		# First compute the log-likelihoods conditioned on class membership
		loglike.vec <- rep( 0.0, mc.mix$n.classes )
		for( k in 1:mc.mix$n.classes ) {
			loglike.vec[k] <- traj.loglike( count.mats[[j]], mc.mix$trans.mats[[k]] ) ;
			loglike.vec[k] <- loglike.vec[k] + log(mc.mix$pi[k])
		}
		
		loglike <- loglike + log.sum.exp( loglike.vec )
	}
	
	# Compute AIC and BIC
	n.transitions <- sum(sapply(unlist(count.mats), sum))
	n.params <- (mc.mix$n.states*(mc.mix$n.states - 1))*mc.mix$n.classes # Trans mats
	n.params <- n.params + (mc.mix$n.classes - 1) # class membership probs
	my.AIC <- 2.0*(n.params - loglike)
	my.BIC <- n.params*n.transitions -2.0*loglike # Is a single transition an observation?
	
	# Assemble the results and return
	result <- list(
		loglike	= loglike,
		AIC		= my.AIC,
		BIC		= my.BIC
	)
	
	return( result )
}

#####################################################################
#	Given a list of state trajectories and a number of classes,
# estimate a mixture of Markov chains that could have generated them.
#####################################################################

MarkovChainMixtureEM <- function(
	traj.list, n.classes, n.states=NA, max.cycles=100, tol=0.000001 
) {
	# If need be, get the number of states.
	if( is.na(n.states) ) {
		max.per.traj <- lapply( traj.list, max ) ;
		n.states <- max( max.per.traj ) ;
	}
	
	# Reduce the trajectories to count matrices
	count.mats <- lapply( traj.list, 
		function( my.mat ) { count.transitions( my.mat, n.states ) } 
	) ;
	
	# Generate a random matrix of class membership probs. This
	# serves a starting guess for gamma.
	n.trajectories <- length( traj.list ) ;
	gamma.mat <- random.gamma.mat( n.trajectories, n.classes ) ;
	
	# Do the thing.
	n.cycles <- 0 ;
	rms.delta.gamma <- NaN ; # This should get a numerical value before we need it.
	while( 
		(n.cycles == 0) ||
		((rms.delta.gamma >= tol) && (n.cycles < max.cycles)) 
	) {
		# Compute maximum-likelihood mixture parameters based on gamma
		pi.vec <- estimate.mixture.weights( gamma.mat ) ;
		trans.mats <- estimate.transition.matrices( gamma.mat, count.mats ) ;
		
		# Estimate a new gamma matrix based on the new parameters
		prev.gamma <- gamma.mat ;
		gamma.mat <- compute.gamma.mat( count.mats, pi.vec, trans.mats ) ;
		
		# Note that we've completed another EM cycle.
		n.cycles <- n.cycles + 1 ;
		
		# Evaluate the convergence criterion
		delta.gamma <- gamma.mat - prev.gamma ;
		dg.squared <- delta.gamma * delta.gamma ; # elementwise multiplication
		rms.delta.gamma <- sqrt( sum(dg.squared)/ (n.trajectories*n.classes) ) ;
	}
		
	# The result has a degeneracy induced by permuting the class labels and so,
	# following advice offered by Michael Betancourt at
	#
	#	https://betanalpha.github.io/assets/case_studies/identifying_mixture_models.html
	#
	# we lift this degeneracy by enforcing the condition that pi.vec is increasing.
	perm <- order( pi.vec ) ;
	
	# Assemble the result
	result <- list(
		# The input data
		n.trajectories	= n.trajectories,
		n.classes		= n.classes,
		n.states		= n.states,
		traj.list 		= traj.list,
		
		# The results
		gamma.mat	= gamma.mat[,perm],
		trans.mats	= trans.mats[perm],
		pi			= pi.vec[perm],
		
		# Details about the computation
		n.cycles		= n.cycles,
		converged		= (rms.delta.gamma < tol),
		rms.delta.gamma = rms.delta.gamma
	)
	
	# Add some goodness-of-fit measures
	gof.result <- goodness.of.fit( result, count.mats )
	result$loglike <- gof.result$loglike
	result$AIC <- gof.result$AIC
	result$BIC <- gof.result$BIC

	return( result )
}
