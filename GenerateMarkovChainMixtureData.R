######################################################
#	Generate data from a mixture of Markov chains
# and write it, along with the paramaters of the
# chains, onto files.
#
# mrm: Man Uni, 16 July 2019
#
# N.B. The convention here is that, in a transition
# matrix P, P_{ij} is the probability of a 
# transition i -> j.
######################################################

######################################################
#	Construct, at random, a mixture of Markov
# chains. Most of the randomly-generated parameters
# are drawn from suitable Dirichlet distributions.
######################################################
library( gtools )  # for rdirichlet()

random.Markov.chain <- function( n.states, alpha.vec=NULL )
{
	# Use a uniform Dirichlet distribution by default.
	if( is.null(alpha.vec) ) {
		alpha.vec <- rep( 1, n.states )
	}
	
	# Generate transition probabilities and a distribution
	# over initial states.
	my.markov.params <- list(
		transition.probs = rdirichlet( n.states, alpha.vec ),
		initial.state.probs = as.vector( rdirichlet( 1, alpha.vec ) )
	)
	
	return( my.markov.params )
}

random.Markov.mixture <- function( n.components, n.states, component.alpha=NULL, state.alpha=NULL  )
{
	# Use uniform distribs by default.
	if( is.null(component.alpha) ) { component.alpha <- rep( 1, n.components ) }
	if( is.null(state.alpha) ) { state.alpha <- rep( 1, n.states ) }
	
	# Get the mixture weights
	my.pi <- as.vector( rdirichlet( 1, component.alpha ) ) ;
	
	# Then construct the chains
	my.chains <- lapply( 1:n.components, 
		function(n){ random.Markov.chain(n.states, state.alpha) ; } 
	) ;
	
	my.Markov.mix <- list( 
		n.components = n.components, 
		n.states = n.states,
		pi = my.pi, 
		chains = my.chains 	
	)
	
	return( my.Markov.mix )
}

######################################################
#	Generate simulated data.
######################################################

sample.one.chain <- function( n.transitions, start.probs, trans.probs ) {	
	# Initialise the result
	state.seq <- rep( 0, n.transitions+1 ) ;
	
	# Choose the initial state
	n.states <- length( start.probs )
	state.seq[1] <- sample( n.states, 1, prob=start.probs ) ;
	
	# Now sample the rest in a Markovian way
	pos <- 1 ;
	while( pos < length(state.seq) ) {
		crnt.state <- state.seq[pos] ;
		next.state <- sample( n.states, 1, prob=trans.probs[crnt.state,] )
		state.seq[pos + 1] <- next.state ;
		pos <- pos + 1 ;
	}
	
	return( state.seq )
}

Markov.chain.mixture.data <- function( n.samples, rMm, target.length=NULL ) 
{
	# Unpack some useful info about the mixture
	n.states <- rMm$n.states ;
	n.comps <- rMm$n.components ;
	
	# If need be, set target.length, then get params for the negative binomial
	# distribution from which we'll draw the lengths.
	if( is.null(target.length) ) {
		target.length <- 10 * n.states * n.states ;
	}
	
	# Choose the desired numbers of observed transitions for the sequences
	# by drawing them from a negative binomial distrib. This is an 
	# arbitrary choice and one could use many other approaches
	nb.size <- 2 ; # Target number of successes for neg. binomial
	n.transitions <- rnbinom( n.samples, nb.size, mu=target.length )
	
	# Choose which components of the mixture to sample from
	my.comp.nums <- sample( n.comps, size=n.samples, replace=TRUE, prob=rMm$pi )
	
	# Finally, generate sequences of observed states 
	my.data <- list( 1:n.samples ) # Initialise the result
	for( j in 1:n.samples ) {
		# Get the parameters of the chain from which this sample
		# will be drawn.
		crnt.comp.num <- my.comp.nums[j] ;
		crnt.chain <- rMm$chains[[crnt.comp.num]] ;
		crnt.start.probs <- crnt.chain$initial.state.probs ;
		crnt.trans.probs <- crnt.chain$transition.probs ;
		crnt.n.trans <- n.transitions[j] ;
		
		# Invoke a separate function to do the actual sampling
		my.data[[j]] <- sample.one.chain(
			crnt.n.trans, crnt.start.probs, crnt.trans.probs
		) ;
	}
	
	# Add some convenient data, including a list to the
	# components from which the sequences were sampled, to the
	# result.
	result = list( 
		n.comps = n.comps,
		n.states = n.states,
		true.component = my.comp.nums,
		state.seqs = my.data
	) ;
	
	return( result )
}

# If the punter wants, do some tests.
if( exists("do.tests") && do.tests ) {
	# Do a small example by way of illustration
	print( "Small example with 2 components and 3 states." )
	rMm <- random.Markov.mixture( 2, 3, c(2,2), c(2,2,2) )
	Markov.chain.mixture.data( 10, rMm )

	# Do a bigger example and write everything to files
	n.comps <- 3 ;
	n.states <- 2 ;
	rMm <- random.Markov.mixture( n.comps, n.states, 
		component.alpha=rep(2, n.comps), state.alpha=rep(2, n.states) )
	saveRDS( rMm, file="MarkovMixtureParams.RData" )

	n.seqs <- 100 ;
	my.data <- Markov.chain.mixture.data( n.seqs, rMm )
	saveRDS( my.data, file="MarkovMixtureData.RData" )
}
