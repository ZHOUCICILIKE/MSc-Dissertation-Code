####################################################################
#	Test a circle of functions that use the EM algorithm to fit
# a mixture of Markov Chains to a collection of state trajectories.
#
# mrm & Rajenki Das: Man Uni, 5 Feb. 2020
####################################################################

rm( list=ls() ) # Wipe the slate clean.
source( "MarkovChainMixtureEM.R" ) # EM algorithm stuff
source( "GenerateMarkovChainMixtureData.R" )

# Set up a test problem proposed by Thomas House. It's a mixture 
# with two components, each of which is a two-state Markov chain
# with one transient state and one absorbing state. The expected
# time spent in the transient state is set by a parameter tau
# below.
#
# The conventions used to specify the mixture are the same as those
# used in GenerateMarkovChainMixtureData.R

tau <- 50 # expected residence in transient state
p.stay <- tau/(1 + tau)
p.go <- 1.0 - p.stay
chain.A <- list(
	transition.probs = matrix( c(p.stay, p.go, 0.0, 1.0), byrow=TRUE, nrow=2 ),
	initial.state.probs = c( 0.99, 0.01 )
)

chain.B <- list(
	transition.probs = matrix( c(1.0, 0.0, p.go, p.stay), byrow=TRUE, nrow=2 ),
	initial.state.probs = c( 0.01, 0.99 )
)

my.Markov.mix <- list( 
	n.components	= 2,
	n.states 		= 2,
	pi 				= c( 0.4, 0.6 ), # should be in increasing order
	chains			= list( chain.A, chain.B ) 	
)

my.Markov.mix # display the test problem

# Generate some random trajectories from the mixture specified above
n.subjects <- 50
target.length <- 3*tau 
synthetic.data <- Markov.chain.mixture.data( n.subjects, my.Markov.mix, target.length )
traj.list <- synthetic.data$state.seqs 

# Try the EM algorithm for a range of numbers of components.
min.n.comp <- max( c(1,my.Markov.mix$n.components-2) )
max.n.comp <- my.Markov.mix$n.components + 2
all.EM.results <- vector(mode='list',  length=4)
for( n.classes in min.n.comp:max.n.comp ) {
	# Invoke the EM algorithm and see how we do.
	all.EM.results[[n.classes]] <- MarkovChainMixtureEM( 
		traj.list, 
		n.classes = n.classes, 
		n.states = my.Markov.mix$n.states,
		max.cycles = 100,
		tol = 1.0e-8
	)
	
	print( all.EM.results[[n.classes]]$BIC )
}

# Look at the results for the correct mixture
EM.result <- all.EM.results[[my.Markov.mix$n.components]]

# Make a confusion matrix
assigned.component <- sapply( 1:n.subjects, 
	function(j ) { which.max(EM.result$gamma.mat[j,]) }
)

confusion.df <- data.frame(
	true = synthetic.data$true.component,
	assigned = assigned.component
) ;

xtabs( ~ assigned + true, data=confusion.df ) 

# Make a plot where the points are coloured by their true class,
# and the vertical axis shows prob of being in class 1.
library( RColorBrewer )

plot.df <- data.frame(
	prob.class.one = EM.result$gamma.mat[,1],
	true = synthetic.data$true.component
) ;

perm <- with( plot.df, order(true, prob.class.one))
plot.df <- plot.df[perm,]

dark.pal <- brewer.pal( 8, "Dark2" )
plot( x=1:n.subjects, y=plot.df$prob.class.one, 
	type="p", pch=19, col=dark.pal[plot.df$true],
	xlab = "Subject number", 
	ylab = "Fitted prob. of being in component 1",
	main = paste( 
		"Mixture of ", my.Markov.mix$n.states,
		"-state Markov chains: ", 
		my.Markov.mix$n.components, 
		" components", sep=""
	)
)

abline( h=0.5, lty="dashed", col="gray75" )

n.classes = my.Markov.mix$n.components ;
legend( 
	x=0.5*n.subjects, y=1.0, 
	title = "True component",
	legend=1:n.classes,
	col=dark.pal[1:n.classes], pch=19
)