################################################################################
# Do EDA on some data about gait reported in 
#
# 	Nathaniel E. Helwig, K. A. Shorter, Ping Ma, E. Hsiao-Wecksler (2016),
# 	Smoothing spline analysis of variance models: A new tool for the analysis 
# 	of cyclic biomechanical data, Journal of Biomechanics, 49(14):3216-3222.
# 	DOI: 10.1016/j.jbiomech.2016.07.035
#
# The data, along with some explanatory amterial about what the levels of some
# of the categorical variables are, is available at
#
#	https://archive.ics.uci.edu/dataset/760/multivariate+gait+data
#
# mrm, Whalley Range,	25 july 2025
################################################################################
rm( list=ls()) # Wipe the slate clean

library(purrr) # list_flatten()
library( RColorBrewer)
source( "MarkovChainMixtureEM.R" ) # EM algorithm stuff

set.seed(123)

pair.pal <- brewer.pal( 12, "Paired" )
dark.pal <- brewer.pal( 8, "Dark2" )

# Read the data, converting some of the variables to factors 
# and giving them meaningful levels
all.data <- read.csv("gait.csv", header=TRUE)
all.data$subject <- factor( all.data$subject )

condition.levels <- c( "unbraced", "knee.brace", "ankle.brace" )
all.data$condition <- factor( condition.levels[all.data$condition] )

leg.levels <- c( "left", "right" )
all.data$leg <- factor( leg.levels[all.data$leg] )

joint.levels <- c( "ankle", "knee", "hip" )
all.data$joint <- factor( joint.levels[all.data$joint] )

summary( all.data )

################################################################################
# Define a function to plot traces for a list of subjects and a particular
# leg, joint and condition
################################################################################

plot.all.reps <-function( subj.nums, cond, lg, jnt, brks="FD" )
{
  # Get all reps for the given combination of categorical variables
  my.df <- subset( all.data, subject%in%subj.nums & condition==cond & leg==lg & joint==jnt )
  
  # Plot a set of empty axes big enough to show all the data
  title.str <- sprintf( "%s, %s leg, %s joint", cond, lg, jnt )
  plot( my.df$time, my.df$angle, type="n",
        xlab="Time (as % of gait cycle)",
        ylab="Angle (degrees)",
        main=title.str
  )
  
  # Now plot all traces with replicates for a given subject shown in the same colour
  n.subjects <- length(subj.nums) 
  for( i in 1:n.subjects ) {
    subj <- subj.nums[i]
    my.df <- subset( all.data, subject==subj & condition==cond & leg==lg & joint==jnt )
    n.reps <- max( my.df$replication)
    for( j in 1:n.reps) {
      tmp.df <- subset( my.df, replication==j)
      lines( tmp.df$time, tmp.df$angle, col=dark.pal[i%%8 + 1] )
    }
  }
  
  legend( "bottomleft", legend=subj.nums, col=dark.pal, lty="solid" )
  
  # Dump the plot to a file
  trace.plot <- recordPlot()
  pdf.path <- sprintf( "Figures/%s_%s_%s_traces.pdf", cond, lg, jnt )
  pdf( file=pdf.path, height=210/25.4, width=297/25.4 )  
  replayPlot(trace.plot)
  dev.off()
  
  # Finally, make a histogram
  my.df <- subset( all.data, subject%in%subj.nums & condition==cond & leg==lg & joint==jnt )
  hist.obj <- hist(my.df$angle, breaks=brks,
                   xlab="Angle (degrees)",
                   main=title.str,
                   col=pair.pal[1],
                   border=pair.pal[2]		
  )
  
  # Save the histogram to a file too
  hist.plot <- recordPlot()
  pdf.path <- sprintf( "Figures/%s_%s_%s_histogram.pdf", cond, lg, jnt )
  pdf( file=pdf.path, height=210/25.4, width=297/25.4 )  
  replayPlot(hist.plot)
  dev.off()
  
  return( hist.obj )
}

################################################################################
# Define a function that accepts a set of bin-boundaries and then digitise 
# traces for a particular collection of subjects and a given leg, joint and 
# condition.
################################################################################

digitise.all.reps <- function( subj.nums, cond, lg, jnt, bb )
{
  # Get all reps for the given combination of categorical variables
  my.df <- subset( all.data, subject%in%subj.nums & condition==cond & leg==lg & joint==jnt )
  
  n.subjects <- length(subj.nums) 
  result <- vector(mode='list', length=n.subjects)
  for( i in 1:n.subjects ) {
    subj <- subj.nums[i]
    my.df <- subset( all.data, subject==subj & condition==cond & leg==lg & joint==jnt )
    n.reps <- max( my.df$replication)
    
    subj.result <- vector(mode='list', length=n.reps)
    for( j in 1:n.reps) {
      tmp.df <- subset( my.df, replication==j)
      subj.result[[j]] <- sapply( tmp.df$angle, function(x){findInterval(x,bb)} )
    }
    
    result[[i]] <- subj.result
  }
  
  return(result)
}

################################################################################
# For every combination of condition, side and joint, fit a mixture of Markov
# chains and record the goodness of fit  stats, as well as the optimal number
# of components.
################################################################################

# Get all combinations of factors
subset.df <- with( all.data, expand.grid(levels(condition), levels(leg), levels(joint)))
names(subset.df) <- c( "condition", "leg", "joint" )

for( j in 1:nrow(subset.df) ) {
  cond <- subset.df$condition[j]
  lg <- subset.df$leg[j]
  jnt <- subset.df$joint[j]
  
  # Draw traces and make a histogram
  hist.obj <- plot.all.reps( 1:10, cond, lg, jnt )
  
  # Digitise the trajectories
  state.traj.lists <- digitise.all.reps( 1:10, cond, lg, jnt, hist.obj$breaks )
  
  aic.vals <- rep(NA, 10)
  EM.result <- vector(mode='list',  length=10)
  for( n.classes in 1:10 ) {
    # Invoke the EM algorithm and see how we do.
    EM.result[[n.classes]] <- MarkovChainMixtureEM( 
      list_flatten(state.traj.lists), 
      n.classes = n.classes, 
      n.states = length(hist.obj$counts),
      max.cycles = 100,
      tol = 1.0e-8
    )
    
    aic.vals[n.classes] <- EM.result[[n.classes]]$AIC 
  }
  
  # Report where the min AIC fell and also seek local minima
  pos.min.aic <- which.min( aic.vals )
  local.minima <- c()
  for( k in 2:9 ) {
    if( (aic.vals[k-1] > aic.vals[k]) && (aic.vals[k+1] > aic.vals[k]) ) {
      local.minima <- c( local.minima, k)
    }
  }
  
  msg <- sprintf( "%s, %s leg, %s joint: min AIC for k=%d", cond, lg, jnt, pos.min.aic )
  if( length(local.minima) == 1 ) {
    msg <- sprintf( "%s, local minimum for k=%d.", msg, local.minima[1] )
  } else if( length(local.minima) > 1 ) {
    tmp.str <- paste( local.minima, collapse=", " )
    msg <- sprintf( "%s, local minima for k in {%s}.", msg, tmp.str )
  }
  
  print( msg )
  
  # Plot an AIC curve
  title.str <- sprintf( "%s, %s leg, %s joint", cond, lg, jnt )
  plot( 1:10, aic.vals,type="n", 
        xlab="Number of classes",
        ylab="AIC",
        main=title.str
  )
  
  pair.num <- 1 + j%%4
  lines( 1:10, aic.vals, lty="solid", col=pair.pal[2*pair.num - 1])
  points( 1:10, aic.vals, pch=19, col=pair.pal[2*pair.num] )
  
  # Save the histogram to a file too
  aic.plot <- recordPlot()
  pdf.path <- sprintf( "Figures/%s_%s_%s_AIC.pdf", cond, lg, jnt )
  pdf( file=pdf.path, height=210/25.4, width=297/25.4 )  
  replayPlot(aic.plot)
  dev.off()
}

##########################################################
#	Look at a case with a promising local minimum.
# Here we build a table in which T_{jk} counts the number
# of subject j's traces that were assigned to class k. 
# We might hope that, for the most part, all 10 end up
# in the same class, but they don't always.
##########################################################

cond <- "knee.brace"
lg <- "right"
jny <- "hip"
hist.obj <- plot.all.reps( 1:10, cond, lg, jnt )
state.traj.lists <- digitise.all.reps( 1:10, cond, lg, jnt, hist.obj$breaks )
EM.result <- MarkovChainMixtureEM( 
  list_flatten(state.traj.lists), 
  n.classes = 6, 
  n.states = length(hist.obj$counts),
  max.cycles = 100,
  tol = 1.0e-8
)

assignment.mat <- matrix( rep(0, 60), nrow=10 )
assignment <- apply( EM.result$gamma.mat, MARGIN=1, FUN=which.max)
for( i in 1:10 ) {
  min.idx <- 1 + 10*(i-1)
  max.idx <- 10*i
  assignment.mat[i,] <- tabulate( assignment[min.idx:max.idx], nbins=6 )
}

assignment.mat
################################################################################

library(ggplot2)
library(dplyr)


plot_data <- all.data %>%
  filter(condition == "knee.brace", leg == "right", joint == "hip")

# 2. Create a unique ID for these 100 trajectories and sort them accordingly.
plot_data <- plot_data %>%
  arrange(subject, replication) %>%
  mutate(trajectory_id = rep(1:100, each = 101))

# 3. Add clustering results (assignment) to the data
cluster_assignments <- rep(assignment, each = 101)
plot_data$cluster <- factor(cluster_assignments) 

# 4. use ggplot2 print diagram
cluster_plot <- ggplot(plot_data, aes(x = time, y = angle, group = trajectory_id, color = cluster)) +
  geom_line(alpha = 0.6) + 
  labs(
    title = "Gait Traces by Identified Cluster (k=6)",
    subtitle = "Condition: Knee Brace, Leg: Right, Joint: Hip",
    x = "Time (% of Gait Cycle)",
    y = "Angle (degrees)",
    color = "Cluster"
  ) +
  theme_minimal() + 
  scale_color_brewer(palette = "Set2") 


print(cluster_plot)


ggsave("Figures/knee_brace_right_hip_by_cluster.pdf", plot = cluster_plot, width = 11, height = 8)

cat("\nTask A completed: Trajectory diagram coloured by cluster category generated and saved\n")
################################################################################


library(tidyr)

# 1. Extract the transfer matrix list from the EM results.
transition_matrices <- EM.result$trans.mats

# 2. Get status count
num_states <- EM.result$n.states

# 3. Iterate through the six categories and create a heat map for each category.
for (i in 1:length(transition_matrices)) {
  
  # Extract the i-th matrix
  T_matrix <- transition_matrices[[i]]
  
  # Convert the matrix into a long-format data frame for easy ggplot plotting.
  T_df <- as.data.frame(T_matrix)
  colnames(T_df) <- 1:num_states
  T_df$From <- 1:num_states
  
  T_df_long <- T_df %>%
    pivot_longer(cols = -From, names_to = "To", values_to = "Probability") %>%
    mutate(To = as.integer(To))
  
  # create heatmap
  heatmap_plot <- ggplot(T_df_long, aes(x = To, y = From, fill = Probability)) +
    geom_tile(color = "white") + # geom_tile
    geom_text(aes(label = round(Probability, 2)), color = "black", size = 3) + # 在格子上显示概率值
    scale_fill_gradient(low = "white", high = "steelblue", name = "Probability") + # 设置颜色梯度
    scale_y_reverse() + 
    coord_fixed() + 
    labs(
      title = paste("Transition Matrix for Cluster", i),
      x = "To State",
      y = "From State"
    ) +
    theme_minimal() +
    theme(axis.ticks = element_blank())
  
  
  print(heatmap_plot)
  
  # ensure create "Figure' files
  ggsave(paste0("Figures/transition_matrix_cluster_", i, ".pdf"), plot = heatmap_plot, width = 8, height = 7)
}

cat("\Task B completed: Heat maps for six transfer matrices have been generated and saved.\n")

