directory <- '~/program/rpeode/package/R/'
source( paste(directory,'generate.linear.R',sep='') )

#library('peode')

# Random Seed

arg <- commandArgs(TRUE)

if ( length(arg) > 0 )
{
  set.seed(arg[1])
} else
{
  set.seed(0)
}

# Parameter

dimension <- as.integer(5)
time_point <- seq ( 0 , 1 , length.out=10 ) #[-1]

orthogonal_transformation <-
  lapply (
    1 : ((dimension-1)/2) , function(index)
    {
      return ( as.integer ( c ( 2*index , 2*index+1 ) ) )
    } )

# Generate model
object <-
  generate.linear (
    dimension
    , time_point
#    , scaling = TRUE
#    , orthogonal_transformation = orthogonal_transformation
#    , row_column_permutation = TRUE
#    , constant = list ( 0 , 1e-2 )
    , noiselevel = 0.0
    , sanitycheck = TRUE
    , output =
        list (
          file = 'data'
#          file = '~/program/peode/test/linear_inverse_eigen/ConstructODE/data'
#          file = paste('~/tmp/bash/data',arg,sep='')
          , sparse = FALSE
          , timefirst = FALSE
        )
  )

# Check curve with "ode" function from "deSolve" package.
#linODE <- function ( time , state , pars )
#{
#  res <- pars[[1]] %*% state
#  if ( ! is.null ( pars[[2]] ) )
#    res <- res + pars[[2]]
#  return ( list(res) )
#}
#
#library('deSolve')
#add0 <- time_point[1]!=0
#if ( add0 )
#{
#  time_point <- c ( 0 , time_point )
#}
#ode_res <- deSolve::ode (
#  object$initial
#  , time_point
#  , linODE
#  , list ( object$linear , object$constant )
#)
#if ( add0 )
#{
#  ode_res <- ode_res[-1,]
#}
#
#difference <- ode_res[,-1] - object$data[,-1]
#rss <- norm ( difference , 'F' )
#cat('rss: ')
#cat(rss)
#cat('\n')
#
## Magnitude check
#temp <- abs(object$linear[which(object$linear!=0)])
#cat('Largest entry: ')
#cat(max(abs(temp)))
#cat('\n')
#cat('Smallest entry: ')
#cat(min(abs(temp)))
#cat('\n')
