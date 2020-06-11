import errors

# Example 1: beta distribution
ex1 = errors.ErrorDistribution(errors.betapdf,[1.,10.],'beta')
errors.plot_expdfs(ex1) # plot some example PDFs
errors.plot_errors(ex1,'eps') # vary epsilon
errors.plot_errors(ex1,'sig') # vary sigma

# Example 2: truncated normal distribution
ex2 = errors.ErrorDistribution(errors.normpdf,[.001,.75],'truncated normal')
errors.plot_expdfs(ex2) # plot some example PDFs
errors.plot_errors(ex2,'eps') # vary errorspsilon
errors.plot_errors(ex2,'sig') # vary sigma


