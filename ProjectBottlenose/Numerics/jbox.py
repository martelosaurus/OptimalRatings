import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
#from error_plots import par_var

#==============================================================================#
# EXAMPLE 1: THE BETA DISTRIBUTION                                             #  
#==============================================================================#
def betapdf(x,eps,sig):
	f =
	
	# solve for the parameter 'a' for which the noise has variance sig0**2
	a0 = optimize.newton(lambda a: BetaVar(a)-sig0**2,x0=2.)
	return stats.beta.pdf(x,a,a,loc=-eps,scale=2*eps)

error_plots.par_var(betapdf,'betapdf')

#==============================================================================#
# EXAMPLE 1: THE TRUNCATED NORMAL DISTRIBUTION                                 #  
#==============================================================================#
def normpdf(x,eps,sig):
	pdf_raw = stats.norm.pdf(x,loc=0,scale=a)
	cdf_dn  = stats.norm.cdf(-eps,loc=0.,scale=a)
	cdf_up  = stats.norm.cdf(+eps,loc=0.,scale=a)
	return pdf_raw/(cdf_up-cdf_dn)

error_plots.par_var(normpdf,'normpdf')


# variance of the symmetric beta distribution with parameter 'a' 
def BetaVar(a): 
	return stats.beta.expect(lambda x: x**2,(a,a),-eps0,2*eps0) 

