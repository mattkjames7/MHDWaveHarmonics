import numpy as np
from .GetFieldLine import GetFieldLine
from scipy.optimize import minimize
from .FindHarmonics import FindHarmonics

def GetMisfitFunction(T,s,halpha,f,Params,Harm=1,df=1.0,Method='Complex',RhoBG=None):
	'''
	Returns function that calculated the difference between the desired frequency and the modelled frequency given a plasma mass density.
	
	Args:
		T: TraceField object.
		s: Array containing distance along a traced field line in km.
		halpha: An array containing h_alpha values for the traced field line in order to solve Singer et al. wave equation, or None to use simple wave equation.
		f: Frequency of wave in mHz.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
				In both cases, the first element is not important, this will be overwritten.
		Harm: Harmonic number to fit to - default = 1
		df: Frequency step for FindHarmonics routine, for exceptionally low frequencies (< ~5 mHz) use a smaller value than the default (ignored for Complex shooting method).
		Method: Which shooting method to use: 'Complex' or 'Simple'
		
	Returns:
		Function which will return a frequency misfit (absolute) given a plasma mass density.
		
	'''	
	niter = 0
	def CalculateMisfit(x):
		global niter
		Par = np.copy(Params)
		if np.size(Par) == 2:
			Par[0] = x
		else:
			Par[0] = x/Par[-1] #convert peq to neq by dividing by average mass
			#Par[2] = x
		fout,_,_ = FindHarmonics(T,s,Par,halpha,RhoBG,[Harm],None,df,Method)
		misfit = np.abs(f-fout[0])
		niter += 1
		if np.size(Par) == 2:
			print('\rIteration: {:5d}, Cost: {:12.8f}, p_eq: {:7.2f}, Power: {:4.1f}'.format(niter,misfit,x[0],Par[1]),end='')
		else:
			print('\rIteration: {:5d}, Cost: {:12.8f}, n_eq: {:7.2f}, alpha: {:4.1f}, a: {:7.2f}, beta: {:4.1f}, m_av0: {:7.2f}'.format(niter,misfit,Par[0],Par[1],Par[2],Par[3],Par[4]),end='')
		return misfit
	return CalculateMisfit

def FitPlasmaToHarmonic(T,s,halpha,f,Params,Harm=1,df=1.0,Method='Complex',RhoBG=None):
	'''
	Numerically find an equatorial plasma mass density that would allow a wave of a given frequency to exist on a field line with a specific power law.
	
	Args:
		T: TraceField object.
		s: Array containing distance along a traced field line in km.
		halpha: An array containing h_alpha values for the traced field line in order to solve Singer et al. wave equation, or None to use simple wave equation.
		f: Frequency of wave in mHz.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
				In both cases, the first element of the arrays will be taken as the initial value for equatorial plasma mass density.
		Harm: Harmonic number to fit to - default = 1
		df: Frequency step for FindHarmonics routine, for exceptionally low frequencies (< ~5 mHz) use a smaller value than the default (ignored for Complex shooting method).
		Method: Which shooting method to use: 'Complex' or 'Simple'
	
	Returns:
		Plasma mass density in amu/cm^3
		
	'''
	Func = GetMisfitFunction(T,s,halpha,f,Params,Harm,df,Method,RhoBG=RhoBG)
	global niter
	niter = 0

	

	res = minimize(Func,Params[0],method='Nelder-Mead',tol=1e-5,options={'maxiter':1000})
	
	if np.isnan(Func(res.x)):
		p_eq = np.nan
	else:	
		p_eq = res.x[0]
	print()
	return p_eq
