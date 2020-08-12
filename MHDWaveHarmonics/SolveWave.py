import numpy as np
from . import Globals
import ctypes
from numpy.ctypeslib import ndpointer



def SolveWave(f,x,B,R=None,Va=None,halpha=None,RhoBG=None,Params=None,InPlanet=None,Method='Complex',Unscale=True):
	'''
	Solves wave equation using the Runge-Kutta-Gill method
	
	Args:
		f: Wave frequency in mHz.
		x: Array of positions along a field line in km, must be monotonic with length N.
		B: Magnetic field magnitude along the field line in nT, length N.
		Va: Alfven speed at each point along the field line in km/s, length N (Not used if Params are specified).
		R: Radial distance from centre of planet along the field line in Rp, length N. Not necessary if supplying Va.
		halpha: h_alpha parameter for each position along field line with length N
				is required to solve the Singer wave equation, set to None for basic wave equation.
		Params: Ignored if Va is supplied to the function.
				For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		Method: Set to 'Complex' or 'Simple' to use either complex or simple shooting method.
		InPlanet: floating point array where each element denotes whether that part of the trace is
				inside the planet (1.0) or not (0.0), values between 0 and 1 can be used to smooth the
				transition. 
		Unscale: By deafault - the complex method of solving the wave equation normalized the components,
				setting this to True will return the array with its original amplitude.
				
	Returns:
		Solution to wave equation as ndarray, length N. If halpha is defined then the output is 
		xi/h_alpha, if halpha=None then the output is simply xi. 
	
	'''
	
	
	_f = np.float32(f)
	_x = np.float32(x)
	if not R is None:
		_R = np.float32(R)
		Rmax = np.nanmax(R) 
	if not Va is None:
		_Va = np.float32(Va)
	_B = np.float32(B)
	if not halpha is None:
		_halpha = np.float32(halpha)
	else:
		_halpha = np.ones(np.size(x),dtype='float32')
	if not Params is None:
		_Params = np.float32(Params)
		_nP = np.size(_Params)
	if InPlanet is None:
		_InPlanet = np.zeros(np.size(x),dtype='float32')
	else:
		_InPlanet = np.float32(InPlanet)
	if not RhoBG is None:
		_RhoBG = RhoBG.astype('float32')
	else:
		_RhoBG = np.zeros(np.size(x),dtype='float32')
	_n = np.int32(np.size(_x))
	_yr = np.zeros(_n,dtype='float32')
	_yi = np.zeros(_n,dtype='float32')
	_phase = np.zeros(_n,dtype='float32')
	_mxr = np.array([0.0],dtype='float32')
	_mxi = np.array([0.0],dtype='float32')

	if Params is None:
		#use predefined Alfven speed array
		if Va is None:
			print("Please either supply an array of Va, or a set of Params")
			return None
		else:
			Globals._CppSolveWaveVa(_f,_B,_Va,_x,_halpha,_InPlanet,_n,_yr)
			if Unscale:
				return _yr
			else:
				return _yr/np.nanmax(np.abs(_yr))
	elif Method is 'Simple':
		if R is None:
			print("Please supply an array for R")
			return None
		else:
			Globals._CppSolveWave(_f,_B,_R,_x,_halpha,_InPlanet,_RhoBG,_n,_Params,_nP,Rmax,_yr)
			if Unscale:
				return _yr
			else:
				return _yr/np.nanmax(np.abs(_yr))
	elif Method is 'Complex':
		if R is None:
			print("Please supply an array for R")
			return None
		else:
			Globals._CppSolveWaveComplex(_f,_B,_R,_x,_halpha,_InPlanet,_RhoBG,_n,_Params,_nP,Rmax,_yr,_yi,_phase,_mxr,_mxi)
			if Unscale:
				return _yr*_mxr[0],_yi*_mxi[0],_phase
			else:
				return _yr,_yi,_phase
	else:
		print("Method should be either 'Complex' or 'Simple'")
		return	None
									
