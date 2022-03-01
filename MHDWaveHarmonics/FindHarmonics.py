import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer
from . import Globals

	
def FindHarmonics(T,s,Params,halpha=None,RhoBG=None,Harmonics=[1,2,3],x0=None,df=1.0,Method='Complex'):
	'''
	Finds harmonic frequencies of waves capable of standing on a given field line.
	
	Args:
		T: TraceField object.
		s: Array storing distance along field line in km.
		halpha: h_alpha array for field line in order to use Singer et al. wave equation, or None for simple wave equation.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		Harmonics: A list or array of the harmonic numbers to be found.
		x0: starting frequency.
		df: Frequency step size in mHz (default=1.0), smaller values should be used when expecting very low frequencies
			this parameter is not used if Method='Complex', but is used for Method='Simple'.
		Method: Set to 'Complex' or 'Simple' to use either complex or simple shooting method.
		
	Returns:
		ndarray containing a list of nh harmonic frequncies in mHz.
		boolean array to say if fit was succesful (only really applicable for Complex method)
		number of iterations used to calculate harmonics.
	
	'''
	
	
	_B = np.sqrt(T.Bx**2.0 + T.By**2.0 + T.Bz**2.0).astype('float32')

	_R = np.sqrt(T.x**2.0 + T.y**2.0 + T.z**2.0).astype('float32')
	_maxR = np.float32(np.nanmax(_R))
	_s = s.astype('float32')
	if not halpha is None:
		_halpha = halpha.astype('float32')
	else:
		_halpha = np.ones(np.size(s),dtype='float32')
	if not RhoBG is None:
		_RhoBG = RhoBG.astype('float32')
	else:
		_RhoBG = np.zeros(np.size(s),dtype='float32')
	_n = T.nstep

	_Params = np.float32(Params)
	_nP = np.int32(np.size(Params))
	
	_HarmInds = np.array(Harmonics,dtype='int32')
	_nh = np.int32(_HarmInds.size)
	_freqs = np.zeros(_nh,dtype='float32')
	
	if x0 is None:
		if Method is 'Complex':
			x0 = np.float32(1.0)
		else:
			x0 = np.float32(0.0)
	
	if Method is 'Complex':
		if np.size(x0) == 1:
			_x0 = np.zeros(_nh,dtype='float32')+x0
		else:
			_x0 = np.array(x0,dtype='float32')
		_nIter = np.zeros(_nh,dtype='int32')
	else:
		_x0 = np.float32(x0)
		_nIter = np.zeros(1,dtype='int32')
	
	_df = np.float32(df)

	
	_Success = np.zeros(_nh,dtype='bool8').astype('bool8')
	if hasattr(T,'InPlanet'):
		_InPlanet = np.float32(T.InPlanet)
	elif hasattr(T,'Rmso'):
		_InPlanet = np.float32(T.Rmso < 1.0)
	else:
		_InPlanet = np.float32(T.R < 1.0)
	
	if Method is 'Complex':
		Globals._CppFindHarmonicsComplex(_B,_R,_s,_halpha,_InPlanet,_RhoBG,_n,_Params,_nP,_maxR,_HarmInds,_nh,_x0,_Success,_nIter,_freqs)
	else:
		Globals._CppFindHarmonics(_B,_R,_s,_halpha,_InPlanet,_RhoBG,_n,_Params,_nP,_maxR,_df,_HarmInds,_nh,_x0,_nIter,_freqs)
					
	
	return _freqs,_Success,_nIter

def FindHarmonicsPMD(B,pmd,s,halpha=None,RhoBG=None,Harmonics=[1,2,3],x0=None,df=1.0,Method='Complex'):
	'''
	Finds harmonic frequencies of waves capable of standing on a given field line.
	
	Args:
		T: TraceField object.
		s: Array storing distance along field line in km.
		halpha: h_alpha array for field line in order to use Singer et al. wave equation, or None for simple wave equation.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		Harmonics: A list or array of the harmonic numbers to be found.
		x0: starting frequency.
		df: Frequency step size in mHz (default=1.0), smaller values should be used when expecting very low frequencies
			this parameter is not used if Method='Complex', but is used for Method='Simple'.
		Method: Set to 'Complex' or 'Simple' to use either complex or simple shooting method.
		
	Returns:
		ndarray containing a list of nh harmonic frequncies in mHz.
		boolean array to say if fit was succesful (only really applicable for Complex method)
		number of iterations used to calculate harmonics.
	
	'''
	
	
	_B = np.array(B).astype('float32')
	_pmd = np.array(pmd).astype('float32')


	_s = s.astype('float32')
	if not halpha is None:
		_halpha = halpha.astype('float32')
	else:
		_halpha = np.ones(np.size(s),dtype='float32')
	if not RhoBG is None:
		_RhoBG = RhoBG.astype('float32')
	else:
		_RhoBG = np.zeros(np.size(s),dtype='float32')
	_n = _B.size
	
	_HarmInds = np.array(Harmonics,dtype='int32')
	_nh = np.int32(_HarmInds.size)
	_freqs = np.zeros(_nh,dtype='float32')
	
	if x0 is None:
		if Method is 'Complex':
			x0 = np.float32(1.0)
		else:
			x0 = np.float32(0.0)
	
	if Method is 'Complex':
		if np.size(x0) == 1:
			_x0 = np.zeros(_nh,dtype='float32')+x0
		else:
			_x0 = np.array(x0,dtype='float32')
		_nIter = np.zeros(_nh,dtype='int32')
	else:
		_x0 = np.float32(x0)
		_nIter = np.zeros(1,dtype='int32')
	
	_df = np.float32(df)

	
	_Success = np.zeros(_nh,dtype='bool8').astype('bool8')
	_InPlanet = np.zeros(_n,dtype='float32')

	
	if Method is 'Complex':
		Globals._CppFindHarmonicsPMDComplex(_B,_s,_halpha,_pmd,_InPlanet,_RhoBG,_n,_HarmInds,_nh,_x0,_Success,_nIter,_freqs)
	else:
		Globals._CppFindHarmonicsPMD(_B,_s,_halpha,_pmd,_InPlanet,_RhoBG,_n,_df,_HarmInds,_nh,_x0,_nIter,_freqs)
					
	
	return _freqs,_Success,_nIter

