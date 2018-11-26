import numpy as np
from . import Globals
import ctypes as ct

def CalcFieldLineVa(T,s,Params,halpha=None):
	'''
	Calculated the Alfven speed along a field line.
	
	Input:
		T: TraceField object.
		s: Distance along the field line in km.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		halpha: halpha value based on the separation of two traces - if not supplied then is filled with ones (simple wave equation).
	
	Returns:
		Array containing Alfven speed.
	'''
	Bm = np.sqrt(T.Bx**2.0 + T.By**2.0 + T.Bz**2.0).astype('float32')

	R = np.sqrt(T.x**2.0 + T.y**2.0 + T.z**2.0).astype('float32')
	maxR = np.float32(np.nanmax(R))
	S = s.astype('float32')
	if not halpha is None:
		Ha = halpha.astype('float32')
	else:
		Ha = np.ones(np.size(s),dtype='float32')
	n = T.nstep

	PAR = np.float32(Params)
	nP = np.int32(np.size(Params))
	
	if hasattr(T,'InPlanet'):
		InPlanet = np.float32(T.InPlanet)
	elif hasattr(T,'Rmso'):
		InPlanet = np.float32(T.Rmso < 1.0)
	else:
		InPlanet = np.float32(T.R < 1.0)
		
	Va = np.zeros(n,dtype='float32')
	Globals._CppCalcFieldLineVa(Bm,R,S,Ha,InPlanet,n,PAR,nP,maxR,Va)
	
	return Va
	
	
def CalcFieldLineVaMid(T,s,Params,halpha=None):
	'''
	Calculated the Alfven speed along a field line.
	
	Input:
		T: TraceField object.
		s: Distance along the field line in km.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		halpha: halpha value based on the separation of two traces - if not supplied then is filled with ones (simple wave equation).
	
	Returns:
		Array containing Alfven speed at the midpoints of s.
	'''	
	Bm = np.sqrt(T.Bx**2.0 + T.By**2.0 + T.Bz**2.0).astype('float32')

	R = np.sqrt(T.x**2.0 + T.y**2.0 + T.z**2.0).astype('float32')
	maxR = np.float32(np.nanmax(R))
	S = s.astype('float32')
	if not halpha is None:
		Ha = halpha.astype('float32')
	else:
		Ha = np.ones(np.size(s),dtype='float32')
	n = T.nstep

	PAR = np.float32(Params)
	nP = np.int32(np.size(Params))
	
	if hasattr(T,'InPlanet'):
		InPlanet = np.float32(T.InPlanet)
	elif hasattr(T,'Rmso'):
		InPlanet = np.float32(T.Rmso < 1.0)
	else:
		InPlanet = np.float32(T.R < 1.0)
	Va = np.zeros(n-1,dtype='float32')
	Globals._CppCalcFieldLineVaMid(Bm,R,S,Ha,InPlanet,n,PAR,nP,maxR,Va)
	
	return Va

