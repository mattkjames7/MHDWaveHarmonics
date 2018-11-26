import numpy as np

def GetSandhuParams(mlt,L):
	'''
	Get default parameters used in Sandhu et al. 2016a,b
	
	Args:
		mlt: Magnetic local time in hours
		L: Farthest point on the field line in Rp
	
	Returns:
		5-element array containing [ne0,alpha,a,beta,mav0]
		ne0: equatorial plasma number density
		alpha: power law for number density
		a: size of Gaussian at the equator
		beta: power law for average ion mass
		mav0: equatorial average ion mass in amu
	
	'''
	ne0 = 35.0 - 3.35*L + (9.38-0.756*L)*np.cos((mlt*15.0 + 76.0)*np.pi/180.0)
	alpha = -0.173 + 0.113*L + 0.412*np.cos((mlt*15.0 + 81.9 + 16.0*L)*np.pi/180.0)
	a = -1.24 + 0.944*L + 2.29*np.cos((mlt*15.0 + 40.0)*np.pi/180.0)
	beta = -2.13 + 0.223*L + (2.26 - 0.218*L)*np.cos((mlt*15.0 + 219.0)*np.pi/180.0)
	mav0 = 16.4 - 1.32*L + (7.21 - 0.665*L)*np.cos((mlt*15.0 + 32.0)*np.pi/180.0)
	
	return np.array([ne0,alpha,a,beta,mav0])
