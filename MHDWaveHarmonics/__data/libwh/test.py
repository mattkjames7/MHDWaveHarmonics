import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

libwh = ctypes.CDLL('./libwh.so')

		
_CppArgSortFlt = libwh.ArgSortFlt
_CppArgSortFlt.argtypes = [ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]	



def TestArgSort(n):
	n = np.int32(n)
	x = np.random.rand(n).astype('float32')
	inds = np.arange(n,dtype='int32')
	
	x_py = np.copy(x)
	i_py = np.argsort(x)
	
	_CppArgSortFlt(n,x,inds)
	
	print('Original: ',x_py)
	print('x now   : ',x)
	
	print('C sorted:')
	print('x: ',x[inds])
	print('i: ',inds)
	
	print('Python Sorted:')
	print('x: ',x_py[i_py])
	print('i: ',i_py)
