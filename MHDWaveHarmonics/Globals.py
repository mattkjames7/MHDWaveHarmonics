import os
import platform
import ctypes
from numpy.ctypeslib import ndpointer

kt17_loaded=False
geopack_loaded = False

DefArgs = { 'KT14':	{	'Params':	[1.42,7.37,2.16],
						'maxlen':	1000,
						'initstep':	0.025,
						'maxstep':	0.025,
						'LimType': 57,
						'FlattenSingleTraces':	True,
						'SmoothSurfaceTransition': False,
						'AssumeOutsidePlanet': True	},
			'KT17':	{	'Params':	[0.387,50.0],
						'maxlen':	1000,
						'initstep':	0.025,
						'maxstep':	0.025,
						'LimType': 57,
						'FlattenSingleTraces':	True,
						'SmoothSurfaceTransition': False,
						'AssumeOutsidePlanet': True	},
			'T89':	{	'CoordIn':	'SM',
						'CoordOut':	'SM',
						'alt':	100.0,
						'MaxLen':	1000,
						'DSMax':	0.1,
						'FlattenSingleTraces': True,
						'Verbose':	False},
			'T96':	{	'CoordIn':	'SM',
						'CoordOut':	'SM',
						'alt':	100.0,
						'MaxLen':	1000,
						'DSMax':	0.1,
						'FlattenSingleTraces': True,
						'Verbose':	False},
			'T01':	{	'CoordIn':	'SM',
						'CoordOut':	'SM',
						'alt':	100.0,
						'MaxLen':	1000,
						'DSMax':	0.1,
						'FlattenSingleTraces': True,
						'Verbose':	False},
			'TS05':	{	'CoordIn':	'SM',
						'CoordOut':	'SM',
						'alt':	100.0,
						'MaxLen':	1000,
						'DSMax':	0.1,
						'FlattenSingleTraces': True,
						'Verbose':	False},
			'Dipole':{	'ds':	0.1,
						'Beq':	-31200.0,
						'nmax':	10000,
						'bounds':	[1.0,100.0] }}
						
						

DataPath = os.path.dirname(__file__)+'/__data/'

libwh = ctypes.CDLL(DataPath+'libwh/libwh.so')

_CppCalcFieldLineVa = libwh.CalcFieldLineVa
_CppCalcFieldLineVa.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]

_CppCalcFieldLineVaMid = libwh.CalcFieldLineVaMid
_CppCalcFieldLineVaMid.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]

_CppCalcdlndx = libwh.Calcdlndx
_CppCalcdlndx.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]


		
_CppSolveWave = libwh.SolveWaveWrapper
_CppSolveWave.argtypes = [ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	

_CppSolveWaveComplex = libwh.SolveWaveComplexWrapper
_CppSolveWaveComplex.argtypes = [ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),							
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	
							

_CppSolveWaveVa = libwh.SolveWaveWrapperVa
_CppSolveWaveVa.argtypes = [ctypes.c_float,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	

_CppFindHarmonics = libwh.FindHarmonicsWrapper
_CppFindHarmonics.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ctypes.c_float,
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	

_CppFindHarmonicsComplex = libwh.FindHarmonicsComplexWrapper
_CppFindHarmonicsComplex.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_bool, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	



_CppGridMisfit = libwh.GridMisfitWrapper
_CppGridMisfit.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_int,
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_float,
							ctypes.c_bool,
							ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]	

_CppInitMisfitObject = libwh.InitMisfitObject
_CppInitMisfitObject.restype = ctypes.c_int
_CppInitMisfitObject.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_int,
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
							ctypes.c_float]

_CppGetMisfit = libwh.GetMisfit
_CppGetMisfit.restype = ctypes.c_float
_CppGetMisfit.argtypes = [ctypes.c_int,
							ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
							ctypes.c_int,
							ctypes.c_bool]
							
_CppDestroyMisfitObject = libwh.DestroyMisfitObject
_CppDestroyMisfitObject.argtypes = [ctypes.c_int]
						
