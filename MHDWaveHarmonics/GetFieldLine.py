import numpy as np
from FieldTracing import RK4 as rk4
import inspect
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import minimize
from .DipoleField import TraceDipoleField

from . import Globals
try:
	import KT17 as kt17
except:
	try:
		import Models.KT17 as kt17
	except:
		print('KT17 module not found')
	else:
		Globals.kt17_loaded = True
else:
	Globals.kt17_loaded = True

try:
	import PyGeopack as gp
except:
	try:
		import Models.PyGeopack as gp
	except:	
		print('Geopack module not found')
	else:
		Globals.geopack_loaded = True
else:
	Globals.geopack_loaded = True


def _SortModelDirection(T):
	if T.z[0] > T.z[T.nstep-1]:
		I = np.arange(T.nstep)
		O = I[::-1]
		T.x[I] = T.x[O]
		T.y[I] = T.y[O]
		T.z[I] = T.z[O]
		
		T.Bx[I] = T.Bx[O]
		T.By[I] = T.By[O]
		T.Bz[I] = T.Bz[O]
		
		if hasattr(T,'Rmso'):
			T.Rmso[I] = T.Rmso[O]
			T.Rmsm[I] = T.Rmsm[O]
		

def GetModelFunction(Model,ModelArgs):
	if Model in ['T89','T96','T01','TS05','T89c','T96c','T01c','TS05c']:
		if not Globals.geopack_loaded:
			print('Geopack module not installed!!!!')
			return None
		def ModelFunc(p):
			T0 = gp.TraceField(p[0],p[1],p[2],ModelArgs[0],ModelArgs[1],Model,ModelArgs[2],ModelArgs[3],ModelArgs[4], ModelArgs[5],ModelArgs[6],
				ModelArgs[7],ModelArgs[8],ModelArgs[9],ModelArgs[10],ModelArgs[11],ModelArgs[12],FlattenSingleTraces=True)
			return T0
	elif Model == 'KT17':
		if not Globals.kt17_loaded:
			print('KT17 module not installed!!!!')
			return None
		def ModelFunc(p):
			T0 = kt17.TraceField(p[0],p[1],p[2],ModelArgs[0],ModelArgs[1],ModelArgs[2],ModelArgs[3],ModelArgs[4],ModelArgs[5],ModelArgs[6],ModelArgs[7])
			return T0
	elif Model == 'KT14':
		if not Globals.kt17_loaded:
			print('KT17 module not installed!!!!')
			return None
		def ModelFunc(p):
			T0 = kt17.TraceField(p[0],p[1],p[2],ModelArgs[0],ModelArgs[1],ModelArgs[2],ModelArgs[3],ModelArgs[4],ModelArgs[5],ModelArgs[6],ModelArgs[7])
			return T0	
	elif Model == 'Dipole':
		def ModelFunc(p):
			T0 = TraceDipoleField(p,ModelArgs[0],ModelArgs[1],ModelArgs[2],ModelArgs[3])
			return T0
	else:
		print('Model not found')
		return None
	return ModelFunc

def GetFieldLine(pos,Model='KT14',MaxStepSize=None,ModelArgs=None,Delta=None,Polarization='none',Core=False):
	'''
	Provides a field line trace with distance along field line and optionally h_alpha.
	
	Args:
		pos: ndarray with cartesian positional vector in Rp.
		Model: Model field name -- 'KT14'|'KT17'|'T89'|'T96'|'T01'|'TS05'|'Dipole'
				or 'T89c'|'T96c'|'T01c'|'TS05c' for use with customised model parameters
		MaxStepSize: Trace step size maximum in Rp (default = None).
		ModelArgs: 	Tuple containing arguments for the magnetic field models, 
					when set to None, a set of default parameters are used.
					'T89'|'T96'|'T01'|'TS05'|'T89c'|'T96c'|'T01c'|'TS05c': 
						ModelArgs = (Date,ut,CoordIn,CoordOut,Alt,MaxLen,DSMax,iopt,parmod,tilt,Vx,Vy,Vz)
						****NOTE****
						When using models 'T89'|'T96'|'T01'|'TS05' - the iopt,parmod,tilt,Vx,Vy,Vz parameters need not
						be specified as they will be calculated automatically within the Geopack model from Omni data
						using the Date and ut parameters
						
						To use those parameters, add 'c' tot he model name string: 'T89c'|'T96c'|'T01c'|'TS05c'
						Then all parameters will be needed
						
						************
						Date:  Date in format yyyymmdd.
						UT: Time in format hh.hh (hh + mm/60.0 + ss/3600.0).
						CoordIn: Coordinate system of input position, by default 'SM' (Solar-Magnetic), can be set to 'GSM' (Geocentric Solar Magnetospheric).
						CoordOut: Coordinate system of output positions and field vectors, by default 'SM', can be set to 'GSM'.
						Alt: Altitude to stop tracing at - default = 100km
						MaxLen: maximum number of trace steps
						DSMax: Maximum step size
						iopt: integer for controlling T89c model
						parmod: 10-element floating point array to control T96c,T01c and TS05c models,
								for T96:
									parmod[0] = Pdyn (nPa)
									parmod[1] = Dst (nT)
									parmod[2] = IMF By (nT)
									parmod[3] = IMF Bz (nT)
								for T01:
									parmod[0] = Pdyn (nPa)
									parmod[1] = Dst (nT)
									parmod[2] = IMF By (nT)
									parmod[3] = IMF Bz (nT)
									parmod[4] = G1 parameter (See Tsyganenko [2001])
									parmod[5] = G2 parameter (See Tsyganenko [2001])
								for TS05:
									parmod[0] = Pdyn (nPa)
									parmod[1] = Dst (nT)
									parmod[2] = IMF By (nT)
									parmod[3] = IMF Bz (nT)
									parmod[4] = W1 parameter (See Tsyganenko and Sitnov [2005])
									parmod[5] = W2 parameter (See Tsyganenko and Sitnov [2005])	
									parmod[6] = W3 parameter (See Tsyganenko and Sitnov [2005])
									parmod[7] = W4 parameter (See Tsyganenko and Sitnov [2005])
									parmod[8] = W5 parameter (See Tsyganenko and Sitnov [2005])
									parmod[9] = W6 parameter (See Tsyganenko and Sitnov [2005])		
						tilt: Geodipole tilt angle - if set to NaN then will be calculated using Date and ut
						Vx,Vy,Vz: IMF velocity components
					KT17:
						ModelArgs = (Rsun,DistIndex,MaxLen,InitStep,MaxStepSize,LimType)
						Rsun: radial distance from the Sun in AU
						DistIndex: Disturbance index (0.0 - 100.0)
						MaxLen: Maximum number of steps for trace
						InitStep: Starting step size.
						MaxStepSize: Maximum step size
						LimType: Integer value to define where to stop the field trace (default 0):
							0: Terminate trace at planet surface and magnetopause
							1: Confine to box -6 < x < 2, -4 < y < 4, -4 < z < 4
							2: Confine to box and terminate at planet
							3: Terminate trace at planet
							4: Trace to MP, Planet and stop at 10Rm
								
					KT14:
						ModelArgs = (MaxLen,InitStep,MaxStep,LimType,Rsm,t1,t2)
						MaxLen: Maximum number of steps for trace
						InitStep: Starting step size.
						MaxStepSize: Maximum step size
						LimType: Integer value to define where to stop the field trace (default 0):
							0: Terminate trace at planet surface and magnetopause
							1: Confine to box -6 < x < 2, -4 < y < 4, -4 < z < 4
							2: Confine to box and terminate at planet
							3: Terminate trace at planet
							4: Trace to MP, Planet and stop at 10Rm
						Rsm: Subsolar magnetopause radius (default=1.42).
						t1: Tail disk current strength (default=7.37).
						t2: Tail quasi-harris current sheet strength (default=2.16).
					Dipole:
						ModelArgs = [Beq]
						Beq: Megnatic field strength at equator in nT (Default=-31200.0).
		Delta: Separation between two traced field lines
		Polarization: 'none'|'toroidal'|'poloidal'
		Core: This only applies to the KT14/KT17 field, as it will include tracing to the core of the planet rather than the surface.
	
	Returns:
		TracedField object
		Array storing distance along the field line.
		if Polarization is 'toroidal' or 'poloidal', then h_alpha array is also returned
	
	'''
	
	#First of all check if the model that has been requested is imported or not!
	if Model in ['KT14','KT17'] and not Globals.kt17_loaded:
		#If the kt17 model hasn't been installed in the expected location
		#then this will happen
		print('KT17 module is required for this model')
		return (None,None,None)
		
	if Model in ['T89','T96','T01','TS05','T89c','T96c','T01c','TS05c'] and not Globals.geopack_loaded:
		#If the geopack model hasn't been installed in the expected location
		#then this will happen
		print('Geopack module is required for this model')
		return (None,None,None)	
	
	#now to trace the model field
	if Model in ['T89','T96','T01','TS05','T89c','T96c','T01c','TS05c']:
		if MaxStepSize is None:
			MaxStepSize = 0.1
		if (np.size(ModelArgs) != 13):
			#Date,ut,CoordIn,CoordOut,Alt,MaxLen,DSMax,iopt,parmod,tilt,Vx,Vy,Vz
			ModelArgs = (20000101,12.0,'SM','SM',100.0,1000,MaxStepSize,None,None,None,None,None,None)
		if Delta is None:
			Delta = 0.05
		Rp = 6380.0
	elif Model == 'KT17':
		if MaxStepSize is None:
			MaxStepSize = 0.025
		if (np.size(ModelArgs) != 8):
			#Rsun,DistIndex,MaxLen,InitStep,MaxStep,LimType
			#old model -> ModelArgs = (0.387,50.0,1000,MaxStepSize,MaxStepSize,0)
			if Core:
				LimType = 57
			else:
				LimType = 15
			ModelArgs = (1000,MaxStepSize,MaxStepSize,LimType,[0.387,50.0],True,False,True)
		if Delta is None:
			Delta = 0.05	
		Rp = 2440.0
	elif Model == 'KT14':
		if MaxStepSize is None:
			MaxStepSize = 0.025
			#MaxStepSize = 0.001
		if (np.size(ModelArgs) != 8):
			#MaxLen,InitStep,MaxStep,LimType,Rsm,T1,T2
			#old - > ModelArgs = (1000,MaxStepSize,MaxStepSize,0,1.42,7.37,2.16)

			if Core:
				LimType = 57
			else:
				LimType = 15
			ModelArgs = (1000,MaxStepSize,MaxStepSize,LimType,[1.42,7.37,2.16],True,False,True)
		if Delta is None:
			Delta = 0.05	
		Rp = 2440.0
	elif Model == 'Dipole':
		if MaxStepSize is None:
			MaxStepSize = 0.1
		if (np.size(ModelArgs) != 4):
			#ds,Beq,nmax,bounds
			ModelArgs = (MaxStepSize,-31200.0,10000,[1.0,100.0])
		if Delta is None:
			Delta = 0.05
		Rp = 6380.0
	else:
		print('Model not recognised')
		return
	
	ModelFunc = GetModelFunction(Model,ModelArgs)
	
	T0 = ModelFunc(pos)
	_SortModelDirection(T0)

	#return T0
	s0 = np.zeros(T0.nstep)
	for i in range(1,T0.nstep):
		s0[i] = s0[i-1] + np.sqrt((T0.x[i]-T0.x[i-1])**2+(T0.y[i]-T0.y[i-1])**2+(T0.z[i]-T0.z[i-1])**2)*Rp


	if Polarization == 'none':
		return (T0,s0)
	elif Polarization == 'poloidal':
		eqpos1 = np.array([-(T0.Lshell-Delta)*np.cos(T0.MltE*15.0*np.pi/180.0),-(T0.Lshell-Delta)*np.sin(T0.MltE*15.0*np.pi/180.0),0.0])
	elif Polarization == 'toroidal':
		dmlt = 2.0*np.arcsin(Delta*0.5/T0.Lshell)*180.0/(np.pi*15.0)
		eqpos1 = np.array([-T0.Lshell*np.cos((T0.MltE+dmlt)*15.0*np.pi/180.0),-T0.Lshell*np.sin((T0.MltE+dmlt)*15.0*np.pi/180.0),0.0])
	T1 = ModelFunc(eqpos1)
	_SortModelDirection(T1)
	s1 = np.zeros(T1.nstep)
	for i in range(1,T1.nstep):
		s1[i] = s1[i-1] + np.sqrt((T1.x[i]-T1.x[i-1])**2+(T1.y[i]-T1.y[i-1])**2+(T1.z[i]-T1.z[i-1])**2)*Rp

	#return T0,T1,s0,s1
	d = np.zeros(T0.nstep,dtype='float32')

	fx = InterpolatedUnivariateSpline(s1,T1.x[0:s1.size])
	fy = InterpolatedUnivariateSpline(s1,T1.y[0:s1.size])
	fz = InterpolatedUnivariateSpline(s1,T1.z[0:s1.size])

	def PosFn(s):
		return fx(s),fy(s),fz(s)
		
	def DistFn(s,j):
		return np.sqrt((fx(s)-T0.x[j])**2.0 + (fy(s)-T0.y[j])**2.0 + (fz(s)-T0.z[j])**2.0)
		
	for i in range(0,T0.nstep):
		res = minimize(DistFn,s0[i],args=(i),method='Nelder-Mead')
		d[i] = DistFn(res.x,i)*Rp
		
	h_alpha = d/(Delta*Rp)
	
	return (T0,s0,h_alpha)
	
	
	
