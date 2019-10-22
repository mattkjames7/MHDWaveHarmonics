import numpy as np
import matplotlib.pyplot as plt
from .GetFieldLine import GetFieldLine

def PlotFieldLineDensity(pos,Params,fig=None,maps=[1,1,0,0],Overplot=False,**kwargs):
	'''
	Simple routine to plot the modelled plasma mass density along a field line.
	
	Args
	=====
		pos: 3 element ndarray containing the cartesian position to trace from in Rp.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		fig: Pyplot instance, useful if plotting as a subplot on a pre-existing figure, by default will create new figure.
		maps: 4 element list containing subplot information [xmaps,ymaps,xmap,ymap] where xmaps and ymaps represent the number of subplots in the x and y directions and xmap and ymap are the specific x and y indices to plot to.
		Overplot: When set to True, this will inhibit the creation of a new subplot, fig must be supplied with a pre-existing figure, where the routine will plot over.

	Keyword Args
	============
		The following keyword argument form **kwargs - they are completely 
		optional and depend on which model is being used.
		
		Model		kwargs
		T89			iopt,Kp,Vx,Vy,Vz,tilt,CoordIn,CoordOut,Alt,MaxLen,DSMax,FlattenSingleTraces,Verbose
		T96			parmod,Pdyn,SymH,By,Bz,Vx,Vy,Vz,tilt,CoordIn,CoordOut,Alt,MaxLen,DSMax,FlattenSingleTraces,Verbose
		T01			parmod,Pdyn,SymH,By,Bz,Vx,Vy,Vz,tilt,CoordIn,CoordOut,Alt,MaxLen,DSMax,FlattenSingleTraces,Verbose
		TS05		parmod,Pdyn,SymH,By,Bz,Vx,Vy,Vz,tilt,CoordIn,CoordOut,Alt,MaxLen,DSMax,FlattenSingleTraces,Verbose
		KT17
		
		
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

		
	Returns
	========
		pyplot instance
	
	'''
	tmp = GetFieldLine(pos,**kwargs)
	Polarization = kwargs.get('Polarization','none')
	if Polarization == 'none':
		T,s = tmp
	else:
		T,s,h = tmp
		
	
	Bm = np.sqrt(T.Bx**2.0 + T.By**2.0 + T.Bz**2.0).astype('float32')
	R = np.sqrt(T.x**2.0 + T.y**2.0 + T.z**2.0)
	maxR = np.float32(R.max())
	
	if np.size(Params) == 2:
		p = Params[0]*(maxR/R)**Params[1]
		label='$\\rho_{eq} = $'+'{:5.1f}'.format(Params[0])+', $m = $'+'{:3.1f}'.format(Params[1])
	else:
		Rnorm = R/maxR
		ne = Params[2]*np.exp(-0.5*((Rnorm-1.0)/0.1)**2) + Params[0]*Rnorm**(-Params[1])
		mav = Params[4]*Rnorm**(-Params[3])
		p = ne*mav
		label='$\\rho_{eq} = $'+'{:5.1f}'.format(Params[0]*Params[4])
	
	if fig is None:
		fig = plt
		fig.figure()
	
	if not Overplot:
		ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	else:
		ax = fig.gca()
			
	ax.plot(s,p,label=label)
	
	return fig 
	
