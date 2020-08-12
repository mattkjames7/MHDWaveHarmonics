import numpy as np
import matplotlib.pyplot as plt
from .GetFieldLine import GetFieldLine
from .FindHarmonics import FindHarmonics
from .SolveWave import SolveWave
#from Plotting.RemoveAxisLabels import RemoveAxisLabels
from scipy.interpolate import InterpolatedUnivariateSpline,interp1d
#from Plotting.AddPlotLabel import AddPlotLabel

def PlotHarmonics(pos,Params,nh=3,df=1.0,Rp=1.0,Colours=None,Method='Complex',RhoBG=None,ReturnData=False,**kwargs):
	'''
	Plots structure of toroidal and poloidal harmonics of a given field line and plasma model.
	
	Inputs
	======
		pos: This is the initial starting positional vector p = np.array([x,y,z]).
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).

		df: Frequency step size in mHz (default=1.0), smaller values should be used when expecting very low frequencies (not used when Complex method is being use).
		nh: Number of harmonics to search for.
		Rp: Planetary radius (can also be set to 'Mercury' of 'Earth')
		Colours: array of colours used for plotting the harmonics.
		Method: Set to 'Complex' or 'Simple' to use either complex or simple shooting method.
		

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
	 
	'''
	Tp,Sp,hp = GetFieldLine(pos,Polarization='poloidal',**kwargs)
	Tt,St,ht = GetFieldLine(pos,Polarization='toroidal',**kwargs)
	Core = kwargs.get('Core',True)

	#find crossing over z = 0
	n = np.where(Tt.z >= 0.0)[0]
	s = np.where(Tt.z < 0.0)[0]
	use = np.append(s[:-2],n[:2])
	if use.size == 4:
		f = InterpolatedUnivariateSpline(Tt.z[use],St[use])
	else:
		f = interp1d(Tt.z[use],St[use])
	Smid = f(0.0)
	
	
	

	B = np.sqrt(Tp.Bx**2 + Tp.By**2 + Tp.Bz**2)
	B = B[np.isfinite(B)]
	R = np.sqrt(Tp.x**2.0 + Tp.y**2.0 + Tp.z**2.0)[:Tp.nstep]

	mu0 = 4*np.pi*1e-7

	if hasattr(Tp,'Rmso'):
		Rgeo = Tp.Rmso
		Rmag = Tp.Rmsm
	else:
		Rgeo = np.sqrt(Tp.x**2.0 + Tp.y**2.0 + Tp.z**2.0)
		Rmag = np.sqrt(Tp.x**2.0 + Tp.y**2.0 + Tp.z**2.0)
	
	if hasattr(Tp,'InPlanet'):
		InPlanet=Tp.InPlanet
	else:
		InPlanet = (Rgeo < 1.0)

	if Rp == 'Mercury':
		Rp = 2440.0
		xlabel = '$x$ ($R_M$)'
	elif Rp == 'Earth':
		Rp = 6380.0
		xlabel = '$x$ ($R_E$)'
	elif Rp == 1.0:
		xlabel = '$x$ (km)'
	else:
		xlabel = '$x$ ($R_P$)'	
	
	if np.size(Params) == 2:
		p = Params[0]*(R.max()/R)**Params[1]
		title = '$\\alpha=$'+'{:3.1f}'.format(Params[1])+'$, \\rho_{eq}=$'+'{:5.1f}'.format(Params[0])+' amu cm$^{-3}$'
	else:
		Rnorm = R/R.max()
		ne = Params[2]*np.exp(-0.5*((Rnorm-1.0)/0.1)**2) + Params[0]*Rnorm**(-Params[1])

		mav = Params[4]*Rnorm**(-Params[3])
		p = ne*mav
		title = '$\\rho_{eq} = $'+'{:5.1f}'.format(Params[0]*Params[4])+' amu cm$^{-3}$'

	
	p = p*1.67377e-27*1e6
	Va = (B*1e-9/np.sqrt(mu0*p))/1000.0
	h2Bp = hp**2 * B
	h2Bt = ht**2 * B

	fp,_,_ = FindHarmonics(Tp,Sp,Params,hp,RhoBG,np.arange(nh)+1,None,df,Method)
	ft,_,_ = FindHarmonics(Tt,St,Params,ht,RhoBG,np.arange(nh)+1,None,df,Method)

	fig = plt
	fig.figure(figsize=(6,8))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	for i in range(0,nh):
		y,_,_= SolveWave(fp[i],Sp,B,Rmag,None,hp,RhoBG,Params,InPlanet)

		y/=hp
		if i == 0:
			ax0.axis([0.0,np.max(Sp/Rp),-np.max(np.abs(y)),np.max(np.abs(y))])
			ax0.plot([0.0,np.max(Sp/Rp)],[0.0,0.0],color=[0.0,0.0,0.0],linestyle='--')
		if Colours is None:
			ax0.plot(Sp/Rp,y,label='f={:5.2f} mHz'.format(fp[i]))
		else:	
			ax0.plot(Sp/Rp,y,label='f={:5.2f} mHz'.format(fp[i]),color=Colours[i])
	R = ax0.axis()
	ax0.vlines(Smid,R[2],R[3],color=[0.0,0.0,0.0],linestyle=':')
	ax0.text(0.2*(R[1]-R[0])+R[0],0.75*(R[3]-R[2])+R[2],'South',ha='center',va='center',color=[0.7,0.7,0.7],zorder=-2.0,fontsize='x-large')
	ax0.text(0.8*(R[1]-R[0])+R[0],0.75*(R[3]-R[2])+R[2],'North',ha='center',va='center',color=[0.7,0.7,0.7],zorder=-2.0,fontsize='x-large')
	ax0.text(Smid,0.5*(R[3]-R[2])+R[2],'Magnetic Equator',va='center',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=-90.0)
	fig.legend(loc='lower right',fontsize=10)
	fig.title(r'Poloidal Resonance '+title)
	fig.xlabel(xlabel)
	fig.ylabel(r'$\xi$ (Arb. Units)')
	
	
	lbl = ax0.get_yticklabels()
	ax0.set_yticklabels(['']*np.size(lbl))	
	R = ax0.axis()
	ax0.text(0.05*(R[1]-R[0])+R[0],0.9*(R[3]-R[2])+R[2],'(a)',ha='center',va='center')
	
	ax1 = fig.subplot2grid((2,1),(1,0))
	for i in range(0,nh):

		y,_,_= SolveWave(ft[i],St,B,Rmag,None,ht,RhoBG,Params,InPlanet)
		y/=ht
		if i == 0:
			ax1.axis([0.0,np.max(St/Rp),-np.max(np.abs(y)),np.max(np.abs(y))])
			ax1.plot([0.0,np.max(St/Rp)],[0.0,0.0],color=[0.0,0.0,0.0],linestyle='--')
		if Colours is None:
			ax1.plot(St/Rp,y,label='f={:5.2f} mHz'.format(ft[i]))
		else:	
			ax1.plot(St/Rp,y,label='f={:5.2f} mHz'.format(ft[i]),color=Colours[i])	
	R = ax1.axis()
	ax1.vlines(Smid,R[2],R[3],color=[0.0,0.0,0.0],linestyle=':')
	ax1.text(0.2*(R[1]-R[0])+R[0],0.75*(R[3]-R[2])+R[2],'South',ha='center',va='center',color=[0.7,0.7,0.7],zorder=-2.0,fontsize='x-large')
	ax1.text(0.8*(R[1]-R[0])+R[0],0.75*(R[3]-R[2])+R[2],'North',ha='center',va='center',color=[0.7,0.7,0.7],zorder=-2.0,fontsize='x-large')
	ax1.text(Smid,0.5*(R[3]-R[2])+R[2],'Magnetic Equator',va='center',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=-90.0)
	fig.legend(loc='lower right',fontsize=10)
	fig.title(r'Toroidal Resonance '+title)
	fig.xlabel(xlabel)
	fig.ylabel(r'$\xi$ (Arb. Units)')	

	
	lbl = ax1.get_yticklabels()
	ax1.set_yticklabels(['']*np.size(lbl))	
	R = ax1.axis()
	ax1.text(0.05*(R[1]-R[0])+R[0],0.9*(R[3]-R[2])+R[2],'(b)',ha='center',va='center')
	
	if Core:

		crust0 = np.where((Rgeo[:-1] < 1.0) & (Rgeo[1:] >= 1.0))[0][0]
		crust1 = np.where((Rgeo[:-1] >= 1.0) & (Rgeo[1:] < 1.0))[0][0]
		S0 = St[crust0]
		S1 = St[crust1]
		R = ax0.axis()

		ax0.vlines(S0,R[2],R[3],color=[0.0,0.0,0.0],linestyle='--')
		ax0.vlines(S1,R[2],R[3],color=[0.0,0.0,0.0],linestyle='--')
		ax0.text(S0,0.5*(R[3]-R[2])+R[2],'Mercury Surface',va='center',ha='right',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=-90.0)
		ax0.text(S1,0.5*(R[3]-R[2])+R[2],'Mercury Surface',va='center',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=90.0)

		R = ax1.axis()
		
		ax1.vlines(S0,R[2],R[3],color=[0.0,0.0,0.0],linestyle='--')
		ax1.vlines(S1,R[2],R[3],color=[0.0,0.0,0.0],linestyle='--')
		ax1.text(S0,0.5*(R[3]-R[2])+R[2],'Mercury Surface',va='center',ha='right',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=-90.0)
		ax1.text(S1,0.5*(R[3]-R[2])+R[2],'Mercury Surface',va='center',color=[0.0,0.0,0.0],zorder=2.0,fontsize='medium',rotation=90.0)
		
	fig.tight_layout()
	
	if ReturnData:
		return fig,Tp,Sp,hp,Tt,St,ht 	
	else:
		return fig


	
