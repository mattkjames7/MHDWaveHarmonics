import numpy as np
from .SolveWave import SolveWave
from scipy.interpolate import InterpolatedUnivariateSpline,interp1d
from .FindHarmonics import FindHarmonics
from scipy.optimize import minimize
import matplotlib.pyplot as plt
#from Plotting.PlotGrid import PlotGrid
import ctypes
from numpy.ctypeslib import ndpointer
from . import Globals


# def GetMisfitFunction2(T,s,halpha,freqs,harms,df=1.0,Params=None,nmax=10000):
	# '''
	# Returns a function which calculates the misfit between modelled harmonics and the supplied frequency peaks given a plasma mass density and power law index.
	
	# Args:
		# T: TraceField object or list of TraceField objects.
		# s: Array of distance along the traced field line, or a list of Arrays.
		# halpha: h_alpha values for each field line.
		# freqs: Array of frequencies to fit to.
		# harms: Array containing harmonic numbers of the frequencies in the freqs array.
		# df: 
		
	# Returns:
		# Function
	# '''
	# def CalculateMisfit(x):
		# misfit = 0.0
		# nT = np.size(T)
		# nf = np.size(freqs)
		# if Params is None:
			# Par = x
		# else:
			# Par = np.copy(Params)
			# Par[0] = x[0]
			# Par[2] = x[1]
		# if nT == 1:
			# fout = FindHarmonics(T,s,Par,halpha,df=df,nh=harms.max(),nmax=nmax)
			# for i in range(0,nf):
				# misfit += ((freqs[i] - fout[harms[i]-1])**2)/(harms[i]**2)
		# else:
			# for i in range(0,nf):
				# fout = FindHarmonics(T[i],s[i],x,halpha[i],df=df,nh=harms[i],nmax=nmax)
				# misfit += ((freqs[i] - fout[harms[i]-1])**2)#/(harms[i]**2)
		# if np.size(Par) == 2:
			# print('\rCost: {:12.8f}, p_eq: {:7.2f}, Power: {:4.1f}'.format(np.sqrt(misfit/nf),x[0],x[1]),end='')
		# else:
			# print('\rCost: {:12.8f}, n_eq: {:7.2f}, a: {:7.2f}'.format(np.sqrt(misfit/nf),x[0],x[1]),end='')
	
		# return np.sqrt(misfit/nf)
	# return CalculateMisfit
	
	
def _GetMisfitFunctionInstance(T,s,halpha,freqs,harms,Params0,ParamFit,df=0.1,Method='Complex'):
	_B,_R,_s,_halpha,_InPlanet,_n,_nTrace,_nsteps,_maxR,_freqs,_harms,_Complex,_df = _ConvertToCppInput(T,s,halpha,freqs,harms,df,Method)
	instance = Globals._CppInitMisfitObject(_B,_R,_s,_halpha,_InPlanet,_n,_nTrace,_nsteps,_maxR,_freqs,_harms,_df)
	_nP = np.int32(np.size(ParamFit))
	_Params = np.zeros(_nP,dtype='float32')
	use = np.where(ParamFit)[0]
	nx = use.size
	def MF(x):
		pos = 0
		for i in range(0,_nP):
			if ParamFit[i]:
				_Params[i] = x[pos]
				pos+=1
			else:
				_Params[i] = Params0[i]
		out = Globals._CppGetMisfit(instance,_Params,_nP,_Complex)
		if _nP == 2:
			print('\rCost: {:12.8f}, p_eq: {:7.2f}, Power: {:4.1f}'.format(out,x[0],x[1]),end='')
		else:
			print('\rCost: {:12.8f}, n_eq: {:7.2f}, alpha: {:7.2f}, a: {:7.2f}, beta: {:7.2f}, mav_eq: {:7.2f}'.format(np.sqrt(misfit/nf),_Params[0],_Params[1],_Params[2],_Params[3],_Params[4]),end='')

		return out
	
	return instance,MF

		

def FitPlasma(T,s,halpha,freqs,harms,Params0,df=1.0,Method='Complex',ParamFit=None):
	'''
	This routine attempts to fit a plasma mass density and power law to a set of frequencies on a single or multiple field lines
	
	
	Args:
		T: TraceField object or list of TraceField objects.
		s: Array of distance along the traced field line, or a list of Arrays.
		halpha: h_alpha values for each field line.
		freqs: Array of frequencies to fit to.
		harms: Array containing harmonic numbers of the frequencies in the freqs array.
		Params0: For power law: 2-element array/list [p_eq,power].
				 For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		
		NOTE: Initial guesses for p_eq0 and power0 should be close to what the actual solution might be, so using GetMisfitGrid might help finding a good starting value.
		Also, it is highly likely that this function won't provide a physical solution due to the insensitivity the the power law index.
	
	Returns:
		Best fitting plasma mass density and power law.
		
	'''
	
	if ParamFit is None:
		ParamFit = np.ones(np.size(Params0),dtype='bool')
	X = np.array(Params0)[ParamFit]
	
	Inst,MF = _GetMisfitFunctionInstance(T,s,halpha,freqs,harms,Params0,ParamFit,df,Method)
	
	res = minimize(MF,X,method='nelder-mead')
	ParamOut = np.copy(Params0)
	ParamOut[ParamFit] = res.x
	print()
	print(res)
	return ParamOut


def _ConvertToCppInput(T,s,halpha,freqs,harms,df=0.1,Method='Complex'):
	if not isinstance(T,list):
		T = [T]
	if not isinstance(s,list):
		s = [s]
	if not isinstance(halpha,list):
		halpha = [halpha]


	nT = len(T)
	nS = 0
	for si in s:
		nS = np.max([nS,np.size(si)])
	
	_nTrace = np.int32(nT)
	_B = np.zeros((nT,nS),dtype='float32')+np.nan
	_R = np.zeros((nT,nS),dtype='float32')+np.nan
	_maxR = np.zeros((nT),dtype='float32')+np.nan
	_s = np.zeros((nT,nS),dtype='float32')+np.nan
	_InPlanet = np.zeros((nT,nS),dtype='float32')+np.nan
	_nsteps = np.zeros(nT,dtype='int32')
	_halpha = np.ones((nT,nS),dtype='float32')

	for i in range(0,nT):
		_nsteps[i] = np.size(s[i])
		Bm = np.sqrt(T[i].Bx**2.0 + T[i].By**2.0 + T[i].Bz**2.0).astype('float32')
		_B[i][:_nsteps[i]] = Bm[:_nsteps[i]]
		Rm = np.sqrt(T[i].x**2.0 + T[i].y**2.0 + T[i].z**2.0).astype('float32')
		_R[i][:_nsteps[i]] = Rm[:_nsteps[i]]
		_maxR[i] = np.float32(np.nanmax(_R[i]))
		_s[i][:_nsteps[i]] = s[i][:_nsteps[i]].astype('float32')
		if not halpha is None:
			_halpha[i][:_nsteps[i]] = halpha[i][:_nsteps[i]]
		if hasattr(T[i],'InPlanet'):
			_InPlanet[i][:_nsteps[i]] = T[i].InPlanet[:_nsteps[i]]
	
	_B = _B.flatten()
	_R = _R.flatten()
	_maxR = _maxR.flatten()
	_s = _s.flatten()
	_halpha = _halpha.flatten()
	
	_n = np.size(_R)

	_df = np.float32(df)
	
	_freqs = np.array(freqs).astype('float32')
	_harms = np.array(harms).astype('int32')
	_Complex = np.bool8(Method is 'Complex')
	
	return _B,_R,_s,_halpha,_InPlanet,_n,_nTrace,_nsteps,_maxR,_freqs,_harms,_Complex,_df

	
def GetMistfitGrid(T,s,halpha,freqs,harms,ParamRange=[[-1.0,2.0],[-6.0,6.0]],dParam=[0.25,1.0],df=1.0,Method='Complex'):
	'''
	Calculates the misfit between modelled and supplied harmonic frequncies for a range 
	of plasma mass densities and power law indices, each on different field traces.
	
	Args:
		T: TraceField object or list of TraceField objects.
		s: Array of distance along the traced field line, or a list of Arrays.
		halpha: h_alpha values for each field line.
		freqs: Array of frequencies to fit to.
		harms: Array containing harmonic numbers of the frequencies in the freqs array.
		logpeqrange: Range for the log_10 of the equatorial plasma mass density.
		dlogpeq: Step size for evaluations within logpeqrange.
		powrange: Range for power law indices to evaluate.
		dpow: Step size for evaluations within powrange.
		df: Frequency step size for FindHarmonics function.
		nmax: Maximum numberof iterations for FindHarmonics function.
		
	Returns:
		logpeqaxis,powaxis,grid,logpeq,power
		logpeqaxis: Bin ranges for log10 of the plasma mass density.
		powaxis: Bin ranges for the power law index.
		grid: 2D array containing misfit values for each logpeq and power combination.
		logpeq: Center of logpeqaxis bins
		power: Center of pow axis bins.
		
		
	'''
	nP = np.size(dParam)
	if (nP != 2) & (nP != 5):
		print("The number of input parameters must equal 2 for the power law model, or 5 for the Sandhu model!")
	
	ParAx = [] #full list of paramter axes
	ParC = [] #central points of the above lists
	ParN = []
	for i in range(0,nP):
		if np.size(ParamRange[i]) == 2:
			n = np.int32(np.round((ParamRange[i][1]-ParamRange[i][0])/dParam[i]))	
			ParN.append(n)
			ParAx.append(np.arange(n+1)*dParam[i] + ParamRange[i][0])
			ParC.append((np.arange(n)+0.5)*dParam[i] + ParamRange[i][0])
		else:
			ParN.append(N)
			ParAx.append(ParamRange[i])
			ParC.append(ParamRange[i])
	
	ParN = np.array(ParN)

	nDim = np.sum(ParN > 1)
	nComb = np.prod(ParN)

	#convert Peq or neq to linear scale
	ParC[0] = 10**ParC[0]
	
	
	M = np.meshgrid(*ParC)
	_Params = np.zeros((nComb,nP),dtype='float32')
	for i in range(0,nP):
		_Params[:,i] = M[i].flatten()
	
#	nT = len(T)
#	nS = 0
#	for si in s:
#		nS = np.max([nS,np.size(si)])
	
#	_nTrace = np.int32(nT)
#	_B = np.zeros((nT,nS),dtype='float32')+np.nan
#	_R = np.zeros((nT,nS),dtype='float32')+np.nan
#	_maxR = np.zeros((nT),dtype='float32')+np.nan
#	_s = np.zeros((nT,nS),dtype='float32')+np.nan
#	_InPlanet = np.zeros((nT,nS),dtype='float32')+np.nan
#	_nsteps = np.zeros(nT,dtype='int32')
#	_halpha = np.ones((nT,nS),dtype='float32')

#	for i in range(0,nT):
#		_nsteps[i] = np.size(s[i])
#		Bm = np.sqrt(T[i].Bx**2.0 + T[i].By**2.0 + T[i].Bz**2.0).astype('float32')
#		_B[i][:_nsteps[i]] = Bm[:_nsteps[i]]
#		Rm = np.sqrt(T[i].x**2.0 + T[i].y**2.0 + T[i].z**2.0).astype('float32')
#		_R[i][:_nsteps[i]] = Rm[:_nsteps[i]]
#		_maxR[i] = np.float32(np.nanmax(_R[i]))
#		_s[i][:_nsteps[i]] = s[i][:_nsteps[i]].astype('float32')
#		if not halpha is None:
#			_halpha[i][:_nsteps[i]] = halpha[i][:_nsteps[i]]
#		if hasattr(T[i],'InPlanet'):
#			_InPlanet[i][:_nsteps[i]] = T[i].InPlanet[:_nsteps[i]]
	
	#_B = _B.flatten()
	#_R = _R.flatten()
	#_maxR = _maxR.flatten()
	#_s = _s.flatten()
	#_halpha = _halpha.flatten()
	
	#_n = np.size(_R)
	_Par = _Params.flatten()
	_nP = np.int32(nP)
	#_df = np.float32(df)
	
	#_freqs = np.array(freqs).astype('float32')
	#_harms = np.array(harms).astype('int32')
	
	_grid = np.zeros(nComb,dtype='float32')
	_nG = np.int32(nComb)

	#_Complex = np.bool8(Method is 'Complex')

	_B,_R,_s,_halpha,_InPlanet,_n,_nTrace,_nsteps,_maxR,_freqs,_harms,_Complex,_df = _ConvertToCppInput(T,s,halpha,freqs,harms,df,Method)


	Globals._CppGridMisfit(_B,_R,_s,_halpha,_InPlanet,_n,_nTrace,_nsteps,_maxR,_freqs,_harms,_Par,_nP,_df,_Complex,_nG,_grid)
	grid = _grid.reshape(tuple(ParN[::-1]))
	#convert Peq or neq to log scale again
	ParC[0] = np.log10(ParC[0])

	return ParAx,ParC,ParN,grid
	
	
def PlotGridMisfit(T,s,halpha,freqs,harms,ParamRange=[[-1.0,2.0],[-6.0,6.0]],dParam=[0.25,1.0],PlotAx=[0,1],AxInds=None,df=1.0,Method='Complex',fig=None,maps=[1,1,0,0],scale=None,nobar=False):
	'''
	Plots misfit grid created by GetMisfitGrid routine.
	
	Args:
		T: TraceField object or list of TraceField objects.
		s: Array of distance along the traced field line, or a list of Arrays.
		halpha: h_alpha values for each field line.
		freqs: Array of frequencies to fit to.
		harms: Array containing harmonic numbers of the frequencies in the freqs array.
		logpeqrange: Range for the log_10 of the equatorial plasma mass density.
		dlogpeq: Step size for evaluations within logpeqrange.
		powrange: Range for power law indices to evaluate.
		dpow: Step size for evaluations within powrange.
		fig: Pyplot instance, useful if plotting as a subplot on a pre-existing figure, by default will create new figure.
		maps: 4 element list containing subplot information [xmaps,ymaps,xmap,ymap] where xmaps and ymaps represent the number of subplots in the x and y directions and xmap and ymap are the specific x and y indices to plot to.
		scale: 2 element list or array containing the limits for the color bar, None will force the routine to calculate this automatically.
		nobar: Inhibits the creation of a color bar.
		df: Frequency step size for FindHarmonics function.
		nmax: Maximum numberof iterations for FindHarmonics function.
			
	
	'''
	ParAx,ParC,ParN,grid =GetMistfitGrid(T,s,halpha,freqs,harms,ParamRange,dParam,df,Method)
	
	Xax = ParAx[PlotAx[0]]
	Yax = ParAx[PlotAx[1]]
	Xc = ParC[PlotAx[0]]
	Yc = ParC[PlotAx[1]]
	nP = np.size(dParam)
	if nP > 2:
		Axes = np.zeros(nP,dtype='int32')
		Axes[PlotAx[0]] = 1
		Axes[PlotAx[1]] = 1
		if AxInds is None:
			AxInds = np.zeros(nP,dtype='int32')
		string='grid = grid['
		for i in range(0,nP):
			if Axes[i] == 0:
				string += '{:d}'.format(AxInds[i])
			else:
				string += ':'
			if i < nP-1:
				string+=','
		string +=']'
		eval(string)

	if PlotAx[0] < PlotAx[1]:
		grid = grid.T
			
	
	mn = np.where(grid == np.nanmin(grid))
	i = mn[0][0]
	j = mn[1][0]
	grid = np.log10(grid)
	if fig is None:
		fig = plt
		fig.figure()
	if scale is None:
		scale = [np.nanmin(grid),np.nanmax(grid)]
	#fig = PlotGrid(fig,grid,Xax,Yax,scale=scale,xtitle='log$_{10}(\\rho_{eq})$ (amu cm$^{-3}$)',ytitle='Power',maps=maps,cmap=plt.cm.get_cmap('gnuplot'))
	ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))

	xrnge=[np.nanmin(Xax),np.nanmax(Xax)]
	yrnge=[np.nanmin(Yax),np.nanmax(Yax)]
	fig.axis([xrnge[0],xrnge[1],yrnge[0],yrnge[1]])
	fig.xlabel('log$_{10}(\\rho_{eq})$ (amu cm$^{-3}$)')
	fig.ylabel('Power')

	cmap=plt.cm.get_cmap('gnuplot')

	xmsh,ymsh=np.meshgrid(Xax,Yax)
	grid=np.transpose(np.array(grid))
	grid=ma.masked_where(np.isnan(grid),grid)
	fig.pcolormesh(xmsh,ymsh,grid,cmap=cmap,vmin=scale[0],vmax=scale[1])
	
	ztitle = 'Error log$_{10}(\Delta f_{RMS})$ (mHz)' 
	if nobar == False:
		sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=scale[0], vmax=scale[1]))
		sm._A = []
		cbar=fig.colorbar(sm)
		cbar.set_label(ztitle)


	return fig

