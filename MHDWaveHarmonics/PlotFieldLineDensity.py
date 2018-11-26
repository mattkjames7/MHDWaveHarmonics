import numpy as np
import matplotlib.pyplot as plt
from .GetFieldLine import GetFieldLine

def PlotFieldLineDensity(pos,Params,Model='KT14',StepSize=None,ModelArgs=None,Delta=None,Polarization='none',fig=None,maps=[1,1,0,0],Overplot=False):
	'''
	Simple routine to plot the modelled plasma mass density along a field line.
	
	Args:
		pos: 3 element ndarray containing the cartesian position to trace from in Rp.
		Params: For power law: 2-element array/list [p_eq,power].
				For Sandhu model: 5-element array/list [n0,alpha,a,beta,mav0] (see GetSandhuParams).
		Model: 'KT14'|'T96'|'Dipole'.
		StepSize: Step length for field line trace.
		ModelArgs: Parameters to supply to field model, see GetFieldLine.
		Delta: Distance between field lines.
		Polarization: Wave polarization 'none'|'toroidal'|'poloidal'.
		fig: Pyplot instance, useful if plotting as a subplot on a pre-existing figure, by default will create new figure.
		maps: 4 element list containing subplot information [xmaps,ymaps,xmap,ymap] where xmaps and ymaps represent the number of subplots in the x and y directions and xmap and ymap are the specific x and y indices to plot to.
		Overplot: When set to True, this will inhibit the creation of a new subplot, fig must be supplied with a pre-existing figure, where the routine will plot over.
		
	Returns:
		pyplot instance
	
	'''
	tmp = GetFieldLine(pos,Model,StepSize,ModelArgs,Delta,Polarization)
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
	
