import numpy as np
from Tracing import RK4 as rk4
from scipy.interpolate import InterpolatedUnivariateSpline

	
def GetDipoleModel(Beq):
	'''
	Create a function which provides magnetic field vector at a given positional vector.
	
	Args:
		Beq: Magnetic field strength at the equator in nT
		
	Returns:
		Function which takes positional vector p and returns a field vector.
	'''
	def DipoleField(p):
		#p must be a vector, or array of vectors which are in Rp
		mu = np.array([0.0,0.0,Beq]) #nT R_p^3
		M = np.linalg.norm(mu)
		mhat = np.array([mu/M])
		if np.size(np.shape(p)) == 1:
			r = np.linalg.norm(p)
			rhat = p/r
			mdotr = np.sum(mhat*rhat)
			B = M*((3.0*mdotr*rhat - mhat)/r**3.0)
		else:
			r = np.sqrt(p[:,0]**2.0 + p[:,1]**2.0 + p[:,2]**2.0)
			rhat = (p.T/r).T
			mdotr = np.sum(mhat*rhat,axis=1)
			B = M*((3.0*mdotr*rhat.T - mhat.T)/r**3.0).T
		return B
		
	return DipoleField
	
class TraceDipoleField(object):
	'''
		Provides the trace along a field line in dipole field model.
		
		Args:
			p: This is the initial starting positional vector p = np.array([x,y,z]).
			ds: Step size in Rp.
			Beq: Equatorial field strength.
			nmax: Maximum number of steps to be taken.
			bounds: Boundaries of the trace, should be a 2-element array or list with inner and outer radii in Rp,
					or a function which provides True when a positional vector is within the accepted bounds and False
					when the position is outside of those bounds.
		
		Returns:
			Class object with the following variables:
				x,y,z: Each stores an array of positions along the field line trace in Rp
				Bx,By,Bz: Each is an array of magnetic field components for each pointalong the field line.
				nstep: number of steps taken in the trace.
				Lshell: This is not exactly the L-shell, unless a perfect dipole is used - this is the radial
						distance of the magnetic equatorial footprint of the field line (where z == 0) in Rp.
				MltE: This is the local time of the equatorial footprint of the field line.
				
				NOTE: lshell and mlte are both NaN if the field line is open. 
			
	'''
	def __init__(self,p,ds,Beq,nmax=10000,bounds=[1.0,100.0]):
		BFn = GetDipoleModel(Beq)
		T = rk4.RK4Trace(np.array([p]),ds,BFn,n=nmax,bounds=bounds,direction='both')
		good = np.where(np.isfinite(T[:,0]))[0]
		T = T[good]
		
		self.x = T[:,0]
		self.y = T[:,1]
		self.z = T[:,2]
		B = BFn(T)

		self.Bx = B[:,0]
		self.By = B[:,1]
		self.Bz = B[:,2]
		
		self.nstep = good.size
		
		n = np.where(self.z > 0)[0]
		s = np.where(self.z < 0)[0]
		

		
		if n.size > 1  and s.size > 1:
			x = np.append(self.x[s[:-2]],self.x[n[:2]])
			y = np.append(self.y[s[:-2]],self.y[n[:2]])
			z = np.append(self.z[s[:-2]],self.z[n[:2]])
			fzx = InterpolatedUnivariateSpline(z,x)
			fzy = InterpolatedUnivariateSpline(z,y)
			e = np.array([fzx(0.0),fzy(0.0),0.0])
			self.Lshell = np.sqrt(e[0]**2.0 + e[1]**2.0)
			self.MltE = (-np.arctan2(e[1],-e[0])*180.0/(np.pi*15.0) + 24.0)%24.0
			
			if inspect.isfunction(bounds) == False:
				R = np.sqrt(self.x**2 + self.y**2 + self.z**2)
				fRx = InterpolatedUnivariateSpline(R[n[::-1]][0:4],self.x[n[::-1]][0:4])
				fRy = InterpolatedUnivariateSpline(R[n[::-1]][0:4],self.y[n[::-1]][0:4])
				fRz = InterpolatedUnivariateSpline(R[n[::-1]][0:4],self.z[n[::-1]][0:4])
			
			
				
				xn = fRx(1.0)
				yn = fRy(1.0)
				zn = fRz(1.0)

				fRx = InterpolatedUnivariateSpline(R[s[0:4]],self.x[s[0:4]])
				fRy = InterpolatedUnivariateSpline(R[s[0:4]],self.y[s[0:4]])
				fRz = InterpolatedUnivariateSpline(R[s[0:4]],self.z[s[0:4]])
				xs = fRx(1.0)
				ys = fRy(1.0)
				zs = fRz(1.0)
				

				gd = np.where(R > 1.0)[0]
				self.x = np.append(xs,np.append(self.x[gd],xn))
				self.y = np.append(ys,np.append(self.y[gd],yn))
				self.z = np.append(zs,np.append(self.z[gd],zn))

				B = BFn(np.array([self.x,self.y,self.z]).T)
				self.Bx = B[:,0]
				self.By = B[:,1]
				self.Bz = B[:,2]
				self.nstep = np.size(self.x)
		else:
			self.Lshell = np.nan
			self.MltE = np.nan
