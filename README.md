# MHDWaveHarmonics
Some simple tools for modelling MHD wave harmonics in an arbitrary magnetic field geometry.

## Installation

Install using `pip3`:

```
pip3 install MHDWaveHarmonics --user
```

The installation will require the following packages:

* `numpy`
* `scipy`
* `matplotlib`
* `PyGeopack`
* `FieldTracing`

all of which will be installed automatically. For modelling waves in
Mercury's magnetosphere you will also require the `KT17` model.

## Usage

### 1. `GetFieldLine`

The `GetFieldLine` function will trace a model field, returning a 
`TraceField` object alongside a `numpy.ndarray`, `s`, which contains the 
distance along the traced field line and optionally `h` if the 
polarization is specified:

```python
import MHDWaveHarmonics as wh
T,s = wh.GetFieldLine(pos,Date=None,ut=None,Model='KT17',Delta=None,Polarization='none',**kwargs)
T,s,h = wh.GetFieldLine(pos,Date=None,ut=None,Model='KT17',Delta=None,Polarization='poloidal',**kwargs)
```

where `h` is provided using a second trace if the `Polarization` parameter
is set to `'poloidal'`, `'toroidal'` or a `float` corresponding to an
angle in degrees anticlockwise from the poloidal direction (the outward
pointing normal direction of the field line). Use ```wh.GetFieldLine?```
to find out more about this procedure from its docstring.

### 2. `SolveWave`

This procedure will attempt to solve the wave equation along an arbitrary
magnetic field trace such as that provided by `GetFieldLine`:

```python
yr = wh.SolveWave(f,x,B,R=None,Va=None,halpha=None,Params=None,InPlanet=None,Method='Simple',Unscale=True)
yr,yi,phase = SolveWave(f,x,B,R=None,Va=None,halpha=None,Params=None,InPlanet=None,Method='Complex',Unscale=True)
```

where `yr` is the real part of the scaled field displacement, &#x03BE;/h<sub>&#x03B1;</sub> , `yi` 
is the imaginary part, and `phase` is a continuous measure of the waves
phase along the trace. See the docstring for more information.

### 3. `FindHarmonics`

This will attempt to solve the wave equation in order to find a number
of harmonics which would be allowed to exist:

```python
freq,success,niter = wh.FindHarmonics(T,s,Params,halpha=None,Harmonics=[1,2,3],x0=None,df=1.0,Method='Complex')
```
where `freq` is an array of allowed frequencies, `success` is a Boolean 
array denoting whether each fit was successful or not and `niter` is an 
array containing the number of iterations reqquired for each fit.

### 4. `CalcFieldLineVa`

This will calculate the Alfven speed at each point along the trace:

```python
va = CalcFieldLineVa(T,s,Params,halpha=None)
```

or at the midpoint between each pair of consecutive steps along the trace:

```python
vamid = CalcFieldLineVaMid(T,s,Params,halpha=None)
```

### 5. `PlotHarmonics`

This will produce a plot of the allowed toroidal and poloidal harmonics 
on a field line given an initial position for the trace and a plasma
mass density profile along the field.

```python
PlotHarmonics(pos,Params,nh=3,df=1.0,Rp=1.0,Colours=None,Method='Complex',**kwargs)
```
### 6. `FitPlasmaToHarmonic`

Attempts to fit a power law or the Sandhu et al model plasma mass density
profile to a given field trace.

```python
p_eq = FitPlasmaToHarmonic(T,s,halpha,f,Params,Harm=1,df=1.0,Method='Complex')
```

### 7. `FitPlasma`

This tries to fit the plasma mass density profile to a set of observed
frequencies (with their assumed harmonic numbers) along a field trace.

```python
params = FitPlasma(T,s,halpha,freqs,harms,Params0,df=1.0,Method='Complex',ParamFit=None)
```

### 8. `GetSandhuParams`

Calculates the parameters for the Sandhu et al electron density and 
average ion mass models.

```python
ne0,alpha,a,beta,mav0 = GetSandhuParams(mlt,L)
```

### 9. `PlotFieldLineDensity`

Plots the modelled density along a field line given a position in which
to start a field trace and a set of plasma profile parameters.

```python
PlotFieldLineDensity(pos,Params,fig=None,maps=[1,1,0,0],Overplot=False,**kwargs)
```


