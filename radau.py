import sys
import numpy as np
from scipy import integrate
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.integrate import odeint
from scipy.integrate import quad, quad_vec, nquad
from scipy.optimize import minimize

i=sys.argv[1]
#print(i)
data_dir = sys.argv[2]
#print(data_dir)

N=np.load(data_dir+'/N.npy')
w=np.load(data_dir+'/w_'+str(i)+'.npy')
g=np.load(data_dir+'/g.npy')
[Nin,Nout,Ruin,dRuin]=np.load(data_dir+'/lims_'+str(i)+'.npy')

square_freq=np.square(w)
theta=np.real(square_freq)
theta_func=interp1d(N,theta)
phi=np.imag(square_freq)
phi_func=interp1d(N,phi)
g_func=interp1d(N,g)
def dx(x,n):
    du=x[0]
    u=x[1]
    dv=x[2]
    v=x[3]
    return [-2*g_func(n)*du-theta_func(n)*u+phi_func(n)*v,
            du,
            -2*g_func(n)*dv-theta_func(n)*v-phi_func(n)*u,
            dv]
def qx(n,x):
    return dx(x,n)
u0=np.real(Ruin)
v0=np.imag(Ruin)
du0=np.real(dRuin)
dv0=np.imag(dRuin)
x0=[du0,u0,dv0,v0]
NIntegrate=np.linspace(Nin,Nout,100)
#sol=odeint(dx,x0,NIntegrate)
#x1=sol[-1]
trim=1e-9
sol=integrate.solve_ivp(qx,(np.real(Nin+trim),np.real(Nout-trim)),x0,method='Radau',rtol=1e-10,atol=1e-13)
x1=sol.y[:,-1]
val1=x1[1]+1j*x1[3]
#print(val1)
np.save(data_dir+'/p_'+str(i),val1)
