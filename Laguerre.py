import random
import time
import numpy as np
start = time.time()

def MAX(a,b):
    if(a>b): return(a)
    return(b)

appo=1.0e-5; 
Nmax=10000

poly = [c0,c1,c2,c3,c4,c5]
N=len(poly)-1### 5
qoly= poly

def poolish(cof, x0):
    n = len(cof)
    b=np.zeros(MAX(1, n-1), dtype=complex)
    for i in range(n-1, 0, -1):
        if(i<n-1):   b[i-1] = cof[i] +  b[i]*x0;
        else:          b[i-1] =cof[i]  ;     
    return(b);
def deriv(cof):
    n=len(cof);
    b=np.zeros(MAX(1, n-1), dtype=complex)
    for i in range(n-2):
        b[i]= cof[i+1]*complex(i+1,0.0)  
    return(b)   
#########################
def evalu(cof, x0):
    n= len(cof); 
    b=np.zeros(MAX(1, n-1), dtype=complex)
    for i in range(n-1, 0, -1):
        if(i<n-1): b[i-1]= cof[i] + b[i]*x0;
        else:        b[i-1]=cof[i]    
    return( cof[0] + b[0] * x0)
def Froots(cof , x0):
    n = len(cof)-1;
    nite=0; 
    Der1 =[]
    Der2 =[]
    while(True): 
        pol= evalu(cof, x0);
        Der1= deriv(cof);
        Der2= deriv(Der1); 
        der1= evalu(Der1, x0); 
        der2= evalu(Der2, x0); 
        if(abs(pol)<float(appo)):  break;
        G= der1/ pol; 
        H= G*G - der2 - pol; 
        dd= np.sqrt( complex(n-1.0, 0.0)*(H*complex(n, 0.0) -G*G)); 
        A1= (G+ dd); 
        A2= (G- dd);
        if(abs(A2)>abs(A1)):  A1=A2;
        A1=complex(n, 0.0)/A1;     
        if(A1==0.0): break; 
        x0 -= A1; 
        nite+=1; 
        if(abs(A1)<appo or nite>Nmax):  break    
    return(x0)
##########################
roots=np.zeros(N,  dtype=complex)
ni=0; 
while(len(qoly)>2):
    xr=random.uniform(0,0.5)
    xi=random.uniform(0,0.5)
    x0=complex(xr,xi)
    x0=Froots(qoly,x0) 
    x0=Froots(poly,x0) 
    qoly=poolish(qoly,x0)
    roots[ni]=x0
    ni+=1;
roots[ni]=(-qoly[0]/qoly[1]) 
print (roots)  
T2.append(time.time() - start)