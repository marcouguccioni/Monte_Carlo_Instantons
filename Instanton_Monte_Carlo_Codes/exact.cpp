#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

	//c----------------------------------------------------------------------------c
    	//c     direct diagonalization of quantum mechanical anharmonic oscillator.    c
    	//c----------------------------------------------------------------------------c
    	//c     hamiltonian m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                  c
    	//c----------------------------------------------------------------------------c
    	//c     harmonic oscillator: H_0=m/2(\dot x)^2+m/2*w^2x^2                      c
    	//c     perturbation:        H_1=A*x^4+B*x^2+C                                 c
    	//c----------------------------------------------------------------------------c
    	
    	
    	
// ---------------------------------------------------------------------------------------  	
void jacobi(double **&a, double *&d, double **&v, int n)
{
	//c----------------------------------------------------------------------------c
	//c     diagonalize real symmetric n*n matrix a using jacobi method.           c
    	//c----------------------------------------------------------------------------c
    	//c     a(n,n)    real, symmetric n*n matrix                                   c
    	//c     n         dimension of matrix                                          c
    	//c     d(n)      eigenvalues                                                  c
    	//c     v(n,n)    eigenvectors v(i,1),...,v(n,i)                               c
    	//c----------------------------------------------------------------------------c
    	//c     note that the matrix a(n,n) is destoyed !                              c
    	//c----------------------------------------------------------------------------c
    	int nrot; 
    	double b[1000],z[1000];
    	double sm, tresh, g, h, t, theta, c, s, tau;
    	
    	
    	for(int p=0; p<n; p++)
    	{
        	for(int q=0; q<n; q++)
        	{
            	v[p][q] = 0.;
        	}
        	v[p][p] = 1.;
    	}
    	for(int p=0; p<n; p++)
    	{
        	b[p] = a[p][p];
        	d[p] = b[p];
        	z[p] = 0.;
    	}
    	nrot = 0;
    	for(int i=0; i<50; i++)
    	{
        	sm = 0.;
        	for(int p = 0; p<n-1; p++)
        	{
        		for(int q = p+1; q<n; q++)
        		{
                		sm += abs(a[p][q]);
                	}
        	}
        	if (sm == 0.) return;
        	if (i < 4) tresh = 0.2*sm/pow(n,2);
        	else tresh = 0.;
        	for(int p = 0; p < n-1; p++)
        	{
            		for(int q = p+1; q<n; q++)
            		{
             			g = 100.*abs(a[p][q]);
                		if (((i > 4)) && ((abs(d[p]) + g == abs(d[p]))) && ((abs(d[q]) + g == abs(d[q]))))
                		{
                    			a[p][q] = 0.;
                		} 
                		else if (abs(a[p][q]) > tresh)
                		{
                    			h = d[q] - d[p];
                    			if (abs(h) + g == abs(h))
                    			{
                        			t = a[p][q]/h;
                    			} 
                    			else
                    			{
                        			theta = 0.5*h/a[p][q];
                        			t = 1./(abs(theta) + sqrt(1. + pow(theta,2)));
                        			if (theta < 0.) t = -t;
                    			}
                    			c = 1./sqrt(1 + pow(t,2));
                    			s = t*c;
                    			tau = s/(1. + c);
                    			h = t*a[p][q];
                    			z[p] -= h;
                    			z[q] += h;
                    			d[p] -= h;
                    			d[q] += h;
                    			a[p][q] = 0.;
                    			for(int j = 0; j<p-1; j++)
                    			{
                        			g = a[j][p];
                        			h = a[j][q];
                        			a[j][p] = g - s*(h + g*tau);
                        			a[j][q] = h + s*(g - h*tau);
                        
                    			}
                    			for(int j = p+1; j<q-1; j++)
                    			{
                        			g = a[p][j];
                        			h = a[j][q];
                        			a[p][j] = g - s*(h + g*tau);
                        			a[j][q] = h + s*(g - h*tau);
                        
                    			}
                    			for(int j = q+1; j<n; j++)
                    			{
                        			g = a[p][j];
                        			h = a[q][j];
                        			a[p][j] = g - s*(h + g*tau);
                        			a[q][j] = h + s*(g - h*tau);
                        
                    			}
                    			for(int j=0; j<n; j++)
                   			{		
                        			g = v[j][p];
                        			h = v[j][q];
                        			v[j][p] = g - s*(h + g*tau);
                        			v[j][q] = h + s*(g - h*tau);
                    			}
                    			nrot = nrot + 1;
                		}
			}
		}
        	for(int p=0; p<n; p++)
        	{
            	b[p] += z[p];
            	d[p] = b[p];
            	z[p] = 0.;
            	}
	}
    	cout<< "too many iterations in jacobi";
    	return;
}

//--------------------------------------------------------------------------------------
void eigensort(double *&d,double **&v,int n)
{
	//c----------------------------------------------------------------------------c
    	//c     order eigenvectors and eigenvalues found by jacobi.                    c
    	//c----------------------------------------------------------------------------c
    	int k;
    	double p;
    	
    	for(int i = 0; i<n-1; i++)
    	{
        	k = i;
        	p = d[i];
        	for(int j= i+1; j<n; j++)
        	{
            		if (d[j] <= p)
            		{
                		k = j;
                		p = d[j];
            		}	
            
        	}
        	if (k != i)
        	{
            		d[k] = d[i];
            		d[i] = p;
            		for(int j =0; j<n; j++)
            		{
                		p = v[j][i];
                		v[j][i] = v[j][k];
                		v[j][k] = p;
                	}
        	}
        
    	}
    	return;
}

//-------------------------------------------------------------------------------------
void hermite(int n, double x, double *&p)
{
    	//c----------------------------------------------------------------------------c
    	//c     rescaled hermite polynomials H_i(x)/2^i/sqrt(i!)                       c
    	//c----------------------------------------------------------------------------c
    	p[0] = 1.0;
    	p[1] = x;
    	for(int i=2; i<n; i++)
    	{
        	p[i] = (x*p[i - 1] - sqrt(i - 1.0)*p[i - 2]/2.0)/sqrt(double(i));
        }
    	return;
}

//-------------------------------------------------------------------------------------
void psi_oscillator(int n, double x, double *&psi, double m, double w)
{
	//c----------------------------------------------------------------------------c
    	//c     harmonic oscillator wave functions psi(i=0,..,n)=psi_i(x)              c
    	//c----------------------------------------------------------------------------c
        double *h=new double[1000];
    	double y = sqrt(m*w)*x;
    	double norm;
        hermite(n, y, h);
    	for(int i=0; i<n; i++)
    	{
		norm = pow((m*w/M_PI),0.25)*pow(2.0,(i/2.0));
        	psi[i] = norm*h[i]*exp(-m*w/2.0*pow(x,2));
        }
    	return;
}

//-------------------------------------------------------------------------------------
int main() 
{

	int nmax=1000;
	double **h=new double *[nmax], *e=new double[nmax];
	double *psi=new double [nmax];
	double **v=new double *[nmax];
	double *rho=new double[nmax], *rho2=new double[nmax],*rho3=new double[nmax];
	double xcor,x2cor,x3cor;
	double m,w;
	
	for(int i=0; i<nmax; i++)
	{
		h[i]=new double [nmax];
		v[i]=new double [nmax];
		
	}
	
	double f; 
	int ndim;
	double w0;
	
	cout<<"parameter (f=1.4): ", cin>>f, cout<<endl;
	cout<<"dimension of matrix (ndim=40): ", cin>>ndim, cout<<endl;
	cout<<"unperturbed oscillator frequency w0 (4*f=5.6): ", cin>>w0, cout<<endl;
	
	double taumax = 2.5;
    	double ntau = 100;
    	double dtau = taumax/double(ntau);
    	
    	double xmax = 2.0*f;
    	double nx = 100;
    	double dx = 2.0*xmax/double(nx);
    	
    	cout<<"exact diagonalization"<<endl;
     	cout<<"---------------------------"<<endl;
    	cout<<" f = "<<f<<" ndim = "<<ndim<<endl; 
    	cout<<" taumax = "<<taumax<<" ntau = "<<ntau<<endl;
   	cout<<" xmax = "<<xmax<<" nx = "<<nx<<endl;
   	cout<<"---------------------------"<<endl;
   	
   	//c----------------------------------------------------------------------------c
    	//c     initialize parameters, clear arrays                                    c
    	//c----------------------------------------------------------------------------c
    	
    	for(int i=0; i<nmax; i++)
    	{
        	for(int j=0; j<nmax; j++)
        	{
            		h[i][j] = 0.0;
            	}
    	}
    
    	m = 0.5;
    	w = w0;
    
    	double a = 1.0;
    	double b = -2.0*pow(f,2) - m*pow(w,2)/2.0;
    	double c = pow(f,4);
        double cw = 1.0/sqrt(m*w);
    
    	double c22 = pow(cw,2)/2.0;
    	double c44 = pow(cw,4)/4.0;
    	
    	//c----------------------------------------------------------------------------c
    	//c     build up h                                                             c
    	//c----------------------------------------------------------------------------c
    	
    	double x4, x2, e0, hh;
    	for(int n=0; n<ndim; n++)
    	{
    		//c----------------------------------------------------------------------------c
        	//c     <n|h|n>                                                                c
        	//c----------------------------------------------------------------------------c
         	x4 = c44*3.0*(pow((n+1),2)+pow(n,2));
         	x2 = c22*(2*n+1);
         	e0 = w*(n+0.5) + c;
         
         	h[n][n] = a*x4 + b*x2 + e0;
         	
         	//c----------------------------------------------------------------------------c
        	//c     <n|h|n+2>                                                              c
        	//c----------------------------------------------------------------------------c
         
        	x4 = c44*sqrt((n + 1.0)*(n + 2))*(4*n + 6);
        	x2 = c22*sqrt((n + 1.0)*(n + 2));
         
        	hh = a*x4 + b*x2;
        	h[n][n + 2]= hh;
        	h[n + 2][n] = hh;
         
        	//c----------------------------------------------------------------------------c
        	//c     <n|h|n+4>                                                              c
        	//c----------------------------------------------------------------------------c
         
        	x4 = c44*sqrt((n + 1.0)*(n + 2)*(n + 3)*(n + 4));
         
        	hh = a*x4;
        	h[n][n + 4] = hh;
       	h[n + 4][n] = hh; 
        }
        
        //c----------------------------------------------------------------------------c
    	//c     diagonalize                                                            c
    	//c----------------------------------------------------------------------------c
    	
    	jacobi(h,e,v,ndim);
    	eigensort(e,v,ndim);
    	
    	//c----------------------------------------------------------------------------c
    	//c     energy eigenvalues and matrix elements <0|x|n>                         c
    	//c----------------------------------------------------------------------------c
    	
    	double dn,cn,en;
    	int kmax1,kmax2,kmax3,kmin1,kmin2,kmin3;
    	
    	cout<<"---------------------------"<<endl;
    	cout<<fixed<<"n"<<"\t"<<"e[n]"<<"\t\t"<<"rho[n]"<<"\t\t"<<"rho2[n]"<<"\t\t"<<"rho3[n]"<<endl;
    	for(int n = 0; n<ndim; n++)
    	{
        	cn = 0.0;
        	dn = 0.0;
        	en = 0.0;
        	for(int k = 0; k<ndim; k++)
        	{
            		kmax3 = max(k - 3, 0);
            		kmax2 = max(k - 2, 0);
            		kmax1 = max(k - 1, 0);
            		kmin1 = min(k + 1, ndim - 1);
            		kmin2 = min(k + 2, ndim - 1);
            		kmin3 = min(k + 3, ndim - 1);
            		cn +=  (sqrt(double(k))*v[kmax1][0] + sqrt(double(k + 1))*v[kmin1][0])*v[k][n];
            		dn += (sqrt(double(k*(k - 1)))*v[kmax2][0] + (2*k + 1)*v[k][0] + sqrt(double((k + 1)*(k + 2)))*v[kmin2][0])*v[k][n];
            		en += (sqrt(double(k*(k - 1)*(k - 2)))*v[kmax3][0] + 3*k*sqrt(double(k))*v[kmax1][0] + 3*(k + 1)*sqrt(double(k + 1))*v[kmin1][0] 
			+ sqrt(double((k + 1)*(k + 2)*(k + 3)))*v[kmin3][0])*v[k][n];
		}
        	rho[n] = pow(cw,2)/2.0*pow(cn,2);
        	rho2[n] = pow(cw,4)/4.0*pow(dn,2);
        	rho3[n] = pow(cw,6)/8.0*pow(en,2);
        	cout<<fixed<<n<<"\t"<<e[n]<<"\t"<<rho[n]<<"\t"<<rho2[n]<<"\t"<<rho3[n]<<endl;       
    	}
    	
    	//c----------------------------------------------------------------------------c
    	//c     groundstate wave function                                              c
    	//c----------------------------------------------------------------------------c
    	double psi0;
	double x;
	
	cout<<"---------------------------"<<endl;
	cout<<fixed<<"x"<<"\t\t"<<"psi(x)"<<"\t\t"<<"psi^2(x)"<<endl;
    	for(int k=0; k<nx; k++)
    	{
        	x = -xmax + k*dx;
        	psi0 = 0.0;
        	psi_oscillator(ndim, x, psi, m, w);
        	for(int j = 0; j<ndim; j++)
        	{
            		psi0 += v[j][0]*psi[j];
        	}
        	cout<<fixed<<x<<"\t"<<psi0<<"\t"<<pow(psi0,2)<<endl;
    	}
    	
    	//c----------------------------------------------------------------------------c
    	//c     x,x^2,x^3 correlators                                                  c
    	//c----------------------------------------------------------------------------c
    	e0 = e[0];
    	double tau;
    	double ej;
    	
    	cout<<"---------------------------"<<endl;
    	cout<<fixed<<"tau"<<"\t\t"<<"<x(0)x(tau)>"<<"\t"<<"<x^2(0)x^2(tau)>"<<"\t"<<"<x^3(0)x^3(tau)>"<<endl;
    	for(int k=0; k<ntau; k++)
    	{	
        	tau = k*dtau;
        	xcor = 0.0;
        	x2cor = 0.0;
        	x3cor = 0.0;
        	for(int j = 1; j<ndim; j++)
        	{
            		ej = e[j];
            		x2cor += rho2[j]*exp(-(ej - e0)*tau);
            		x3cor +=  rho3[j]*exp(-(ej - e0)*tau);
            		xcor +=  rho[j]*exp(-(ej - e0)*tau);
            	}
            	cout<<fixed<<tau<<"\t"<<xcor<<'\t'<<x2cor<<'\t'<<x3cor<<endl; 
       }
       
        //c----------------------------------------------------------------------------c
    	//c     partition function/free energy                                         c
    	//c----------------------------------------------------------------------------c
        double beta_max = 100.0;
    	double beta_min = 0.1;
    	double beta,log_beta,t,z,free;
    	double dlog = (log(beta_max) - log(beta_min))/double(50);
    	
    	cout<<"---------------------------"<<endl;
    	cout<<fixed<<"T"<<"\t\t"<<"beta"<<"\t\t"<<"F"<<endl;
    	for(int i=0; i<50; i++)
    	{
        	log_beta = log(beta_min) + i*dlog;
        	beta = exp(log_beta);
        	t = 1.0/beta;
       	z = 1.0;
        	for(int j = 0; j<ndim; j++)
        	{
            		z+=exp(-(e[j] - e[0])*beta);
            
        	}
        	free = t*log(z) - e[0];
        	cout<<fixed<<t<<"\t"<<beta<<"\t"<<free<<endl; 
        }



}



