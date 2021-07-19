#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>


using namespace std;

//c----------------------------------------------------------------------------c
//c     random instanton calculation in quantum mechanics.                     c
//c     in this version we add gaussian fluctuations to the classical          c
//c     instanton congiguration. this is done by performing a few monte        c
//c     carlo (heating) sweeps in the gaussian effective potential.            c
//c----------------------------------------------------------------------------c
//c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
//c----------------------------------------------------------------------------c


//------------------------------------------------------------------------
void instanton_location(double *&tau,int n_insta,double tau_max)
{
//c----------------------------------------------------------------------c
//c     initialize instanton configuration                               c
//c----------------------------------------------------------------------c

	for(int i=0; i<n_insta; i++)
	{
	double rnd =((double) rand())/((double)RAND_MAX);
	tau[i]=rnd*tau_max;
	}
	sort(tau,tau+n_insta); //include algorithm
	return;
}

//-------------------------------------------------------------------
double x_sum_ansatz( double *tau,int n_insta, double f, double t)
{
//c----------------------------------------------------------------------c
//c     sum ansatz path                                                  c
//c----------------------------------------------------------------------c
	double x_sum=-f;
	for(int i=0; i<n_insta-1; i++)
	{
		x_sum+= f*tanh(2.0*f*(t-tau[i])) - f*tanh(2.0*f*(t-tau[i+1]));
	}
	return x_sum;
}



// -----------------------------------------------------------------	
int main()
{
	srand(time(NULL)); //set seed
	double rnd;
	
	double *x=new double[1000];
	double *x_hot=new double[1000];
	double *w=new double[1000];
	double *tau= new double[1000];
	double x_cor,x2_cor,x3_cor;
	double *x_cor_sum=new double[100];
	double *x_cor_av=new double[100];
	double *x2_cor_sum=new double[100]; 
	double *x2_cor_av=new double[100];
	double *x3_cor_sum=new double[100];
	double *x3_cor_av=new double[100];
	
	double f,a;
	int n;
	
	ofstream path;
	
	path.open("heated_path.dat");
	
	cout<<"separation of wells f (f=1.4): ", cin>>f, cout<<endl;
	cout<<"lattice sites n<1000 (n=800): ", cin>>n, cout<<endl;
	cout<<"lattice spacing a (a=0.05): ", cin>>a, cout<<endl;
	
	double tau_max=n*a;
	double action0=4.0/3.0*pow(f,3);
	int n_insta= int( tau_max*8*sqrt(2.0/M_PI)*pow(f,2.5)*exp(-action0-71.0/72.0/action0) +0.5);
	int  n_mc, n_p ,n_c, k_p, n_heat;
	double delta_x;
	
	cout<<"number of instantons (N="<<n_insta<<" even): ", cin>>n_insta, cout<<endl;
	cout<<"monte carlo sweeps (10^4): ", cin>>n_mc, cout<<endl;
	cout<<"update x (delta_x=0.5): ", cin>>delta_x, cout<<endl;
	cout<<"number of points in correlator (n_p=20): ", cin>>n_p, cout<<endl;
	cout<<"number of correlator measurments per config (n_c=5): ", cin>>n_c, cout<<endl;
	cout<<"write every k-th configuration: (k_p=100): ", cin>>k_p, cout<<endl;
	cout<<"heating sweeps (n_heat=10): ", cin>>n_heat, cout<<endl;
	
	
	cout<<"Heating method"<<endl;
     	cout<<"---------------------------"<<endl;
    	cout<<" f = "<<f<<" n = "<<n<<" a = "<<a<<endl; 
    	cout<<" n_mc = "<<n_mc<<" n_insta = "<<n_insta<<endl;
   	cout<<" n_p = "<<n_p<<" n_c = "<<n_c<<endl;
   	cout<<" delta_x = "<<delta_x<<" n_heat = "<<n_heat<<endl;
   	cout<<"---------------------------"<<endl;
   	
   	//c----------------------------------------------------------------------------c
	//c     clear summation arrays                                                 c
	//c----------------------------------------------------------------------------c
	for(int k=0; k<100; k++)
        {
        	x_cor_sum[k]=0.0;
        	x2_cor_sum[k]=0.0;
        	x3_cor_sum[k]=0.0;
        }	
        
        //c----------------------------------------------------------------------------c
	//c     loop over configurations                                               c
	//c----------------------------------------------------------------------------c
	
	int n_accepted = 0;
	int n_hit = 0;
	int n_config = 0;
	int n_cor = 0;
	
	double x_new;
	double s_old,s_new,delta_s;
	double t;
	int p0;
	
	for (int i=0; i<n_mc; i++)
	{
		n_config += 1;
		instanton_location( tau,n_insta, tau_max);

		//c----------------------------------------------------------------------------c
		//c     new configuration                                                      c
		//c----------------------------------------------------------------------------c
		for( int j=0; j<n; j++)
		{
			t=(double) j*a;
            		x[j] = x_sum_ansatz(tau,n_insta,f,t);
           	}
           	x[n-1] = x[0];
		x[n] = x[1];   
		
		//c----------------------------------------------------------------------------c
		//c     heat configuration: start from classical path                          c
		//c----------------------------------------------------------------------------c
		for(int k=0; k<n; k++)
		{
			x_hot[k] = x[k];
			w[k]=-4.0*(pow(f,2)-3.0*pow(x[k],2));
		}

		//c----------------------------------------------------------------------------c
		//c     heating sweeps                                                         c
		//c----------------------------------------------------------------------------c

         	for(int i_h=0; i_h<n_heat; i_h++)
         	{
         		for(int j=1; j<n; j++)
         		{	
         			rnd=((double) rand())/((double)RAND_MAX);
         			s_old=a*(1.0/4.0*(pow((x_hot[j]-x_hot[j-1])/a,2)+pow((x_hot[j+1]-x_hot[j])/a,2)) +0.5*w[j]*pow(x_hot[j]-x[j],2));
         			x_new = x_hot[j] + delta_x*(2*rnd-1);
         			s_new=a*(1.0/4.0*(pow((x_new-x_hot[j-1])/a,2)+pow((x_hot[j+1]-x_new)/a,2)) +0.5*w[j]*pow(x_new-x[j],2));
           			delta_s = s_new - s_old;
           			rnd=((double) rand())/((double)RAND_MAX);
        			delta_s = min(delta_s, 70.0);
        			delta_s = max(delta_s, -70.0);
        			n_hit+=1;
         			if (exp(-delta_s)  > rnd) 
        			{
            				x_hot[j] = x_new;
            				n_accepted += 1;
        			}
         		}
         		x_hot[n-1] = x_hot[0];
			x_hot[n] = x_hot[1];
		}
		//c----------------------------------------------------------------------------c
		//c     configuration                                                          c
		//c----------------------------------------------------------------------------c
		if ((i % k_p)  ==  0) 
    		{
    			for(int k=0; k<n; k++)
        		{
            			path<<fixed<<k*a<<"\t"<<x[k]<<"\t"<<x_hot[k]<<"\t\t"<<n_config<<endl;
            		}
        	}
        	
        	//c----------------------------------------------------------------------------c
    		//c     correlation function                                                   c
    		//c----------------------------------------------------------------------------c
    
    		for(int c=0; c<n_c; c++)
    		{
         		n_cor += 1;
         		rnd=((double) rand())/((double)RAND_MAX);
        		p0 = (int)((n - n_p)*rnd);
         		for(int p=0; p<n_p; p++)
        		{
				x_cor = x_hot[p0]*x_hot[p0+p];
            			x2_cor = pow(x_cor,2);
            			x3_cor = pow(x_cor,3);
             			x_cor_sum[p] += x_cor;
            			x2_cor_sum[p] += x2_cor;
            			x3_cor_sum[p] += x3_cor;
            		}
         	}
         	//c----------------------------------------------------------------------------c
    		//c     next configuration                                                     c
    		//c----------------------------------------------------------------------------c
  	}
  	//c----------------------------------------------------------------------------c
	//c     averages                                                               c
	//c----------------------------------------------------------------------------c
	for(int p=0; p<n_p; p++)
	{	
		x_cor_av[p]=x_cor_sum[p]/(double)n_cor;
    		x2_cor_av[p]=x2_cor_sum[p]/(double)n_cor;
    		x3_cor_av[p]=x3_cor_sum[p]/(double)n_cor;
    		
    	}
    	//c----------------------------------------------------------------------------c
	//c     x,x^2,x^3 correlation functions                                        c
	//c----------------------------------------------------------------------------c
	
	cout<<"------------------------------------------------"<<endl;
	cout<<fixed<<"tau"<<"\t\t"<<"<x(0)x(tau)>"<<"\t"<<"<x^2(0)x^2(tau)>"<<"\t"<<"<x^3(0)x^3(tau)>"<<endl;
	for(int p = 0; p<n_p; p++)
	{
   		cout<<fixed<<p*a<<'\t'<<x_cor_av[p]<<'\t'<<x2_cor_av[p]<<'\t'<<x3_cor_av[p]<<endl;  
    	}   	
   	
}
