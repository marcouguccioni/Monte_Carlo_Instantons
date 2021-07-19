#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>


using namespace std;


	//c----------------------------------------------------------------------------c
	//c     lattice calculation in quantum mechanics.                              c
	//c----------------------------------------------------------------------------c
	//c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
	//c----------------------------------------------------------------------------c
	//c     lattice x(i=0,n), periodic b.c. x(0)=x(n)                               c
	//c----------------------------------------------------------------------------c
	//c     n_mc    number of mc sweeps                                             c
	//c     n_c     number of correlator measurements in a single configuration     c
	//c     n_p     number of points on which correlators are measured              c
	//c     n_eq    number of equlibration sweeps                                   c
	//c     k_p     number of sweeps between writeout of complete configuration     c
	//c----------------------------------------------------------------------------c

	
// ---------------------------------------------------------------------------
void histo_array(double &a,double a_min, double delta, int n,int  *&hist)
{
    	//c------------------------------------------------------------------------c
    	//c     include value a in histogram array ist(n)                          c
    	//c------------------------------------------------------------------------c
    	//c     a      value to be added to histogram array                        c
    	//c     a_min   minimum value in histogram                                 c
    	//c     delta     bin width                                                c
    	//c     n      number of bins                                              c
    	//c     hist(n) histogram array                                            c
    	//c------------------------------------------------------------------------c
    
    	//DIMENSION ist(150)
    	int j = (a - a_min)/delta + 1.000001;
    	if (j < 1)j = 1;
    	if (j > n)j = n;
    	hist[j] += 1;
    	return;
}

	
// -----------------------------------------------------------------	
int main()
{
	srand(time(NULL)); //set seed
	double rnd;
	
	double *x=new double[1000];
	int *histo_x=new int[100];
	double x_cor,x2_cor,x3_cor;
	double *x_cor_sum=new double[100];
	double *x_cor_av=new double[100];
	double *x2_cor_sum=new double[100]; 
	double *x2_cor_av=new double[100];
	double *x3_cor_sum=new double[100];
	double *x3_cor_av=new double[100];
	
	double f,a;
	int n, n_eq, n_mc, n_p ,n_c, k_p;
	bool hot;
	double delta_x, tau_max;
	
        
	ofstream action,path;
	
	action.open("action.dat");
	path.open("euclidean_path.dat");
	
	cout<<"separation of wells f (f=1.4): ", cin>>f, cout<<endl;
	cout<<"lattice sites n<1000 (n=800): ", cin>>n, cout<<endl;
	cout<<"lattice spacing a (a=0.05): ", cin>>a, cout<<endl;
	cout<<"cold/hot start (0,1): ", cin>>hot, cout<<endl;
	cout<<"equilibration sweeps before measurement (n_eq=100): ", cin>>n_eq, cout<<endl;
	cout<<"monte carlo sweeps (10^5): ", cin>>n_mc, cout<<endl;
	cout<<"update x (delta_x=0.5): ", cin>>delta_x, cout<<endl;
	cout<<"number of points in correlator (n_p=20): ", cin>>n_p, cout<<endl;
	cout<<"number of correlator measurments per config (n_c=5): ", cin>>n_c, cout<<endl;
	cout<<"write every k-th configuration: (k_p=1000): ", cin>>k_p, cout<<endl;
	
	tau_max=n*a;
	
	cout<<"Monte Carlo on lattice"<<endl;
     	cout<<"---------------------------"<<endl;
    	cout<<" f = "<<f<<" n = "<<n<<" a = "<<a<<endl; 
    	cout<<" n_mc = "<<n_mc<<" n_eq = "<<n_eq<<endl;
   	cout<<" n_p = "<<n_p<<" n_c = "<<n_c<<endl;
   	cout<<" delta_x = "<<delta_x<<" hot = "<<hot<<endl;
   	cout<<"---------------------------"<<endl;
   	
   	
   	//c----------------------------------------------------------------------------c
	//c     parameters for histograms                                              c
	//c----------------------------------------------------------------------------c
   	
   	int n_bins = 50;
	double x_hist_min = -1.5*f;
	double delta_bin = 3.0*f/double(n_bins);
	
	//c----------------------------------------------------------------------------c
	//c     clear summation arrays                                                 c
	//c----------------------------------------------------------------------------c
	for(int k=0; k<100; k++)
        {
        	x_cor_sum[k]=0.0;
        	x2_cor_sum[k]=0.0;
        	x3_cor_sum[k]=0.0;
        	histo_x[k]=0.0;
        }	
	
	
        //c----------------------------------------------------------------------------c
	//c     initialize                                                             c
	//c----------------------------------------------------------------------------c
	
	
	for(int i=0; i<n; i++)
	{
		
		rnd=((double) rand())/((double)RAND_MAX); //Uniform(0,1)
    		if (hot  == 0) 
    		{
        		x[i] = -f;
    		} 
    		else
    		{
        		x[i] = f*(2*rnd-1); //Uniform(-f,f)
    		}
    	}
    	
    	//c----------------------------------------------------------------------------c
	//c     periodic boundary conditions                                           c
	//c----------------------------------------------------------------------------c
	x[n-1] = x[0];
	x[n] = x[1];
	
	//c----------------------------------------------------------------------------c
	//c     initial action                                                         c
	//c----------------------------------------------------------------------------c
	double s=0,s0=0;
	
	for(int i=0; i<n; i++)
	{
		s0=(1.0/4.0*pow((x[i+1]-x[i])/a,2)+pow(pow(x[i],2)-pow(f,2),2))*a;
    		s += s0;
    	}
    	
    	//c----------------------------------------------------------------------------c
	//c     monte carlo sweeps                                                     c
	//c----------------------------------------------------------------------------c
	
	int n_accepted = 0;
	int n_hit = 0;
	int n_config = 0;
	int n_cor = 0;
	
	double x_new;
	double s_old,s_new,delta_s;
	int p0;

	
	action<<fixed<<"config"<<"\t\t"<<"acceptance rate"<<"\t\t"<<"S"<<endl;
	path<<fixed<<"tau"<<"\t\t"<<"x[tau]"<<"\t\t"<<"config"<<endl;
    			
	
	for(int i=0; i<n_mc; i++)
	{
		n_config += 1;
		if(i==n_eq)
		{
        		n_config = 0;
        		n_cor = 0;
        		for(int k=0; k<100; k++)
        		{
        		 	x_cor_sum[k]=0.0;
        		 	x2_cor_sum[k]=0.0;
        		 	x3_cor_sum[k]=0.0;
        		 	histo_x[k]=0.0;
        		 }	
    		}
    		//c----------------------------------------------------------------------------c
    		//c     one sweep thorough configuration                                       c
    		//c----------------------------------------------------------------------------c
    		
    		for(int j = 1; j<n; j++)
    		{
    			rnd=((double) rand())/((double)RAND_MAX);
         		n_hit +=  1;
         		x_new = x[j] + delta_x*(2*rnd-1);  //Uniform(x[j]-deltax,x[j]+deltax)
         		s_old=a*(1.0/4.0*(pow((x[j]-x[j-1])/a,2)+pow((x[j+1]-x[j])/a,2)) +pow((pow(x[j],2) - pow(f,2)),2));
         		s_new=a*(1.0/4.0*(pow((x_new-x[j-1])/a,2)+pow((x[j+1]-x_new)/a,2)) +pow((pow(x_new,2) - pow(f,2)),2));
        		rnd=((double) rand())/((double)RAND_MAX);
        		delta_s = s_new - s_old;
        		delta_s = min(delta_s, 70.0);
        		delta_s = max(delta_s, -70.0);
         		if (exp(-delta_s)  > rnd) 
        		{
            			x[j] = x_new;
            			n_accepted += 1;
        		}
         	}
         	x[n-1] = x[0];
		x[n] = x[1];
		
		//c----------------------------------------------------------------------------c
    		//c     calculate action etc.                                                  c
    		//c----------------------------------------------------------------------------c
    		s=0;
    		for(int j=0; j<n; j++)
		{
			s0=(1.0/4.0*pow((x[j+1]-x[j])/a,2)+pow(pow(x[j],2)-pow(f,2),2))*a;
    			s += s0;
    		}
    		
    		if ((i % k_p)  ==  0) 
    		{
    			action<<n_config<<'\t'<<double(n_accepted)/double(n_hit)<<'\t'<<s<<endl;
    			for(int k=0; k<n; k++)
        		{
            			path<<fixed<<k*a<<"\t"<<x[k]<<"\t\t"<<n_config<<endl;
            		}
        	}
        	//c----------------------------------------------------------------------------c
    		//c     populate histogram                                                     c
    		//c----------------------------------------------------------------------------c
    		
    		for(int k=0; k<n; k++)
    		{
        		histo_array(x[k], x_hist_min, delta_bin, n_bins, histo_x);
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
				x_cor = x[p0]*x[p0+p];
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
    	
    	//c----------------------------------------------------------------------------c
	//c     ground state wave function histogram                                   c
	//c----------------------------------------------------------------------------c
	cout<<"---------------------------------------------"<<endl;
	cout<<"x"<<"\t\t"<<"P(x)"<<endl;
	double x_norm = 0,xx;
	for(int i=0; i<n_bins; i++)
	{
    		x_norm += histo_x[i]*delta_bin;
    	}
	for(int i=0; i<n_bins; i++)
	{
    		xx = x_hist_min + i*delta_bin;
    		cout<<fixed<<xx<<'\t'<< histo_x[i]/x_norm<<endl;
   	}


  
	
	
	
	
	
	
}
