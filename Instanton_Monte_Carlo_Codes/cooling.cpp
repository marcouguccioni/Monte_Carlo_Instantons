#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>


using namespace std;


//c----------------------------------------------------------------------------c
//c     lattice calculation in quantum mechanics.                              c
//c----------------------------------------------------------------------------c
//c     this version applies cooling to MC configurations.                     c
//c----------------------------------------------------------------------------c
//c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
//c----------------------------------------------------------------------------c
//c     lattice x(i=0,n), peridic b.c. x(0)=x(n)                               c
//c----------------------------------------------------------------------------c
//c     n_mc    number of mc sweeps                                            c
//c     n_c     number of correlator measurements in a single configuration    c
//c     n_p     number of points on which correlators are measured		  c
//c     n_eq    number of equlibration sweeps                                  c
//c     k_p     number of sweeps between cooling                               c
//c     n_cool  number of cooling sweeps in a single configuration             c
//c     k_p2    number of sweeps between writeout of complete configuration    c
//c----------------------------------------------------------------------------c


 void cooling(double *&x, int n , double f, double a, double delta_x)
 {
//c----------------------------------------------------------------------------c
//c     local cooling algorithm                                                c
//c----------------------------------------------------------------------------c
    
	double s_old,s_new;
	double x_new;
	double rnd;

      	for(int i=1; i<n; i++)
      	{ 

        	s_old=a*(1.0/4.0*(pow((x[i]-x[i-1])/a,2)+pow((x[i+1]-x[i])/a,2)) +pow((pow(x[i],2) - pow(f,2)),2));
		for(int j=0; j<10; j++)
		{
			rnd=((double) rand())/((double)RAND_MAX);
            		x_new = x[i] + delta_x*0.1*(2.0*rnd-1.0);
            		s_new=a*(1.0/4.0*(pow((x_new-x[i-1])/a,2)+pow((x[i+1]-x_new)/a,2)) +pow((pow(x_new,2) - pow(f,2)),2));
             		if (s_new < s_old) x[i] = x_new;
		}
	}
	return;
} 


//--------------------------------------------------------------------
void instanton(double *x,int n, int &n_i, int &n_a, double *&x_i, double *&x_a,double a)
{
//c------------------------------------------------------------------------c
//c     return number and location of (anti) instantons                    c
//c------------------------------------------------------------------------c
      	n_i = 0;
      	n_a = 0;
      
	int i_x;
	int i_xp;
	if (x[0]<0) i_x=-1;
	else i_x=1; 
	
 	for(int i=1; i<n; i++)
      	{ 
      		if(x[i]<0) i_xp=-1;
      		else i_xp=1;
         	if(i_xp > i_x) 
         	{            
         		n_i += 1;
            		x_i[n_i] = i*a;
            	}
         	else if(i_xp < i_x)
         	{
            		n_a += 1;
            		x_a[n_a] = i*a;
    		}
         	i_x = i_xp;
	}
      	return;
}


//-----------------------------------------------------------------------------------
int main()
{
	srand(time(NULL)); //set seed
	double rnd;
	
	double *x=new double[1000];
	double *x_cool=new double[1000];
	double *x_i=new double[1000];
	double *x_a= new double[1000];
	double *z=new double[1000];
	double x_cor,x2_cor,x3_cor;
	double *x_cor_sum=new double[100];
	double *x_cor_av=new double[100];
	double *x2_cor_sum=new double[100]; 
	double *x2_cor_av=new double[100];
	double *x3_cor_sum=new double[100];
	double *x3_cor_av=new double[100];
	
	int n_i,n_a;
	int n_insta;
	double  *n_insta_sum=new double[1000];
	double *n_insta_av=new double[1000];
	
	
	double f,a;
	int n, n_eq, n_mc, n_p ,n_c, k_p, n_cool, k_p2;
	bool hot;
	double delta_x, tau_max;
	
	ofstream path;
	
	path.open("cooled_path.dat");
	
	cout<<"separation of wells f (f=1.4): ", cin>>f, cout<<endl;
	cout<<"lattice sites n<1000 (n=800): ", cin>>n, cout<<endl;
	cout<<"lattice spacing a (a=0.05): ", cin>>a, cout<<endl;
	cout<<"cold/hot start (0,1): ", cin>>hot, cout<<endl;
	cout<<"equilibration sweeps before measurement (n_eq=100): ", cin>>n_eq, cout<<endl;
	cout<<"monte carlo sweeps (10^4): ", cin>>n_mc, cout<<endl;
	cout<<"update x (delta_x=0.5): ", cin>>delta_x, cout<<endl;
	cout<<"number of points in correlator (n_p=20): ", cin>>n_p, cout<<endl;
	cout<<"number of correlator measurments per config (n_c=5): ", cin>>n_c, cout<<endl;
	cout<<"write every k-th configuration: (k_p2=1000): ", cin>>k_p2, cout<<endl;
	cout<<"number of sweeps between cooling (k_p=20): ", cin>>k_p, cout<<endl;
	cout<<"number of cooling sweeps (n_cool=100): ", cin>>n_cool, cout<<endl;
	
	tau_max=n*a;
	
	double action0=4.0/3.0*pow(f,3);
	double density=8*sqrt(2.0/M_PI)*pow(f,2.5)*exp(-action0-71.0/72.0/action0);
	
	cout<<"Monte Carlo on lattice with cooling"<<endl;
     	cout<<"---------------------------"<<endl;
    	cout<<" f = "<<f<<" n = "<<n<<" a = "<<a<<endl; 
    	cout<<" n_mc = "<<n_mc<<" n_eq = "<<n_eq<<endl;
   	cout<<" n_p = "<<n_p<<" n_c = "<<n_c<<endl;
   	cout<<" delta_x = "<<delta_x<<" hot = "<<hot<<" n_cool = "<<n_cool<<endl;
   	cout<<" S0 = " <<action0<<" density = "<<density<<endl;
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
        for(int k=0; k<1000; k++)
        {
        	n_insta_sum[k]=0.0;
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
	//c     monte carlo sweeps                                                     c
	//c----------------------------------------------------------------------------c
	
	int n_accepted = 0;
	int n_hit = 0;
	int n_config = 0;
	int n_cor = 0;
	int n_cool_config = 0;
	
	double x_new;
	double s_old,s_new,delta_s;
	int p0;

	path<<fixed<<"tau"<<"\t\t"<<"x[tau]"<<"\t\t"<<"x_cool[tau]"<<"\t\t"<<"config"<<endl;
    			
	
	for(int i=0; i<n_mc; i++)
	{
		n_config += 1;
		if(i==n_eq)
		{
        		n_config = 0;
        		n_cool_config = 0;
        		n_cor =0;
        		for(int k=0; k<100; k++)
        		{
        		 	x_cor_sum[k]=0.0;
        		 	x2_cor_sum[k]=0.0;
        		 	x3_cor_sum[k]=0.0;
        		 }	
        		 for(int k=0; k<1000; k++)
        		{
        			n_insta_sum[k]=0.0;
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
		
		for(int j=0; j<n; j++) x_cool[j]=x[j];
		x_cool[n]=x[n];
		x_cool[n+1]=x[n+1];
		
		//c----------------------------------------------------------------------------c
		//c     cooling and topological charge                                         c
		//c----------------------------------------------------------------------------c
		
		if(i % k_p == 0)
		{
			n_cool_config += 1;
			instanton(x_cool,n,n_i,n_a,x_i,x_a,a);
            		n_insta = n_i+n_a;
            		n_insta_sum[0] += n_insta;
			for (int i_cool=1; i_cool<n_cool;i_cool++)
         		{
         			cooling(x_cool, n, f, a, delta_x);
            			instanton(x_cool,n,n_i,n_a,x_i,x_a,a);
            			n_insta = n_i+n_a;
         			n_insta_sum[i_cool] += n_insta;
         		}
         		
         		//c----------------------------------------------------------------------------c
			//c     cooled correlator                                                      c
			//c----------------------------------------------------------------------------c
			
			for(int c=0; c<n_c; c++)
    			{
         			n_cor += 1;
         			rnd=((double) rand())/((double)RAND_MAX);
        			p0 = (int)((n - n_p)*rnd);
         			for(int p=0; p<n_p; p++)
        			{
					x_cor = x_cool[p0]*x_cool[p0+p];
            				x2_cor = pow(x_cor,2);
            				x3_cor = pow(x_cor,3);
             				x_cor_sum[p] += x_cor;
            				x2_cor_sum[p] += x2_cor;
            				x3_cor_sum[p] += x3_cor;
            			}
         		}
         	}

		//c----------------------------------------------------------------------------c
		//c     write configuration                                                    c
		//c----------------------------------------------------------------------------c
   
         	if (i % k_p2 == 0) 
         	{
         		for (int k=0; k<n; k++)
         		{
         			path<<fixed<<k*a<<'\t'<<x[k]<<'\t'<<x_cool[k]<<'\t'<<n_config<<endl;
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
    	
    	for(int q=0; q<n_cool; q++)
    	{
    		n_insta_av[q]=n_insta_sum[q]/(double)n_cool_config;
    	
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
	//c     instanton density                                                      c
	//c----------------------------------------------------------------------------c
	
	cout<<"------------------------------------------------"<<endl;
	cout<<fixed<<"n_cool"<<"\t"<<"instanton density n(I+A)"<<endl;
	for(int p = 0; p<n_cool; p++)
	{
   		cout<<fixed<<p<<'\t'<< n_insta_av[p]/tau_max<<endl;  
    	}
    	
}
