#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

	//c----------------------------------------------------------------------------c
    	//c     lattice calculation in quantum mechanics.                              c
    	//c----------------------------------------------------------------------------c
    	//c     calculate partition function from adiabatic switching.                 c
    	//c----------------------------------------------------------------------------c
    	//c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
    	//c----------------------------------------------------------------------------c
    	
   
// -----------------------------------------------------------------	
int main()
{
	srand(time(NULL)); //set seed
	double rnd;
	
	double *x=new double[1000];
	double s_alpha_sum=0.0;
	double s_alpha_av=0.0;
	double *s_alpha=new double [1000];
	
	double f,a;
	int n, n_eq, n_mc, k_p, n_alpha;
	bool hot;
	double delta_x, w;
	
	ofstream action;
	
	action.open("action_alpha.dat");
	
	
	cout<<"separation of wells f (f=1.4): ", cin>>f, cout<<endl;
	cout<<"lattice sites n<1000 (n=400): ", cin>>n, cout<<endl;
	cout<<"lattice spacing a (a=0.05): ", cin>>a, cout<<endl;
	cout<<"cold/hot start (0,1): ", cin>>hot, cout<<endl;
	cout<<"equilibration sweeps before measurement (n_eq=100): ", cin>>n_eq, cout<<endl;
	cout<<"monte carlo sweeps (10^4): ", cin>>n_mc, cout<<endl;
	cout<<"update x (delta_x=0.5): ", cin>>delta_x, cout<<endl;
	cout<<"write every k-th configuration: (k_p=1000): ", cin>>k_p, cout<<endl;
	cout<<"oscillator constant of reference potential w (4f=5.6): ", cin>>w, cout<<endl;
	cout<<"number of steps in switching process (n_alpha=20): ", cin>>n_alpha, cout<<endl;
	
	double beta=n*a;
	double d_alpha=1.0/(double)n_alpha;
	
	cout<<"Monte Carlo adiabatic switching"<<endl;
     	cout<<"---------------------------"<<endl;
    	cout<<" f = "<<f<<" n = "<<n<<" a = "<<a<<endl; 
    	cout<<" nmc = "<<n_mc<<" neq = "<<n_eq<<endl;
   	cout<<" w = "<<w<<" n_alpha = "<<n_alpha<<endl;
   	cout<<" delta_x = "<<delta_x<<" hot = "<<hot<<endl;
   	cout<<"---------------------------"<<endl;
   	
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
		s0=(1.0/4.0*pow((x[i+1]-x[i])/a,2)+1.0/4.0*pow(w*x[i],2))*a;
    		s += s0;
    	}
    	
 	//c----------------------------------------------------------------------------c
    	//c     loop over coupling constant                                            c
    	//c----------------------------------------------------------------------------c
    	
    	double alpha;
    	int n_accepted, n_hit, n_config;
    	
    	double x_new;
	double s_old,s_new,delta_s;
	double t,v0,v;
	double delta_s_2;
	
	double f0 = 1.0/beta*log(2.0*sinh(w/2.0*beta)); //free energy of HO
	double f_i=f0;
	double delta_f;
	
	action<<fixed<<"config"<<"\t\t"<<"alpha"<<"\t\t"<<"acceptance rate"<<"\t\t"<<"S"<<endl;
	cout<<fixed<<"alpha"<<"\t\t"<<"T"<<"\t\t"<<"F(alpha)"<<"\t\t"<<"dF"<<"\t\t"<<"F0"<<endl;
    	
    	for(int i_alpha = 0; i_alpha<2*n_alpha; i_alpha++) 
    	{
    		if (i_alpha  <= n_alpha) //Hysteresis effect
        	{
            		alpha = (i_alpha)*d_alpha;
        	} else
        	{
            		alpha = 2. - (i_alpha)*d_alpha;
        	}
        	n_accepted = 0;
		n_hit = 0;
		n_config = 0;
		
		//c----------------------------------------------------------------------------c
		//c     monte carlo sweeps                                                     c
		//c----------------------------------------------------------------------------c
		
		for(int i=0; i<n_mc; i++)
		{
			n_config += 1;
			if(i==n_eq)
			{
        			n_config = 0;
        			s_alpha_sum=0,0;
        		}
        		
        		//c----------------------------------------------------------------------------c
    			//c     one sweep thorough configuration                                       c
    			//c----------------------------------------------------------------------------c
    		
    			for(int j = 1; j<n; j++)
    			{
    				rnd=((double) rand())/((double)RAND_MAX);
         			n_hit +=  1;
         			x_new = x[j] + delta_x*(2*rnd-1);  //Uniform(x[j]-deltax,x[j]+deltax)
         			t = 1.0/4.0*(pow((x[j]-x[j-1])/a,2)+pow((x[j+1]-x[j])/a,2));
         			v0 = 1.0/4.0*pow(w,2)*pow(x[j],2);
         			v = v0 + alpha*(pow((pow(x[j],2) - pow(f,2)),2)-v0);
         			s_old=a*(t+v);
         			
         			t = 1.0/4.0*(pow((x_new-x[j-1])/a,2)+pow((x[j+1]-x_new)/a,2));
         			v0 = 1.0/4.0*pow(w,2)*pow(x_new,2);
         			v = v0 + alpha*(pow((pow(x_new,2) - pow(f,2)),2)-v0);
         			s_new=a*(t+v);
         			
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
    			delta_s_2=0;
    			for(int j=0; j<n; j++)
			{
				t = 1.0/4.0*pow((x[j+1]-x[j])/a,2);
				v0 = 1.0/4.0*pow(w,2)*pow(x[j],2);
         			v = v0 + alpha*(pow((pow(x[j],2) - pow(f,2)),2)-v0);
         			s0=a*(t+v);
    				s += s0;
    				delta_s_2 += a*(pow((pow(x[j],2) - pow(f,2)),2)-v0);
    			}
    			if ((i % k_p)  ==  0) 
    			{
    				action<<fixed<<n_config<<'\t'<<alpha<<'\t'<<double(n_accepted)/double(n_hit)<<'\t'<<s<<endl;
        		}
        		
        		s_alpha_sum += delta_s_2/beta;
        		
        		//c----------------------------------------------------------------------------c
    			//c     next configuration                                                     c
    			//c----------------------------------------------------------------------------c
  		}
  		//c----------------------------------------------------------------------------c
        	//c     averages                                                               c
        	//c----------------------------------------------------------------------------c
        	s_alpha_av = s_alpha_sum/(double)n_config;
        	s_alpha[i_alpha]=s_alpha_av;
        	
        	if ((i_alpha % 2*n_alpha)  ==  0) //hysteresis
        	{
            		delta_f = d_alpha/4.0*s_alpha_av;
        	} else
        	{
            		delta_f = d_alpha/2.0*s_alpha_av;
        	}
       	f_i += delta_f;
       	
       	cout<<fixed<<alpha<<"\t"<<1.0/beta<<"\t"<<f_i<<"\t"<<delta_f<<"\t"<<f0<<endl;
       	    	
        	//c----------------------------------------------------------------------------c
        	//c     end of loop over coupling constants                                    c
        	//c----------------------------------------------------------------------------c 
        }
        //c----------------------------------------------------------------------------c
    	//c     final estimate of integral over coupling                               c
    	//c----------------------------------------------------------------------------c
    	double f_up = 0.0;
    	double f_down = 0.0;
    	
    	
    	for(int i=0; i<n_alpha; i++) //hysteresis
    	{
        	if (i% n_alpha ==  0)
        	{
        		f_up += d_alpha/4.0*s_alpha[i];
        		f_down +=  d_alpha/4.0*s_alpha[i+n_alpha];
        	} else
        	{
        		f_up += d_alpha/2.0*s_alpha[i];
        		f_down +=  d_alpha/2.0*s_alpha[i+n_alpha];
        	}
        }
        delta_f = f_up + f_down;
        f_i = f0 + delta_f;
        cout<<"------------------------------------------FINAL-----------------"<<endl;
        cout<<fixed<<"T"<<"\t\t"<<"F"<<"\t\t"<<"delta_F"<<"\t\t"<<"F0"<<endl;
	cout<<fixed<<1.0/beta<<"\t"<<f_i<<"\t"<<delta_f<<"\t"<<f0<<endl;
   	
}
