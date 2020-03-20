

//R_A_N_2_P_A_R_A_M_E_T_E_R_S___________________________________________________

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 3.14159
//#1 DECLARATION CALLING SEQUENCE------------------------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
//#include <conio.h>
#include <curses.h>
#include <stdlib.h>    
#include <stdio.h>

using namespace std;


#define T 1000*3600*10
#define T_window 1000
#define t_sampling 200
#define Sub_sampling_size T/t_sampling
#define sample_size 10
#define N 96
#define Dt 0.001 // one step is equal to 0.1 
#define Pi 3.14159
#define omega_o 0.65//0.8 //this is 1 cycle per 2*Pi rad (sec)
#define dt 0.1
#define c_o 3//3
#define c_min 3.0//3//0.8//4.0
#define c_absolute_min 3
#define c_max 100
#define g 1
#define epsilon 0.2
//myelin plasticity 
#define alpha_c 0.001//0.001
#define c_decay -0.001

#define t_on 0
#define t_off T//3*T/4

#define Damage 0.8
#define t_perturbation_on T/2//T/2+1
#define t_perturbation_off t_perturbation_on+2


float ran2(long *idum);//random number generator - initial cond.
double M(double some_minimal_speed);
double Heaviside(double input);
void four1(double data[], unsigned long nn, int isign);

double theta[N][T_window];
int tau[N][N];
double W_temp[N][N];
double W[N][N];
double W_damage[N][N];
double L[N][N];
double Gamma[N][N];
int tau_o[N][N];
double 	delay[N][N];
double delay_o[N][N];
double l[N][N];
double c[N][N];
double r[Sub_sampling_size];
double mean_theta[Sub_sampling_size];
int sample_i[sample_size];
int sample_j[sample_size];
double sample_speed_1[Sub_sampling_size];
double sample_speed_2[Sub_sampling_size];
double sample_speed_3[Sub_sampling_size];
double sample_speed_4[Sub_sampling_size];
double sample_speed_5[Sub_sampling_size];
double sample_speed_6[Sub_sampling_size];
double sample_speed_7[Sub_sampling_size];
double sample_speed_8[Sub_sampling_size];
double sample_speed_9[Sub_sampling_size];
double sample_speed_10[Sub_sampling_size];
/*
double sample_delays_1[T];
double sample_delays_2[T];
double sample_delays_3[T];
double sample_delays_4[T];
double sample_delays_5[T];
double sample_delays_6[T];
double sample_delays_7[T];
double sample_delays_8[T];
double sample_delays_9[T];
double sample_delays_10[T];
*/
int main()
{
srand((unsigned)time(NULL));//CPU time reader call - for random seeds

		
			cout<<"Attention: Phase shuffle"<<endl;

		
			cout<<"Importing Weights..."<<endl;
		    ifstream pfile("../data/Weights96.txt"); 
		     for (int i = 0 ; i < N ; i++) 
		     { 
		        for (int j = 0 ; j < N ; j++) 
		        { 
		            pfile>>W_temp[i][j];
		            //cout<<W_temp[i][j]<<endl;            
		        } 
		     } 
		     
		     
		      for (int i=0;i<N;i++)
				 {
				 	for (int j=0;j<N;j++)
				 	{
												
						W[i][j] = Heaviside(W_temp[i][j]);
				 	    
				
				 	}
				 }
				 for (int i=0;i<N;i++)
				{
							for (int j=0;j<N;j++)
							{
										long d=rand(); 
										//long seed3= (long) 89*i+454+j+25*i*+50*i*d+4*234;
										long seed3= (long) 89*i+454+j+25*i+50*i*d+4*234;
										if(ran2(&seed3)<Damage)
										{
										        	
										        	W_damage[i][j]=0;
										}
										else
										{
												    W_damage[i][j]=W[i][j];
										}
											
							}
				}
				 
				 
			 ifstream sfile("../data/Lengths96.txt"); 
			 double max_l=0;
		     for (int k = 0 ; k < N ; k++) 
		     { 
		        for (int j = 0 ; j < N ; j++) 
		        { 
		            sfile>>L[k][j];
					 if(L[k][j]>max_l)
					 {
					 	max_l=L[k][j];
					 }
		       
			    } 
		     } 
		     
		     for (int k = 0 ; k < N ; k++) 
		     { 
		        for (int j = 0 ; j < N ; j++) 
		        { 
		            
					 Gamma[k][j]=0.01*L[k][j]/max_l;
			    } 
		     } 
		
				//initialization of network structure + delays
				for (int i=0;i<N;i++)
				{
					for (int j=0;j<N;j++)
					{	
						
				        delay_o[i][j] = (L[i][j])/(c_o+0.0001);
			   			c[i][j] = c_o;	
			   			delay[i][j] = delay_o[i][j];
			   			tau_o[i][j]=(int) floor( delay_o[i][j]);
			   			tau[i][j] = tau_o[i][j]; 

					}
					
				}
				//initial conditions
				for(int t=0;t<T_window;t++)
			 	{
			 		
			 		for (int i=0;i<N;i++)
			   		{
			   			long d=rand(); 
				        long seed1= (long) 4*t+t+i+5*d+65+46+1;
			   			theta[i][t]=2*Pi*ran2(&seed1);
			   			
					}
					
				}
				
				//define sample
				int s=0;
				 for (int i=2;i<N;i++)
				 {
				 	for(int j=2;j<N;j++)
				 	{
				 		if(W_damage[i][j]>0)
				 		{
				 			if(s<sample_size)
				 			{
							 
				 				sample_i[s]=i;
								sample_j[s]=j;	
								s++;
							}
						}
					 }
				 }
			
			

//initialize slidding window
for (int t=0;t<T_window-1;t++)
{
			cout<<t<<endl;
			double mean=0;
			double A=0;
			double B=0;
			for (int i=0;i<N;i++)
			{
					double sum=0;
				
					for (int j=0;j<N;j++)
					{
							sum=sum+g*1/((double)N)*W[i][j]*sin(theta[j][t-tau[i][j]]-theta[i][t]);
					}
					theta[i][t+1] = theta[i][t]+dt*(omega_o+sum);
					double temp = theta[i][t+1]/(2*Pi);
					theta[i][t+1]=2*Pi*(temp-floor(temp));   
					mean = mean + 1/(double(N))*theta[i][t+1];
					A=A+1/((double) N)*cos(theta[i][t+1]);
					B=B+1/((double) N)*sin(theta[i][t+1]);		
			}
			double sr_ratio=0;//or 0.2,0.4,0.6,0.8,1.0;
			for (int i=0;i<N;i++)
			{
				for (int j=0;j<N;j++)
				{
											double phase_sync = sin(-(theta[j][t]-theta[i][t])); 
											double Mij;
											if(phase_sync>0)
											{
												Mij=phase_sync;
											}
											else
											{
												Mij=sr_ratio*phase_sync;
											}
											double grad=-Gamma[i][j]*(c[i][j]-c_min)+epsilon*W[i][j]*Mij; 
										    c[i][j] = c[i][j]+dt*alpha_c*(grad);
										       				    
										    if (c[i][j]>c_max)
										  	{
										       	c[i][j]=c_max;
											}
											if (c[i][j]<c_min)
										    {
										    	c[i][j]=c_min;
											}				    
										    delay[i][j]=(L[i][j])/c[i][j];
										    tau[i][j] = (int) floor(delay[i][j]);
				}
			}


			

}// T_window loop

//run system
int counter=0;
int index=0;
for (int time_loops=0;time_loops<T-T_window;time_loops++)
{
			
			for (int s=0;s<T_window-1;s++)
			{
				for (int i=0;i<N;i++)
				{
						theta[i][s] = theta[i][s+1]	;
				}	
			}	
			if(time_loops>t_perturbation_on && time_loops<t_perturbation_off)
			{
								//apply damage
								for (int i=0;i<N;i++)
								{
									for (int j=0;j<N;j++)
									{
																				        	
										        	W[i][j]=W_damage[i][j];
										
											
									}
								}
								//shuffle
								for (int i=0;i<N;i++)
								{
									
												long d=rand()%100; 
										        theta[i][T_window-1]=(double) d/((double)100)*2*Pi;
									   			
								}
								
				
			}
		
			
			
			double mean=0;
			double A=0;
			double B=0;
			for (int i=0;i<N;i++)
			{
					double sum=0; 
					for (int j=0;j<N;j++)
					{
							sum=sum+g*1/((double)N)*W[i][j]*sin(theta[j][T_window-1-tau[i][j]]-theta[i][T_window-1]);
					}
					theta[i][T_window-1] = theta[i][T_window-1]+dt*(omega_o+sum);
					double temp = theta[i][T_window-1]/(2*Pi);
					theta[i][T_window-1]=2*Pi*(temp-floor(temp)); 
					mean = mean + 1/(double(N))*theta[i][T_window-1];
					A=A+1/((double) N)*cos(theta[i][T_window-1]);
					B=B+1/((double) N)*sin(theta[i][T_window-1]);	
			}
			counter++;
			if (counter==t_sampling)
			{
					mean_theta[index] = mean;
					r[index] = sqrt(A*A+B*B);	
					sample_speed_1[index]=c[sample_i[1]][sample_j[1]];
					sample_speed_2[index]=c[sample_i[2]][sample_j[2]];
					sample_speed_3[index]=c[sample_i[3]][sample_j[3]];
					sample_speed_4[index]=c[sample_i[4]][sample_j[4]];
					sample_speed_5[index]=c[sample_i[5]][sample_j[5]];
					sample_speed_6[index]=c[sample_i[6]][sample_j[6]];
					sample_speed_7[index]=c[sample_i[7]][sample_j[7]];
					sample_speed_8[index]=c[sample_i[8]][sample_j[8]];
					sample_speed_9[index]=c[sample_i[9]][sample_j[9]];
					sample_speed_10[index]=c[sample_i[10]][sample_j[10]];
					index++;
					counter=0;
			}
			cout<<(double)time_loops/((double)T)*100<<"		"<<counter<<"	"<<r[index-1]<<endl;
			
				double sr_ratio=0;//or 0.2,0.4,0.6,0.8,1.0;
			for (int i=0;i<N;i++)
			{
				for (int j=0;j<N;j++)
				{
											double phase_sync =sin(-(theta[j][T_window-1]-theta[i][T_window-1])); 
											double Mij;
											if(phase_sync>0)
											{
												Mij=phase_sync;
											}
											else
											{
												Mij=sr_ratio*phase_sync;
											}
											double grad=-Gamma[i][j]*(c[i][j]-c_min)+epsilon*W[i][j]*Mij; 
										    c[i][j] = c[i][j]+dt*alpha_c*(grad);
										       				    
										    if (c[i][j]>c_max)
										  	{
										       	c[i][j]=c_max;
											}
											if (c[i][j]<c_min)
										    {
										    	c[i][j]=c_min;
											}				    
										    delay[i][j]=(L[i][j])/c[i][j];
										    tau[i][j] = (int) floor(delay[i][j]);
				}
			}
				
			
				
			
			
}//T loop	

			
		ofstream outfile; 
      
        outfile.open("..data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__SubSampledDynamics.txt", ios::out);
         for(int t=0;t<(int)T/t_sampling;t++)
         {
       			
				 	
       				  outfile<<t*Dt*t_sampling<<"	"<<mean_theta[t]<<endl;
         	 	
         } 
        outfile.close();  
       
		
		  outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__SubSampledOrderParamaterPLI.txt", ios::out);
         for(int t=0;t<(int)T/t_sampling;t++)
         {
       			
				 	
       				  outfile<<t*Dt*t_sampling<<"	"<<r[t]<<endl;
         	 	
         } 
        outfile.close();  
 
        outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__SubSampledConductionSpeeds.txt", ios::out);
         for(int t=0;t<(int)(T)/t_sampling;t++)
         {
       			
				 	
       				 	outfile<<t*Dt*t_sampling<<"	"<<sample_speed_1[t]<<"	"<<sample_speed_2[t]<<"	"<<sample_speed_3[t]<<"	"<<sample_speed_4[t]<<"	"<<sample_speed_5[t]<<"	"<<sample_speed_6[t]<<"	"<<sample_speed_7[t]<<"	"<<sample_speed_8[t]<<"	"<<sample_speed_9[t]<<"	"<<sample_speed_10[t]<<endl;
         } 
        outfile.close();  
        
         outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__SubSampledConductionDelays.txt", ios::out);
         for(int t=0;t<(int)T/t_sampling;t++)
         {
       			
				 	
       				   	outfile<<t*Dt*t_sampling<<"	"<<L[sample_i[1]][sample_j[1]]/sample_speed_1[t]<<"	"<<L[sample_i[2]][sample_j[2]]/sample_speed_2[t]<<"	"<<L[sample_i[3]][sample_j[3]]/sample_speed_3[t]<<"	"<<L[sample_i[4]][sample_j[4]]/sample_speed_4[t]<<"	"<<L[sample_i[5]][sample_j[5]]/sample_speed_5[t]<<"	"<<L[sample_i[6]][sample_j[6]]/sample_speed_6[t]<<"	"<<L[sample_i[7]][sample_j[7]]/sample_speed_7[t]<<"	"<<L[sample_i[8]][sample_j[8]]/sample_speed_8[t]<<"	"<<L[sample_i[9]][sample_j[9]]/sample_speed_9[t]<<"	"<<L[sample_i[10]][sample_j[10]]/sample_speed_10[t]<<endl;
         } 
        outfile.close();  
      
        
     
        
          outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__ConductionSpeedMatrix.txt", ios::out);
          for(int i=0;i<N;i++)
         {
       		for (int j=0;j<N;j++)
       		{	 	
       				  outfile<<i<<"	"<<j<<"	"<<c_o<<"	"<<W[i][j]*c[i][j]<<endl;
         	}
         } 
        outfile.close();  
        
        outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__MetabolicDragMatrix.txt", ios::out);
          for(int i=0;i<N;i++)
         {
       		for (int j=0;j<N;j++)
       		{	 	
       				  outfile<<i<<"	"<<j<<"	"<<W[i][j]*Gamma[i][j]<<endl;
         	}
         } 
        outfile.close();             
        
         outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__ConnectivityMatrix.txt", ios::out);
          for(int i=0;i<N;i++)
         {
       		for (int j=0;j<N;j++)
       		{	 	
       				  outfile<<i<<"	"<<j<<"	"<<W[i][j]<<endl;
         	}
         } 
        outfile.close();   
		
		  outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__LengthsMatrix.txt", ios::out);
          for(int i=0;i<N;i++)
         {
       		for (int j=0;j<N;j++)
       		{	 	
       				  outfile<<i<<"	"<<j<<"	"<<W[i][j]*L[i][j]<<endl;
         	}
         } 
        outfile.close();            
        
        
		  outfile.open("../data/sim_results/Fig6__Kuramoto_e0p2_a0p001_damage__DelaysMatrix.txt", ios::out);
          for(int i=0;i<N;i++)
         {
       		for (int j=0;j<N;j++)
       		{	 	
       				  outfile<<i<<"	"<<j<<"	"<<W[i][j]*delay[i][j]<<endl;
         	}
         } 
        outfile.close();
		 
		
			
			
			
		


			
			
			
			
			
			
			
			
      
return 0;    
}

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}


double M(double some_minimal_speed)
{
	double output;
	if(some_minimal_speed>c_min)
	{
		output=1;
	}
	else
	{
		output=0;
	}
	return output;
}


/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

   
double Heaviside(double input)
{
	double output=0;
	if (input>0)
	{
		output=1;
	}
	return output;
}

