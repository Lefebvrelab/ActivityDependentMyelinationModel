

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

#define Q 21
#define TRIALS 10
#define K 10
#define T 1024*12
#define t_transient 512
#define N 96
#define delay_graining 100
#define speed_graining 100
#define phase_graining 100
#define max_phase Pi
#define max_speed 100
#define max_delay 250				
#define Dt 0.001 // one step is equal to 0.1 
#define Pi 3.14159
#define omega_o 0.8//0.8 //this is 1 cycle per 2*Pi rad (sec)
#define dt 0.1
#define c_o 3
#define c_min 3
#define c_max 100

//myelin plasticity 
#define alpha_c 1//0.001
#define c_decay -0.001
//#define epsilon 30

#define t_on 0
#define t_off T

#define t_perturbation_on T//T/2//T/2+1
#define t_perturbation_off t_perturbation_on+2



float ran2(long *idum);//random number generator - initial cond.
double M(double some_minimal_speed);
double Heaviside(double input);
void four1(double data[], unsigned long nn, int isign);

double theta[N][T];
int tau[N][N];
double W_temp[N][N];
double W[N][N];
double W_o[N][N];
double L[N][N];
int tau_o[N][N];
int tau_half_time[N][N];
double 	delay[N][N];
double delay_o[N][N];
double l[N][N];
double c[N][N];
double Gamma[N][N];

double Phase_differences[N][N];
double Phase_differences_half_time[N][N];
double mean_phase_diff[N][N];
double mean_delay_in_time[N][N];
//double delays_in_time[N][N][T];
double var_tau_in_time[N][N];
double var_phase_diff[N][N];
double  Delay_distribution[T][delay_graining];
double Speed_distribution[T][speed_graining];
double Phase_distribution[T][phase_graining];
double Phase_distribution_theo[phase_graining];
double mean_theta[T];
double r[T];
double mean_delta[T]; 
double mean_delay[T];
double mean_speed[T];
double mean_phase_difference[T];
double var_phase_difference[T];
double var_delay[T];
double var_speed[T];
double single_connection_delays[K][T];
double single_phase_difference[K][T];
double single_connection_speeds[K][T];
double single_phases[K][T];
double Plasticity[T];
double PSD_1[(int)T];//Power Spectral Density 
double PSD_2[(int)T];//Power Spectral Density 
double PSD_3[(int)T];//Power Spectral Density 
double PSD_4[(int)T];
double Freq[T/2];//Array of output frequencies for PSD display

double INSULT[Q];
double OMEGA[Q];
double G[Q];
double EPSILON[Q];
double MEAN_DELAY[Q][Q];
double SIGMA[Q][Q];
double MEAN_R[Q][Q];
double PEAK_POWER[Q][Q];
double PEAK_FREQUENCY[Q][Q];

int sample_i[K];
int sample_j[K];

int main()
{
srand((unsigned)time(NULL));//CPU time reader call - for random seeds

		
			cout<<"Attention: Phase shuffle"<<endl;

		
			cout<<"Importing Weights..."<<endl;
		    ifstream pfile("../data/connectivity/Weights96.txt"); 
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
												
						W_o[i][j] = Heaviside(W_temp[i][j]);
				 	  
				
				 	}
				 }
				 
			 ifstream sfile("../data/connectivity/Lengths96.txt"); 
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
		
		
     
for (int q1=0;q1<Q;q1++)
{
		for (int q2=0;q2<Q;q2++)
		{
	 			INSULT[q1]=0.0;//+1.0*q/((double)Q-1);
	 			EPSILON[q1]=0.2;//0.08;//40;//100;//60;//50;//50;//0+100*q/((double)Q);;//;30
	 			
				G[q1]=0.0+3.0*q1/((double)Q-1);
	 			OMEGA[q2]=0.1+2.9*q2/((double)Q-1);;//omega_o;
       //initial conditions   
       
       
       	for (int trials =0;trials<TRIALS;trials++)
       	{
		 	//reset IC
		 	
		 	for (int t=0;t<T;t++)
		 	{
		 		for (int i=0;i<N;i++)
				{
						theta[i][t] =0;
				}		 		
				mean_theta[t]=0;
		 		
			 }
		 	
	
		 
		 
			 	for(int t=0;t<T;t++)
			 	{
			 		for (int i=0;i<N;i++)
			   		{
			   			long d=rand(); 
				        long seed1= (long) 4*t+t+i+5*d*q1+q1+65*q2+trials*46+1;
			   			theta[i][t]=2*Pi*ran2(&seed1);
			   			
					}
					if (t>t_on&&t<t_off)
					{
						Plasticity[t] = 1;
					}
					else
					{
						Plasticity[t]=0;
					}
					
					 
			   
				}
				//delays
				
		
				 
				
				for (int i=0;i<N;i++)
				{
					for (int j=0;j<N;j++)
					{	
						W[i][j]=W_o[i][j];
				        delay_o[i][j] = (L[i][j])/(c_o+0.0001);
			   			c[i][j] = c_o;	
			   			delay[i][j] = delay_o[i][j];
			   			tau_o[i][j]=(int) floor( delay_o[i][j]);
			   			tau[i][j] = tau_o[i][j]; 

					}
					
				}
				
					 	for(int t=t_transient;t<T;t++)
					    {
					  		for (int i=0;i<N;i++)
					 		{
					
					       		double sum=0;
					       		for (int j=0;j<N;j++)
					       		{
					       				sum=sum+G[q1]/((double)N)*W[i][j]*sin(theta[j][t-tau[i][j]]-theta[i][t]);
					       				mean_phase_diff[i][j] = mean_phase_diff[i][j]+1/((double)T)*(theta[j][t]-theta[i][t]);
					       				mean_delay_in_time[i][j] = mean_delay_in_time[i][j]+1/((double)T)*tau[i][j]; 
								}
					       	
							   
							   	theta[i][t+1]=theta[i][t]+dt*(OMEGA[q2]+sum);
								double temp = theta[i][t+1]/(2*Pi);
								theta[i][t+1]=2*Pi*(temp-floor(temp));   
								mean_theta[t] = mean_theta[t]+1/((double)N)*theta[i][t];
								
					
					       }
					       double A=0;
					       double B=0;
					       for (int i=0;i<N;i++)
							{
									A=A+1/((double) N)*cos(theta[i][t]);
									B=B+1/((double) N)*sin(theta[i][t]);
							}
					        r[t] = sqrt(A*A+B*B);					       	
						
					
							if(t>t_perturbation_on && t<t_perturbation_off)
							{
								
								for (int i=0;i<N;i++)
								{
									for (int j=0;j<N;j++)
									{
										long d=rand(); 
										long seed3= (long) 89*i+454*q1+j+t+q2*87+25*i*trials+50*i*d+4*q2*trials*234;
										 if(ran2(&seed3)<INSULT[q1])
										        {
										        	
										        	W[i][j]=0;
												}
											
									}
								}
								
											 		for (int i=0;i<N;i++)
											   		{
											   			long d=rand(); 
												        long seed1= (long) 4*t+t+i+5*d+6565*q2+98*q1+trials*46+1;
											   			theta[i][t+1]=2*Pi*ran2(&seed1);
											   			//cout<<	theta[i][t]<<endl;
													}
									      
										       
								    	
									
							}
					
						
			double sr_ratio=1;//or 0.2,0.4,0.6,0.8,1.0;
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
												double grad=-Gamma[i][j]*(c[i][j]-c_min)+Plasticity[t]*EPSILON[q1]*W[i][j]*Mij; 
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
			}// tloop
								
							double mean_r=0;
							for (int t=3*T/4;t<T;t++)
							{
								mean_r=mean_r+1/((double) T/4)*(r[t]);
							}
			
			
								MEAN_R[q1][q2]= MEAN_R[q1][q2]+1/((double)TRIALS)*mean_r;
			
		/*
	 
				for (int k=0;k<T/2;k++)
					 {
						PSD_1[k] = mean_theta[k+T/2];	
						PSD_2[k] = 0;
						PSD_3[k] = 0; 
						PSD_4[k] = 0;		 
					 }
					 					
					 unsigned long nn=T/4;
					 four1(PSD_1-1, nn,1);
					 

					 for (int k=0;k<T/2;k++)
					 {
						PSD_1[k] = 1/((double)T/2*T/2)*(fabs(PSD_1[k])*fabs(PSD_1[k])+fabs(PSD_1[(int)T/2-k])*fabs(PSD_1[(int)T/2-k]));
					 }
					 
					four1(PSD_2-1, nn,1);
					 

					 for (int k=0;k<T/2;k++)
					 {
						PSD_2[k] = 1/((double)T/2*T/2)*(fabs(PSD_2[k])*fabs(PSD_2[k])+fabs(PSD_2[(int)T/2-k])*fabs(PSD_2[(int)T/2-k]));
					 }
					 
					 four1(PSD_3-1, nn,1);
					 

					 for (int k=0;k<T/2;k++)
					 {
						PSD_3[k] = 1/((double)T/2*T/2)*(fabs(PSD_3[k])*fabs(PSD_3[k])+fabs(PSD_3[(int)T/2-k])*fabs(PSD_3[(int)T/2-k]));
					 }
					 
					  four1(PSD_4-1, nn,1);
					 

					 for (int k=0;k<T/2;k++)
					 {
						PSD_4[k] = 1/((double)T/2*T/2)*(fabs(PSD_4[k])*fabs(PSD_4[k])+fabs(PSD_4[(int)T/2-k])*fabs(PSD_4[(int)T/2-k]));
					 
					 
					 }
					
					//compute peak power
					double max_f=0;
					double max_power=0;
					for(int k =2*1*(T/4)*Dt; k<2*100*(T/4)*Dt;k++)
					{
					if (	PSD_1[k] >max_power)
					{
									max_power = PSD_1[k];
									max_f=k*1/((double) 2*(T/2)*Dt);
					}	
					}
								
										 
					PEAK_POWER[q1][q2]=PEAK_POWER[q1][q2]+1/((double)TRIALS)*max_power;
					
					PEAK_FREQUENCY[q1][q2] =PEAK_FREQUENCY[q1][q2]+1/((double)TRIALS)*max_f; 
		*/							
			}//trials
				cout<<G[q1]<<"	"<<OMEGA[q2]<<"	"<<MEAN_R[q1][q2]<<endl;

	}// q1 loop
}//q2-loop
    
    	for(int k=0;k<(int) T/4;k++)
					
					{
								Freq[k] = k*1/((double) 2* (T/2)*Dt);
							
					}

       
          ofstream outfile; 
          
          
        	outfile.open("../data/sim_results/Fig7__Epsilon0p2_sr1p0__OPvsGandOMEGAvsOP0.2.txt", ios::out);
		
	
				for (int q1=0;q1<Q;q1++)
				{
					for (int q2=0;q2<Q;q2++)
					{
					 
					 	
			         					outfile<<G[q1]<<"	"<<OMEGA[q2]<<"	"<<MEAN_R[q1][q2]<<endl;
			         		 
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

