#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <omp.h>
#include <immintrin.h>
#include "inc/misc.h"

#define vecsize 8

#define ALIGNED __attribute__((aligned(64)))
#define N 100000
#define Nvec 5200 //tunable, multiple of 8
#define GUIDED_CHUNK 1 //tunable

const double PI = 3.14159265358979323846; /* pi */
const double X0 = 10; //price of risky asset at t0
const double SIGMA = 0.5; //volatility of risky asset
const double K = 12; //strike price of the option
const double T = 1.0; //muaturity time 
const unsigned long long M = 1; //Monte Carlo Simulation
//const double EPSILON = X0*1.0E-2; //threshold value
const double prob = 0.95;

double vNormalIntegral(double b);

int main(int argc, char *argv[])
{
  unsigned long long count = 0;
  double EPSILON = X0*1.0E-2;
  double err = 0.0;
  const double dt = T/N;
  const double rootdt = sqrt((double)T/N);
  int nCal = N/Nvec;
  const int left = N%Nvec;
  VSLStreamStatePtr stream; 
  int errcode = vslNewStream(&stream, VSL_BRNG_MT2203, 0);//seed=0

  start_timer();
  for (unsigned long long m = 0; m < M; ++m){
    // one-time MC simulation
#pragma omp parallel default(none) shared(err, rootdt, dt, nCal)
    {
      double NRV[Nvec] ALIGNED; // normal distribution random vector
      double BM[Nvec] ALIGNED; // brownian motion
      double PX[Nvec+1] ALIGNED; // price 
      double errloc = 0.0;
      double BMlast = 0.0;
      double PXlast = 0.0;
      double tmp, Tj, upbd;
      int lastk = -1;
      int i,j,k;
      VSLStreamStatePtr stream;
      int errcode = vslNewStream(&stream, VSL_BRNG_MT2203, 0);//seed=0
#pragma omp for schedule(guided,GUIDED_CHUNK) nowait
      for (k = 0; k < nCal; ++k){
	//update BM to the questioned position
	tmp = 0.0;
	for (i = 0; i < k-lastk-1; ++i){
	  vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, Nvec, NRV, 0.0f, 1.0f);
#pragma simd reduction(+:tmp) vectorlengthfor(double) assert
	  for (j = 0; j < Nvec; ++j){
	      tmp += NRV[j];
	  }//loop Nvec reduction
	}//loop k-lastk-1
	if(i!=0){
	  BMlast += tmp*rootdt - NRV[Nvec-1]*rootdt;
	  PXlast = X0*exp(-0.5*SIGMA*SIGMA*(k*Nvec)*dt + SIGMA*BMlast);
	  BMlast += NRV[Nvec-1]*rootdt;
	}
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, Nvec, NRV, 0.0f, 1.0f);
	BM[0] = BMlast + rootdt*NRV[0];
	PX[0] = PXlast;
	for (i = 1; i < Nvec; ++i){
	  tmp = BM[0];
	  for (int j = 1; j <= i; ++j){
	    tmp += NRV[j];
	  }
	  BM[i] = tmp*rootdt;
	  PX[i+1] = X0*exp(-0.5*SIGMA*SIGMA*(k*Nvec+i+1)*dt+SIGMA*BM[i]);
	}//loop Nvec
	PX[1] = X0*exp(-0.5*SIGMA*SIGMA*(k*Nvec+1)*dt+SIGMA*BM[0]);
	for (i = 0; i < Nvec; ++i){
	  j = k*Nvec + i;
	  Tj= j*(double)T/N;
	  upbd = (log(PX[i]/K)+0.5*SIGMA*SIGMA*(T-Tj))/(SIGMA*sqrt(T-Tj));
	  errloc += -1/(sqrt(2*PI))*(PX[i+1]-PX[i])*vNormalIntegral(upbd);
	}
	BMlast = BM[Nvec-1];
	PXlast = PX[Nvec];
	lastk = k;
      }//loop nCal
#pragma omp single 
      {
	tmp = 0.0;
	for (i = 0; i < k-lastk-1; ++i){
	  vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, Nvec, NRV, 0.0f,1.0f);
#pragma simd reduction(+:tmp) vectorlengthfor(double) assert
	  for (j = 0; j < Nvec; ++j){
	    tmp += NRV[j];
	  }//loop Nvec reduction
	}//loop k-lastk-1
	if(i!=0){
	  BMlast += tmp*rootdt - NRV[Nvec-1]*rootdt;
	  PXlast = X0*exp(-0.5*SIGMA*SIGMA*(k*Nvec)*dt + SIGMA*BMlast);
	  BMlast += NRV[Nvec-1]*rootdt;
	}
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, left, NRV, 0.0f, 1.0f);
	BM[0] = BMlast + rootdt*NRV[0];
	PX[0] = PXlast;
	for (i = 1; i < left; ++i){
	  tmp = BM[0];
	  for (int j = 1; j <= i; ++j){
	    tmp += NRV[j];
	  }
	  BM[i] = tmp*rootdt;
	  PX[i+1] = X0*exp(-0.5*SIGMA*SIGMA*(nCal*Nvec+i+1)*dt+SIGMA*BM[i]);
	}//loop Nvec
	PX[1] = X0*exp(-0.5*SIGMA*SIGMA*(nCal*Nvec+1)*dt+SIGMA*BM[0]);
	for (i = 0; i < left; ++i){
	  j = nCal*Nvec + i;
	  Tj= j*(double)T/N;
	  upbd = (log(PX[i]/K)+0.5*SIGMA*SIGMA*(T-Tj))/(SIGMA*sqrt(T-Tj));
	  errloc += -1/(sqrt(2*PI))*(PX[i+1]-PX[i])*vNormalIntegral(upbd);
	}
      }//single
#pragma omp atomic
      err += errloc;
    }//parallel
    
  }//MC simulations  
  printf ("time %g ms\n", stop_timer());
  printf("err=%.10lf\n",err);
  printf("count=%llu, M=%llu\n", count, M);
  printf("%.5g\n", (double)count/(double)M);

  return 0;
}

#ifdef __MIC__
double vNormalIntegral(double b)
{
  __declspec(align(64)) __m512d vec_cf0, vec_cf1, vec_cf2, vec_s, vec_stp, vec_exp; 
  
  const int NN = 1000; //has to be the multiple of 8
  //const int vecsize = 8; 
  const int nCal = (NN/2-1)/vecsize;
  //const int left = NN%vecsize;
  double a = 0.0f;
  double s, h, sum = 0.0f;
  h = (b-a)/NN;
  // add in the first few terms 
  sum += exp(-a*a/2.0) + 4.0*exp(-(a+h)*(a+h)/2.0);
  // and the last one
  sum += exp(-b*b/2.0);

  vec_cf0 = _mm512_set1_pd(a);
  vec_cf1 = _mm512_set1_pd(2*h);
  vec_cf2 = _mm512_set1_pd(-0.5);

  vec_s   = _mm512_set_pd(8,7,6,5,4,3,2,1);//vectorize
  vec_s   = _mm512_mul_pd(vec_s, vec_cf1);//(16h,14h,..,2h)
  vec_s   = _mm512_add_pd(vec_cf0, vec_s);//(a+16h,..,a+2h)
  
  vec_stp = _mm512_set1_pd(2*h*vecsize-h);
  vec_cf0 = _mm512_set1_pd(h);
  
  for (int i = 0; i < nCal; ++i){
    vec_exp = _mm512_mul_pd(vec_s, vec_s);
    vec_exp = _mm512_mul_pd(vec_exp, vec_cf2);
    vec_cf1 = _mm512_exp_pd(vec_exp);//vec_cf1->sum
    sum    += 2.0*_mm512_reduce_add_pd(vec_cf1);

    vec_s   = _mm512_add_pd(vec_s, vec_cf0);//s+=h
    vec_exp = _mm512_mul_pd(vec_s, vec_s);
    vec_exp = _mm512_mul_pd(vec_exp, vec_cf2);
    vec_cf1 = _mm512_exp_pd(vec_exp);
    sum    += 4.0*_mm512_reduce_add_pd(vec_cf1);
    
    vec_s   = _mm512_add_pd(vec_s, vec_stp);
  }

  sum = 0.5*sqrt(2*PI) + h*sum/3.0;
  return sum;
}
#else
double vNormalIntegral(double b)
{
  const int NN = 1000;
  double a = 0.0f;
  double s, h, sum = 0.0f;
  h = (b-a)/NN;
  // add in the first few terms 
  sum += exp(-a*a/2.0) + 4.0*exp(-(a+h)*(a+h)/2.0);
  // and the last one
  sum += exp(-b*b/2.0);
  
#pragma vector always 
  for (int i = 1; i < NN/2; ++i){
    s = a + 2*i*h;
    sum += 2.0*exp(-s*s/2.0);
    s += h;
    sum += 4.0*exp(-s*s/2.0);
  }
  
  sum = 0.5*sqrt(2*PI) + h*sum/3.0;
  return sum;
}
#endif

