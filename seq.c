#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include "inc/misc.h"

#define real double
#define ALIGNED __attribute__((aligned(64)))
#define N 37
#define NN 2962 //integral interval

real NRV[N] ALIGNED; //normal distribution random vector
real BM[N] ALIGNED; //brownian motion
real PX[N+1] ALIGNED; //price

const real PI = 3.14159265358979323846; /* pi */
const real X0 = 12; //price of risky asset at t0
const real SIGMA = 0.5; //volatility of risky asset
const real K = 10; //strike price of the option
const real T = 1.0; //muaturity time 
const unsigned long long M = 2; //Monte Carlo Simulation
const real relEPSILON = 1.0E-1; // threshold value
const real prob = 0.95;

real NormalIntegral(real);

int main(int argc, char *argv[])
{
  unsigned long long count = 0;
  const real EPSILON = X0*relEPSILON;
  real err;
  real upbd1, upbd2;
  const real dt = T/N;
 
  //////read from terminal
  //const real X0 = atof(argv[1]);
  //const real K = atof(argv[2]);
  //const real EPSILON=X0*1.0e-2;

  //For each thread, initialize a random number stream
  VSLStreamStatePtr stream; //stream for random numbers
  //int seed = omp_get_thread_num();
  int seed = 2;
  int errcode = vslNewStream(&stream, VSL_BRNG_MT2203, seed);

  //start_timer();
  for (unsigned long long i = 0; i < M; ++i){
    err = 0.0;
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, N, NRV, 0.0f, 1.0f);
    BM[0] = sqrt((real)T/N)*NRV[0];
    for (int j = 1; j < N; ++j){
      BM[j] = BM[j-1] + sqrt((real)T/N)*NRV[j];
    }

    PX[0] = X0;
    for (int j = 1; j < N+1; ++j){
      PX[j] = X0*exp(-0.5*SIGMA*SIGMA*j*dt+SIGMA*BM[j-1]);
    }

    for (int j = 0; j < N; ++j){
      real Tj = j*(real)T/N;
      real upbd = (log(PX[j]/K)+0.5*SIGMA*SIGMA*(T-Tj))/(SIGMA*sqrt(T-Tj));
      err -= 1/(sqrt(2*PI))*(PX[j+1]-PX[j])*NormalIntegral(upbd); 
      //printf("de=%.10lf\n",1/(sqrt(2*PI))*(PX[j+1]-PX[j])*NormalIntegral(upbd));
    }
    if (PX[N] > K)
      err += PX[N] - K;

    upbd1 = (log(X0/K) + 0.5*SIGMA*SIGMA*T)/(SIGMA*sqrt(T));
    upbd2 = (log(X0/K) - 0.5*SIGMA*SIGMA*T)/(SIGMA*sqrt(T));
    err += K/(sqrt(2*PI))*NormalIntegral(upbd2) - X0/(sqrt(2*PI))*NormalIntegral(upbd1);
    err = fabs(err);
    if(err < EPSILON)
      count++;
    printf("i=%llu: err=%.20lf\n", i, err);
  }
  printf("last time err=%.20lf\n",err);
  //printf ("time %g ms\n",stop_timer());
  printf("count=%llu, M=%llu\n",count, M);
  printf("%.5g\n", (real)count/(real)M);
  
  vslDeleteStream(&stream);
  return 0;
}

real NormalIntegral(real b)
{
  //if(b < -10.0) return  0.0f;
  //if(b > 10.0) return 1.0f;
  //const int NN = 322;//correspond to vNormalIntegral
  real a = 0.0f;
  real s, h, sum = 0.0f;
  h = (b-a)/NN;
  // add in the first few terms 
  sum += exp(-a*a/2.0) + 4.0*exp(-(a+h)*(a+h)/2.0);
  // and the last one
  sum += exp(-b*b/2.0);

  for (int i = 1; i < NN/2; ++i){
    s = a + 2*i*h;
    sum += 2.0*exp(-s*s/2.0);
    s += h;
    sum += 4.0*exp(-s*s/2.0);
  }
  
  sum = 0.5*sqrt(2*PI) + h*sum/3.0;
  return sum;
}


