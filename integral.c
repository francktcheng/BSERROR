#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include "inc/misc.h"

#define NN 322
const double PI = 3.14159265358979323846; /* pi */

double NormalIntegral(double b);
double vNormalIntegral(double b);

int main(int argc, char *argv[])
{
  double integral;
  start_timer();
  integral = NormalIntegral(-0.1146);
  printf ("%.20lf\nseq:%.3lf ms\n",integral,stop_timer());

  start_timer();
  integral = vNormalIntegral(-0.1146);
  printf ("%.20lf\nvec:%.3lf ms\n",integral,stop_timer());
  return 0;
}

#ifdef __MIC__
double vNormalIntegral(double b)
{
  __declspec(align(64)) __m512d vec_cf0, vec_cf1, vec_cf2, vec_s, vec_stp, vec_exp; 
  //NN/2-1 has to be the multiple of 8
  //NN = (8*LV+1)*2, LV = 20 -> NN = 322
  //const int NN = 322; 
  const int vecsize = 8; 
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
  //const int NN = 322;
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


double NormalIntegral(double b)
{
  //if(b < -10.0) return  0.0f;
  //if(b > 10.0) return 1.0f;
  //const int NN = 322;
  double a = 0.0f;
  double s, h, sum = 0.0f;
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



