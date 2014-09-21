#include <sys/time.h>

//Timer Implementation
static struct timeval tstart, tstop;

void start_timer(){
  gettimeofday(&tstart, NULL);
}

double stop_timer(){
  gettimeofday(&tstop, NULL);
  return (tstop.tv_sec*1000.0 + tstop.tv_usec/1000.0) - (tstart.tv_sec*1000.0 + tstart.tv_usec / 1000.0); //in ms
}

void print_vec(char* name, double *ptr, const int n)
{
  int i;
  for (i = 0; i < n; ++i){
    printf("%s(%d)=%.4lf ",name,i,ptr[i]);
  }
  printf("\n\n");
}
