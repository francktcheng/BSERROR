#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <mpi.h>
#include <omp.h>
#include <immintrin.h>
#include "inc/misc.h"
#include "inc/onetimesimu.h"

#define M 3 //Monte Carlo simulation

// constants for MPI communication
const int bossRank = 0;
const int msgReportLenghth = 2; //Mloc, countloc
int msgReportTag = 2; //variable to control the reception of different iteration
const int msgSchedLength = 2; //Nloc, go 
const int msgSchedTag = 1; 
const int msgNameTag = 0;
const int hostNameLen = 128;
typedef char HostNameType[128];


int main(int argc, char *argv[])
{
  const double Prob = 0.95;
  const int Ninit = initialN();
  int countloc, Mloc, seed, len;
  int stop = 0;
  int rptBuf[msgReportLenghth]; //Mloc, countloc
  int shdBuf[msgSchedLength]; //Nloc, stop
  char myName[hostNameLen];
  
  //Initialize MPI
  MPI_Status mpiStatus;
  MPI_Request req;
  int myRank, mpiWorldSize;
  MPI_Request *requestList, requestNull;


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
  
  if(myRank == bossRank){
    int sumM, sumCount;
    int Nglb = Ninit >> 1, Nnew; 
    int threshold = Ninit/1000;
    double t0, t1;
    requestList = (MPI_Request*)malloc((mpiWorldSize-1)*sizeof(MPI_Request));
    

    MPI_Get_processor_name(&myName[0], &len);
    printf("boss:\t%16s has joined.\n",myName);
    HostNameType workerName[mpiWorldSize-1];
    for (int src = 1; src < mpiWorldSize; ++src){
      MPI_Recv(workerName[src-1], hostNameLen, MPI_CHAR, src, msgNameTag, MPI_COMM_WORLD, &mpiStatus);
      printf("worker:\t%16s has joined.\n",workerName[src-1]);
    }
    
    printf("Ninit = %d\n", Ninit);
    t0 = MPI_Wtime();
    while(!stop){
      sumM = 0; sumCount = 0;
      for (int src = 1; src < mpiWorldSize; ++src){
	MPI_Irecv(rptBuf, msgReportLenghth, MPI_INT, src, msgReportTag, MPI_COMM_WORLD, &(requestList[src-1]));  
      }
      
      countloc = 0;
      Mloc = M%(mpiWorldSize-1);
      seed = M-Mloc;
      for (int m = 0; m < Mloc; ++m){
	oneTimeSimu_OMP2(&countloc, seed+m, Nglb);
      }//loop MC

      sumM = Mloc; sumCount = countloc; 
      //printf("init count=%d, init M=%d\n", sumCount, Mloc);
      for (int i = 1; i < mpiWorldSize; ++i){
	int index;
	MPI_Waitany(mpiWorldSize-1, requestList, &index, &mpiStatus);
	sumM += rptBuf[0];
	sumCount += rptBuf[1];
	//printf("From the process %d: count = %d, M = %d\n", index+1, rptBuf[1], rptBuf[0]);
      }
      if ((double)sumCount/sumM > Prob){
	Nnew = Nglb >> 1;
	if(abs(Nglb-Nnew) < threshold){
	  printf("Nglb=%d, Nnew=%d, diff=%d\n", Nglb, Nnew, abs(Nglb-Nnew));
	  stop = 1;
	}
      }
      else{
	Nnew = (Nglb+Ninit)/2;
	if(abs(Nnew-Nglb) < threshold){
	  printf("Nglb=%d, Nnew=%d, diff=%d\n", Nglb, Nnew, abs(Nglb-Nnew));
	  stop = -1;
	}
      }
      Nglb = Nnew;
      shdBuf[0] = Nglb;
      shdBuf[1] = stop;
      for (int dest = 1; dest < mpiWorldSize; ++dest){
	MPI_Isend((void*)&shdBuf[0], msgSchedLength, MPI_INT, dest, msgSchedTag, MPI_COMM_WORLD, &requestNull);
	//MPI_Request_free(&requestNull);
      }
      printf("iteration %d: N = %d, count=%d, M=%d, prob=%.5lf\n", msgReportTag-1, Nglb, sumCount, sumM, (double)sumCount/sumM);
      msgReportTag++;
    }//while
    t1 = MPI_Wtime();
    printf("Time: %.6lfs", t1-t0);
    if(stop==1)
      printf("The N should be at least %d\n", Nglb);
    else if(stop==-1)
      printf("Something must be wrong. We don't find any approriate N.\n");
  }						
  else{//worker   
    int Nloc = Ninit >> 1, Nnew; 
    const int share = M/(mpiWorldSize-1);
    int nslices = share>=10? 10 : 1;
    int flag = 0;
    int m;
    
    MPI_Get_processor_name(&myName[0], &len);
    MPI_Send((void*)&myName,hostNameLen,MPI_CHAR,bossRank,msgNameTag,MPI_COMM_WORLD);
    
    Mloc = share;
    countloc = 0;
    while(1){
      if(Mloc == share)
	seed = (myRank-1)*share;
      else
	seed = (myRank-1)*share + share - Mloc;
      for (m = 0; m < Mloc; ++m){
	// one-time MC simulation
	oneTimeSimu_OMP2(&countloc, seed+m, Nloc);
      }//loop MC
      rptBuf[0] = share;
      rptBuf[1] = countloc;
      printf("iteration %d From %s: count=%d, remaining Mloc=%d\n", msgReportTag-1, myName, countloc, Mloc);
      MPI_Isend(rptBuf, msgReportLenghth, MPI_INT, bossRank, msgReportTag, MPI_COMM_WORLD, &requestNull);
      //MPI_Request_free(&requestNull);
      
      countloc = 0;
      Nloc = Nloc >> 1;
      seed = (myRank-1)*share;
      for (m = 0; m < share/nslices; ++m){
	for (int mm = 0; mm < nslices; ++mm){
	  oneTimeSimu_OMP2(&countloc, seed+m*nslices+mm, Nloc);
	}
	MPI_Iprobe(bossRank, msgSchedTag, MPI_COMM_WORLD, &flag, &mpiStatus);
	if(flag)
	  break;
      }
      MPI_Recv(shdBuf, msgSchedLength, MPI_INT, bossRank, msgSchedTag, MPI_COMM_WORLD, &mpiStatus);
      if(shdBuf[1]==1 || shdBuf[1]==-1)
	break;
      else{
	Nnew = shdBuf[0];
	if(Nnew == Nloc){
	  if(m!=share/nslices)
	    Mloc = share - (m+1)*nslices;
	  else
	    Mloc = 0;
	}
	else
	  Mloc = share;
      }//eles
      msgReportTag++;
    }
  }//else worker 
  
  printf("hostname(rank=%d):\t%16s has finished.\n", myRank, myName);

  MPI_Finalize();
  return 0;
}
