#!/bin/sh

#THis can be included in the .bashrc to source all needed to run with ofed

# or just copy and paste when needed.

#if [ `uname -m` = "k1om" ]
#then
#    export LD_LIBRARY_PATH=/opt/intel/composerxe/lib/mic:$LB_LIBRARY_PATH
#    source /opt/intel/impi/5.0.1.035/mic/bin/mpivars.sh
#else
    source /opt/intel/composerxe/bin/compilervars.sh intel64
    source /opt/intel/impi/5.0.1.035/bin64/mpivars.sh
    export I_MPI_MIC=1
    export I_MPI_FALLBACK=1
    export I_MPI_DAPL_PROVIDER_LIST=ofa-v2-mcm-1,ofa-v2-scif0
    #export I_MPI_DAPL_PROVIDER_LIST=ofa-v2-scif0
    export I_MPI_FABRICS=dapl

    export IMBHOST=${I_MPI_ROOT}/bin64/IMB-MPI1
    export IMBMIC=${I_MPI_ROOT}/mic/bin/IMB-MPI1
    ulimit -s unlimited
    #mpirun -np 2   -host  ph-ivyknc2  ${IMBHOST} PingPong 
    #mpirun -np 1  -host  ph-ivyknc2  ${IMBHOST} PingPong : -np 1 -host mic0   ${IMBMIC} PingPong
    #mpirun -np 1  -host  ph-ivyknc2  ${IMBHOST} PingPong : -np 1 -host mic1 ${IMBMIC} PingPong
#fi
#mpiicc -Wall -std=c99 -mkl -openmp par_mpi.c -O2 -o mpi.mic -mmic 
#mpiicc -Wall -std=c99 -mkl -openmp par_mpi.c -O2 -o mpi.cpu 

    mpirun -np 1 -host localhost ./mpi.cpu : -np 1 -host mic0 ./mpi.mic : -np 1 -host mic1 ./mpi.mic
