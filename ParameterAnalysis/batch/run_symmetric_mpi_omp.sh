#!/bin/bash

#=====================================================
# table of contents
#
# variable definition
#
# system idenficiation
#
# Intel tools settings
#
# Code related settings
#
# ====================================================


#=====================================================
# variable definition
#
declare -x system_name
declare -x sockets
declare -x proc_family
declare -x proc_model
declare -x proc_cores
declare -x proc_frequency
declare -x system_memory
declare -x memory_speed
declare -x ht_enabled
declare -x turbo_enabled
declare -x bios_version

declare -x application_name
declare -x application_version
declare -x application_vertical
declare -x threads
declare -x workload
declare -x timestamp
declare -x performance_metric
declare -x performance_value
declare -x higher_is_better

# ====================================================
# system idenficiation
#
declare cpu_list="ivyknc"
declare mic_list=`cd /var/mpss; /bin/ls -1d mic?`
declare socket_cores=`cat /proc/cpuinfo | grep "physical id" | grep "0" | wc -l`
declare sys_freq_path="/sys/devices/system/cpu/cpu0/cpufreq"

declare system_name=`/bin/hostname`
declare proc_cores=`cat /proc/cpuinfo | grep 'processor' | wc -l`
declare physical_cores=`cat /proc/cpuinfo | grep "physical id" | grep "0" | wc -l` 
declare line=`cat /proc/meminfo | grep 'MemTotal' | tr -s " "`

let sockets=${proc_cores}/${socket_cores}
let cores_per_socket=${proc_cores}/${sockets}

declare yesiam=`whoami`
declare theboss="root"

# proc frequency
if [ -f ${sys_freq_path}/cpuinfo_max_freq ] ; then
    proc_frequency=`cat ${sys_freq_path}/cpuinfo_max_freq | cut -c 1-4`
else
   proc_frequency="need root access" 
fi

# memory frequency
if [ ${yesiam} = 'root' ] ; then 
memory_speed=`dmidecode --type memory | grep "Speed" | grep "MHz" | head -1 | cut -d ":" -f 2 | tr -d " "`
else
memory_speed="need root access"
fi

# memory size
mem_size=`echo "${line}" | cut -d " " -f 2`
let mem_size=mem_size+1048575
let mem_size=mem_size/1048576
system_memory="${mem_size} GB"

# ht
if [ ${cores_per_socket} = ${physical_cores} ] ; then
    ht_enabled=n
    max_cpu_nb_threads=$proc_cores
else
    ht_enabled=y
    max_cpu_nb_threads=`expr $proc_cores \* 2`
fi

# turbo
if [ -f ${sys_freq_path}/scaling_available_frequencies ] ; then
    core_freq=`cat ${sys_freq_path}/scaling_available_frequencies`
    top_bin=`echo ${core_freq} | cut -d " " -f 1`
    top_bin_minus_1=`echo ${core_freq} | cut -d " " -f 2`
    let freq_diff=${top_bin}-${top_bin_minus_1}
    if [ ${freq_diff} = "1000" ] ; then
	turbo_enabled=y
    else
	turbo_enabled=n
    fi
else
    turbo_enabled=n
fi

# system
sys_name=`/bin/hostname `
system=`echo ${cpu_list} | grep "$sys_name"`

if [ "$system" = "${cpu_list}" ] ; then
    proc_codename="ivytown-ep"
    proc_family=`cat /proc/cpuinfo | grep "model name" | head -1 | cut -d ":" -f 2 | tr -s " "`
    proc_model=`echo ${proc_family} | cut -d " " -f 4`
if [ ${yesiam} = 'root' ] ; then  
    bios_version=`sudo dmidecode --type bios | grep "version" | cut -d ":" -f 2 | tr -s " "`
else
    bios_version="need root access"
fi

else
    proc_codename="xeon-phi"
    proc_family="many integrated cores"
#    proc_family=`cat /proc/cpuinfo | grep "model name" |  head -1 | cut -d ":" -f 2 | tr -s " "`
if [ ${yesiam} = 'root' ] ; then  
    bios_version=`dmidecode --type bios | grep "Release Date" | cut -d ":" -f 2 | tr -s " "`
else
    bios_version="need root access"
fi
printf " ***********************\n"    
printf "you better run micinfo \n"
printf ""    
fi

printf "# ================================================ \n"
printf "# system name:          ${system_name}\n"
printf "# memory speed:         ${memory_speed} \n"
printf "# memory size:          ${system_memory} \n"
printf "# hyper threading:      ${ht_enabled} \n"
printf "# turbo enabled:        ${turbo_enabled} \n"
printf "# sockets:              ${sockets} \n"
printf "# cores per socket:     ${cores_per_socket} \n"
printf "# total cores:          ${proc_cores} \n"
printf "# total threads:        ${max_cpu_nb_threads} \n"
printf "# processor frequency:  ${proc_frequency} \n"
printf "# processor codename:   ${proc_codename} \n"
printf "# processor family:     ${proc_family} \n"
printf "# processor model:      ${proc_model} \n"
printf "# bios version:         ${bios_version} \n"
printf "# ================================================ \n"

# ====================================================  
# Intel tools settings
#
export PATH_TO_INSTALL="/opt/intel"
export COMPILERDIR="${PATH_TO_INSTALL}/composerxe"
export MPI="${PATH_TO_INSTALL}/impi/4.1.1"

source ${COMPILERDIR}/bin/compilervars.sh intel64
source ${MPI}/intel64/bin/mpivars.sh
#source ${COMPILER}/mkl/bin/mklvars.sh intel64 lp64
source ${PATH_TO_INSTALL}/itac/8.1.3/bin/itacvars.sh 

IMPIDIR=$I_MPI_ROOT
ITACDIR=$VT_ROOT

export MIC_LD_LIB_DIR="$LD_LIBRARY_PATH:${COMPILERDIR}/lib/mic"
export MPIBINDIR="$IMPIDIR/mic/bin"
export MPILIBDIR="$IMPIDIR/mic/lib"
export ITACLIBDIR="$ITACDIR/mic/itac/slib_impi4"
export COMPILERLIB="$COMPILERDIR/lib/mic"
export I_MPI_MIC=1
export I_MPI_MIC_POSTFIX=
#.mic
export I_MPI_MIC_PREFIX=
export I_MPI_ENV_PREFIX_LIST=
#knc:PHI 
export TRACE_OPTION=
#-t

# ==================================================== 
# Code related settings
#
# ====================================================
# ==================================================== 

#------------------
# Xeon Phi settings
export MIC_AFFINITY="granularity=fine,balanced"
#export MIC_PLACEMENT=61Cx4T
export MIC_LIBRARY=throughput
#export MIC_NUM_THREADS=
export max_mic_nb_cores=60
export max_mic_nb_thread=`expr 4 \* $max_mic_nb_cores`
export factor_speed_xphi=1
num_mic=`echo $mic_list | wc -w`
export mic_thread_list='1 2 3 4 5 6 8 10 12 15 16 20 24 30 40 48 60 80 120 240'
#export mic_thread_list='80'

printf "# ================================================ \n"
echo "Xeon Phi settings"
echo $MIC_AFFINITY $MIC_LIBRARY $max_mic_nb_cores $max_mic_nb_thread $factor_speed_xphi
echo $mic_thread_list

#-----------------
# Cpu settings
export HOST_AFFINITY=compact
export HOST_LIBRARY=throughput
export factor_speed_xeon=5
export cpu_thread_list='1 2 3 4 6 8 12 16'
echo "Cpu settings "
echo $HOST_AFFINITY $HOST_LIBRARY  $max_cpu_nb_threads $factor_speed_xeon
echo $cpu_thread_list
export method_list='rtm fwi'
# !!
export cpu_bin_name=../../bin/fwi2d_4_2_v12_avx.exe
export mic_bin_name=../../bin/fwi2d_4_2_v12_mic.exe

echo $mic_bin_name  
echo $cpu_bin_name

export myhost=172.31.1.254

# -------------------------
# Size and load balance!!
how_many_Bytes_per_MIC=`expr 13 \* 1024 \* 1024 \* 1024`

echo "Max MIC memory ${how_many_Bytes_per_MIC} bytes"
n3_mic_base=`expr $how_many_Bytes_per_MIC \/ 4`  # single precision is 4 bytes
n3_mic_base=`expr $n3_mic_base \/ 3`  # divided by number of arrays (prev, next, vel)
#echo "n3_mic_base=$n3_mic_base"
HOSTS="$myhost"
#echo "HOSTS=$HOSTS"
mic_n3_contrib=`expr $num_mic \* $n3_mic_base`  # MIC contribution
cpu_n3_contrib=`expr $sockets \* \( \( $n3_mic_base \* $factor_speed_xeon \) \/ $factor_speed_xphi \)` 
n3_base=`expr $cpu_n3_contrib \+ $mic_n3_contrib`
common_denon=`expr \( $sockets \* $factor_speed_xeon \) \+ \( $num_mic \* $factor_speed_xphi \)`
n3_base=`expr \( $n3_base \/ $common_denon \) \* $common_denon`

#echo "n3_base=$n3_base"
# ca c'est faux : extraite le nombre de mic
num_hosts=`echo $HOSTS | wc -w`
n3=`expr $num_hosts \* $n3_base`
#echo $n3
#echo $mic_list
#echo $mic_thread_list

export  global="-genv I_MPI_DAPL_PROVIDER_LIST ofa-v2-mlx4_0-1u,ofa-v2-scif0,ofa-v2-mcm-1"
export global="$global -genv I_MPI_FALLBACK 0"
export global="$global -genv I_MPI_RDMA_MAX_MSG_SIZE 1048576"
export global="$global -genv I_MPI_DAPL_DIRECT_COPY_THRESHOLD 4096"
export global="$global -genv I_MPI_DAPL_TRANSLATION_CACHE_MAX_ENTRY_NUM 16384"
export global="$global -genv I_MPI_EAGER_THRESHOLD 4096"
export WORKDIR=/home/fye/phil/validation_test/marm_interp

# =====================================================
# real runs
# =====================================================


export i_cpu_mpi=4
export i_mic_mpi=4
export omp_num_threads=6
export mic_omp_num_threads=60
export processor=fanknc
export I_MPI_PIN_DOMAIN=auto
set -x
# single card
export sym_log_name="mod_symmetric_${processor}_${i_cpu_mpi}_${i_mic_mpi}mpi_${omp_num_threads}_${mic_omp_num_threads}omp.log"

totsrc=`expr $i_cpu_mpi \+ $i_mic_mpi`
totsrc=20
sed -i "s/^number of sources           =.*/number of sources           = ${totsrc}/" par_list.in
grep 'number of sources' par_list.in 
sed -i "s/^modelling                          =.*/modelling                          = yes/" par_list.in
grep 'modelling  ' par_list.in 
sed -i "s/^run fwi                            =.*/run fwi                            = no/" par_list.in
grep 'run fwi  ' par_list.in 

    mpiexec.hydra ${TRACE_OPTION} -np ${i_cpu_mpi} -host ${processor} \
	-env KMP_AFFINITY $HOST_AFFINITY -env KMP_LIBRARY $HOST_LIBRARY \
	-wdir ${WORKDIR} -env OMP_NUM_THREADS ${omp_num_threads} ${cpu_bin_name} : \
	-np ${i_mic_mpi} -host mic0 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH ${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}  : \
	-np ${i_mic_mpi} -host mic1 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH ${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}  : \
	-np ${i_mic_mpi} -host mic2 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH ${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name} : \
	-np ${i_mic_mpi} -host mic3 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH ${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}  
 # &> ${sym_log_name}


exit 0

export sym_log_name="fwi_symmetric_${processor}_${i_cpu_mpi}_${i_mic_mpi}mpi_${omp_num_threads}_${mic_omp_num_threads}omp.log"
sed -i "s/^modelling                          =.*/modelling                          = no/" par_list.in
sed -i "s/^run fwi                            =.*/run fwi                            = yes/" par_list.in

    mpiexec.hydra ${TRACE_OPTION} \
	-np ${i_cpu_mpi} -hosts ${processor} \
	-env KMP_AFFINITY $HOST_AFFINITY -env KMP_LIBRARY $HOST_LIBRARY  \
	-wdir ${WORKDIR} -env OMP_NUM_THREADS  ${omp_num_threads} ${cpu_bin_name} : \
	-np ${i_mic_mpi} -hosts mic0 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}   &> ${sym_log_name}


# double card
export sym_log_name="bi_mod_symmetric_${processor}_${i_cpu_mpi}_${i_mic_mpi}mpi_${omp_num_threads}_${mic_omp_num_threads}omp.log"

totsrc=`expr $i_cpu_mpi \+ $i_mic_mpi \+ $i_mic_mpi`
sed -i "s/^number of sources           =.*/number of sources           = ${totsrc}/" par_list.in
grep 'number of sources' par_list.in 

sed -i "s/^modelling                          =.*/modelling                          = yes/" par_list.in
sed -i "s/^run fwi                            =.*/run fwi                            = no/" par_list.in

    mpiexec.hydra ${TRACE_OPTION} \
	-np ${i_cpu_mpi} -hosts ${processor} \
	-env KMP_AFFINITY $HOST_AFFINITY -env KMP_LIBRARY $HOST_LIBRARY  \
	-wdir ${WORKDIR} -env OMP_NUM_THREADS  ${omp_num_threads} ${cpu_bin_name} : \
	-np ${i_mic_mpi} -hosts mic0 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name} : \
	-np ${i_mic_mpi} -hosts mic1 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}   &> ${sym_log_name}


export sym_log_name="bi_fwi_symmetric_${processor}_${i_cpu_mpi}_${i_mic_mpi}mpi_${omp_num_threads}_${mic_omp_num_threads}omp.log"
sed -i "s/^modelling                          =.*/modelling                          = no/" par_list.in
sed -i "s/^run fwi                            =.*/run fwi                            = yes/" par_list.in

    mpiexec.hydra ${TRACE_OPTION} \
	-np ${i_cpu_mpi} -hosts ${processor} \
	-env KMP_AFFINITY $HOST_AFFINITY -env KMP_LIBRARY $HOST_LIBRARY  -env I_MPI_PIN_DOMAIN=auto \
	-wdir ${WORKDIR} -env OMP_NUM_THREADS  ${omp_num_threads} ${cpu_bin_name} : \
	-np ${i_mic_mpi} -hosts mic0 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name} : \
	-np ${i_mic_mpi} -hosts mic1 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}   &> ${sym_log_name}

#!!

export TRACE_OPTION=-t


export sym_log_name="bi_fwi_symmetric_${processor}_${i_cpu_mpi}_${i_mic_mpi}mpi_${omp_num_threads}_${mic_omp_num_threads}omp.log"
sed -i "s/^modelling                          =.*/modelling                          = no/" par_list.in
sed -i "s/^run fwi                            =.*/run fwi                            = yes/" par_list.in

    mpiexec.hydra ${TRACE_OPTION} \
	-np ${i_cpu_mpi} -hosts ${processor} \
	-env KMP_AFFINITY $HOST_AFFINITY -env KMP_LIBRARY $HOST_LIBRARY   \
	-wdir ${WORKDIR} -env OMP_NUM_THREADS  ${omp_num_threads} ${cpu_bin_name} : \
	-np ${i_mic_mpi} -hosts mic0 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name} : \
	-np ${i_mic_mpi} -hosts mic1 -env KMP_AFFINITY $MIC_AFFINITY -env KMP_LIBRARY $MIC_LIBRARY \
	-env LD_LIBRARY_PATH=${MIC_LD_LIB_DIR}  -wdir ${WORKDIR}  -env OMP_NUM_THREADS ${mic_omp_num_threads} ${mic_bin_name}   &> ${sym_log_name}


exit 0


