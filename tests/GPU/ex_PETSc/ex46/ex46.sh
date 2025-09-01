#!/bin/bash
# Create a PETSc matrix from test case and run ex10 on GPU in parallel
echo "Usage: $0 [-n NP] [petsc_options]"

ex=ex46
ROOT=`pwd`
jdd=`basename $ROOT`"_BENCH"
exec=`pwd`/$ex
if [ ! -f $exec ]
then
   [ "$TRUST_CC_BASE_EXTP" != "" ] && export MPICH_CC=$TRUST_cc_BASE_EXTP && export OMPI_CC=$TRUST_cc_BASE_EXTP
   make $ex 1>/dev/null 2>&1 || exit -1
fi
NP=2 && [ "$1" = -n ] && shift && NP=$1 && shift
petsc_options="-da_grid_x 1000 -da_grid_y 1000 -ksp_monitor -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_strong_threshold 0.7 -dm_mat_type mpiaijcusparse -dm_vec_type mpicuda"
# gamg OK amgx OK !
#petsc_options="-da_grid_x 1000 -da_grid_y 1000 -ksp_monitor -ksp_type cg -pc_type amgx -dm_mat_type mpiaijcusparse -dm_vec_type mpicuda"
petsc_options=$petsc_options" "$*
cmd="-ksp_view -log_view $petsc_options"
echo "mpirun -np $NP ./$ex $cmd" | tee $ex.log
touch dummy.data 
trust -nsys dummy $NP $cmd 2>&1 | tee -a $ex.log


