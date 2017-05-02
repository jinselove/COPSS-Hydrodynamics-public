#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh



# 1. *****************************************************************************************
# ========================== FieldSplit Schur Complement Solver ====================
# with schur_user_ksp=true
# mesh          :  50x10x50                                 100x20x100(5+M dofs)
# solver time   :  82.4527 s                                4520.9760 s
# iteration step:  42 KSP Residual norm 1.343672186970e-06  44 KSP Residual norm 2.435167699399e-06
# This works very fast only for schur_user_ksp = true, but needs more iterative steps to converge;
# if false, the iteration is less for convergence, but total time used is more!
# For some special mesh configuration, e.g. 40X16X16 or coarser, it fails to converge! (if true)
#-malloc_dump # outputs all unfreed memory after PetscFinalize() to see the leak.
#-memory_view -log_view

 echo ------------------------------------------------------------------
 echo ------------------------ Use 1 processors ------------------------
 mpirun -n 1 ./example-opt \
 -ksp_type gmres  -ksp_gmres_restart 100   \
 -pc_type fieldsplit                       \
 -pc_fieldsplit_type schur                 \
 -pc_fieldsplit_schur_fact_type lower      \
 -fieldsplit_0_ksp_type cg            \
 -fieldsplit_0_pc_type bjacobi               \
 -fieldsplit_1_pc_type bjacobi             \
 -fieldsplit_1_ksp_type cg              \
 -malloc_dump -malloc_debug -ksp_monitor -fp_trap \
 -memory_view -log_view 2>&1 | tee log_history.txt
###-start_in_debugger noxterm
###-fieldsplit_0_ksp_monitor
###-fieldsplit_0_ksp_max_it 100
#valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./example-opt 


