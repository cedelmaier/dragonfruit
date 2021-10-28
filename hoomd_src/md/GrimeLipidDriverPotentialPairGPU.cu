// Copyright (c) 2009-2021 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

/*! \file GrimeLipidDriverPotentialPairGPU.cu
    \brief Defines the driver functions for computing all types of pair forces on the GPU
*/

#include "AllDriverPotentialPairGPU.cuh"
#include "EvaluatorPairGrimeLipid.h"

hipError_t gpu_compute_grimelipid_forces(const pair_args_t& args,
                                         const EvaluatorPairGrimeLipid::param_type* d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairGrimeLipid>(args, d_params);
    }
