// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_BORNMAYERHUGGINS_H__
#define __PAIR_EVALUATOR_BORNMAYERHUGGINS_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairLJ.h
    \brief Defines the pair evaluator class for BornMayerHuggins potentials
    \details As the prototypical example of a MD pair potential, this also serves as the primary
   documentation and base reference for the implementation of pair evaluators.
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hoomd
    {
namespace md
    {
//! Class for evaluating the LJ pair potential
/*! <b>General Overview</b>

    Put something here

*/
class EvaluatorPairBornMayerHuggins
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar A;
        Scalar sigma;
        Scalar rhoinv;
        Scalar C;
        Scalar D;

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : A(0), sigma(0), rhoinv(0), C(0), D(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            auto rho(v["rho"].cast<Scalar>());

            A = v["A"].cast<Scalar>();
            sigma = v["sigma"].cast<Scalar>();
            rhoinv = Scalar(1.0) / rho;
            C = v["C"].cast<Scalar>();
            D = v["D"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["A"] = A;
            v["sigma"] = sigma;
            v["rho"] = Scalar(1.0) / rhoinv;
            v["C"] = C;
            v["D"] = D;
            return v;
            }
#endif
        }
#ifdef SINGLE_PRECISION
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorPairBornMayerHuggins(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), A(_params.A), sigma(_params.sigma), rhoinv(_params.rhoinv),
          C(_params.C), D(_params.D)
        {
        }

    //! LJ doesn't use diameter
    DEVICE static bool needsDiameter()
        {
        return false;
        }
    //! Accept the optional diameter values
    /*! \param di Diameter of particle i
        \param dj Diameter of particle j
    */
    DEVICE void setDiameter(Scalar di, Scalar dj) { }

    //! LJ doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional diameter values
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r2inv = Scalar(1.0) / rsq;
            Scalar r6inv = r2inv * r2inv * r2inv;
            Scalar r = fast::sqrt(rsq);

            Scalar rexp = fast::exp((sigma - r)*rhoinv);

            force_divr = (A*rhoinv*rexp*r - Scalar(6.0)*C*r6inv - Scalar(8.0)*D*r6inv*r2inv) * r2inv;
            pair_eng = A*rexp - C*r6inv - D*r6inv*r2inv;
            return true;
            }
        else
            return false;
        }

    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("bornmayerhuggins");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar A;      //!< lj1 parameter extracted from the params passed to the constructor
    Scalar sigma;  //!< lj1 parameter extracted from the params passed to the constructor
    Scalar rhoinv; //!< lj1 parameter extracted from the params passed to the constructor
    Scalar C;      //!< lj1 parameter extracted from the params passed to the constructor
    Scalar D;      //!< lj1 parameter extracted from the params passed to the constructor
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_BORNMAYERHUGGINS_H__
