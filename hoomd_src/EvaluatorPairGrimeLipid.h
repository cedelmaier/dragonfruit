// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_GRIMELIPID_H__
#define __PAIR_EVALUATOR_GRIMELIPID_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorPairGrimeLipid.h
    \brief Defines the pair evaluator class for GrimeLipid potentials
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
//! Class for evaluating the GrimeLipid pair potential
/*! <b>General Overview</b>

    EvaluatorPairGrimeLipid is a thing. 

*/
class EvaluatorPairGrimeLipid
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar A;
        Scalar B;
        Scalar r0;
        Scalar rc;

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
        param_type() : A(0), B(0), r0(0), rc(0) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            A = v["A"].cast<Scalar>();
            B = v["B"].cast<Scalar>();
            r0 = v["r0"].cast<Scalar>();
            rc = v["rc"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["A"] = A;
            v["B"] = B;
            v["r0"] = r0;
            v["rc"] = rc;
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
    DEVICE EvaluatorPairGrimeLipid(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), A(_params.A), B(_params.B), r0(_params.r0), rc(_params.rc)
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
            Scalar a = Scalar(M_PI)/(Scalar(2.0)*r0);
            Scalar b = Scalar(M_PI)/(rc-r0);

            // Unfortunately have to take the square root, as it appears inside of the functions
            Scalar rinv = fast::rsqrt(rsq);
            Scalar r = Scalar(1.0) / rinv;

            // Two regimes of the potential
            if (r <= r0)
                {
                Scalar C = (A/a) - Scalar(2.0)*(B/b);  // Constant of integration
                force_divr = A * fast::cos(r * a) * rinv;
                pair_eng = -(A/a) * fast::sin(r * a) + C;
                }
            else if ((r > r0) && (r <= rc))
                {
                Scalar C = -(B/b);
                force_divr = B * fast::cos(Scalar(M_PI)/Scalar(2.0) + (r - r0)*b) * rinv;
                pair_eng = -(B/b) * fast::sin(Scalar(M_PI)/Scalar(2.0) + (r - r0)*b) + C;
                }
            else
                {
                force_divr = 0.0;
                pair_eng = 0.0;
                }
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
        return std::string("grimelipid");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar A;      //!< A parameter extracted from the params passed to the constructor
    Scalar B;      //!< B parameter extracted from the params passed to the constructor
    Scalar r0;     //!< r0 parameter extracted from the params passed to the constructor
    Scalar rc;     //!< rc parameter extracted from the params passed to the constructor
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_LJ_H__
