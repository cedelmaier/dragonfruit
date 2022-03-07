   /* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(grime,PairGrime)

#else

#ifndef LMP_PAIR_GRIME_H
#define LMP_PAIR_GRIME_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGrime : public Pair
{
  public:

    PairGrime(class LAMMPS *);
    virtual ~PairGrime() override;

    void compute(int, int) override;

    void settings(int, char **) override;
    void coeff(int, char **) override;
    double init_one(int, int) override;
    double single(int, int, int, int, double, double, double, double &) override;

  protected:

    double** A_vec;
    double** B_vec;
    double** r0_vec;
    double** rc_vec;

    void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
