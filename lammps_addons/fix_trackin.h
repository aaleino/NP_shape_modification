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
/* ----------------------------------------------------------------------
   Swift heavy energy deposition module

   Based on the following:

   Electronic stopping power
   Contributing authors: K. Avchaciov and T. Metspalu
   Information: k.avchachov@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(trackin,FixTrackIn)

#else

#ifndef LMP_FIX_TRACKIN_H
#define LMP_FIX_TRACKIN_H

#include "fix.h"

#include <vector>


namespace LAMMPS_NS {

class FixTrackIn : public Fix {
 public:
  FixTrackIn(class LAMMPS *, int, char **);
  ~FixTrackIn();
  int setmask();
  void init();
  void post_force(int);
  void init_list(int, class NeighList *);
  double compute_scalar();

 private:
  void read_table(const char *);
  void grow_table();

  int iregion;      // region index if used, else -1

  std::vector <std::vector <double>> track;
  std::vector<double> boxmid;

  bool was_deposited; // energy deposition done

  bool initialized = false;
  double trackmaxr;
  double trackdeltar;
  double centerx,centery;

  class NeighList *list;
  class Compute *c_ke;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
