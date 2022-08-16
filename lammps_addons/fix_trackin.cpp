/* ----------------------------------------------------------------------
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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_trackin.h"
#include "mpi.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "fix.h"
#include "compute.h"
#include "modify.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "stringlib.h"


using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024
#define MINNEIGHB 5

/* ---------------------------------------------------------------------- */

FixTrackIn::FixTrackIn(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;  // Has compute_scalar
  nevery = 1;       // Run fix every step

  // args: 0 = fix ID, 1 = group ID,  3 = track.in file
  //       4 = centerx 5 = centery   
  // optional rest: "region" <region name>

  if (narg < 6)
    error->all(FLERR, "Illegal fix trackarray command: too few arguments");

#if 0 

  int iarg = 4;
  iregion = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iregion >= 0)
         error->all(FLERR, "Illegal fix elstop command: region given twice");
      if (iarg+2 > narg)
        error->all(FLERR, "Illegal fix elstop command: region name missing");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion < 0)
        error->all(FLERR, "region ID for fix elstop does not exist");
      iarg += 2;
    }
    else error->all(FLERR, "Illegal fix elstop command: unknown argument");
  }

  // need an occasional full neighbor list
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 0;
  neighbor->requests[irequest]->occasional = 0;
  // Read the input file for energy ranges and stopping powers.
  // First proc 0 reads the file, then bcast to others.
  const int ncol = atom->ntypes + 1;
  if (comm->me == 0) {
    maxlines = 300;
    memory->create(elstop_ranges, ncol, maxlines, "elstop:tabs");
    read_table(arg[4]);
  }

  MPI_Bcast(&maxlines, 1 , MPI_INT, 0, world);
  MPI_Bcast(&table_entries, 1 , MPI_INT, 0, world);

  if (comm->me != 0)
    memory->create(elstop_ranges, ncol, maxlines, "elstop:tabs");

  MPI_Bcast(&elstop_ranges[0][0], ncol*maxlines, MPI_DOUBLE, 0, world);
#endif
   // read track.in
   track = split_array_2<double>(read_textfile(arg[3]));	
   boxmid.reserve(3);
   for (int i = 0; i < 3; i++) boxmid[i] = (domain->boxhi[i] + domain->boxlo[i])/2; 
   was_deposited = false;
   trackmaxr = track[track.size()-1][0];
   centerx = atof(arg[4]);
   centery = atof(arg[5]);
//   cout << "centerx " << centerx << " centery " << centery << "\n";

#if 0 
  

  cout << "trackmaxr " << trackmaxr << " trackdeltar" << trackdeltar << "\n";
   // print the data
   for(int i = 0; i < array_timesteps.size(); i++) {
 
      cout << "timestep [" <<  i  << "] = " << std::scientific << array_timesteps[i] << "\n";
      for(int j = 0; j < trackarray[i].size(); j++) {
        printf("%j %18.10e  %18.10e \n", j, trackarray[i][j][0] , trackarray[i][j][1] );
//        cout << "data [" <<  j  << "] = " << std::scientific << trackarray[i][j][0] << " " << trackarray[i][j][1] << "\n";
      }
   }
#endif

}

/* ---------------------------------------------------------------------- */

FixTrackIn::~FixTrackIn()
{
 // memory->destroy(elstop_ranges);
  // delete id_ke_atom;
}

/* ---------------------------------------------------------------------- */

int FixTrackIn::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}



/* ---------------------------------------------------------------------- */

void FixTrackIn::init()
{
   if(initialized) return; 
   initialized=true;
   
}
/* ---------------------------------------------------------------------- */

void FixTrackIn::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double current_time = 0;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> myranf(0, 1);


void FixTrackIn::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double dt = update->dt;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double i_time = 0; // which 'frame' in array_timesteps to use

  if(was_deposited) return; // do this only once
  was_deposited = true; 

  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit)) continue;

    // find the radius of the particle, shifted
    double r_particle = sqrt((x[i][0]-centerx-boxmid[0])*(x[i][0]-centerx-boxmid[0])+(x[i][1]-centery-boxmid[1])*(x[i][1]-centery-boxmid[1]));

    // do not try to deposit beyond maximum radius
    if(r_particle >= trackmaxr) continue;

    // printf("%i: %e\n", i, r_particle);

    
    int i_shell = 0;

    // find out to which shell this atom belongs to
    // TODO: replace with binary search
    while(i_shell < track.size() && track[i_shell][0] <= r_particle) i_shell++;
    if(i_shell > 0) {
      i_shell--; 
    }
    else {
      cerr << "track.in radii should start from 0\n";
      exit(1);
    }

    double v2 = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
    double r_relative = (r_particle - track[i_shell][0]) / (track[i_shell+1][0]-track[i_shell][0]);

    int atom_column = type[i];
    double e_add = (track[i_shell][atom_column]*(1 - r_relative) + track[i_shell][atom_column]*(r_relative));

    // get a random direction in 3d
    double phi = myranf(gen) * 2 * M_PI;
    double theta = acos(cos(M_PI) + myranf(gen) * (1 - cos(M_PI)));
     
    // adding this velocity gives on average the right kinetic energy change
    // proof: compute the magnitude of \vec \delta v, 2 cos(alpha) term is on average 0 
    // sqrt(eV / u)= 98.2269479 angstrom / ps
    double v_new = sqrt(2*e_add/(atom->mass[type[i]]))*98.2269479; 
    
    v[i][0] += v_new*sin(theta)*cos(phi);
    v[i][1] += v_new*sin(theta)*sin(phi);
    v[i][2] += v_new*cos(theta);

  }
  current_time += dt;
}

/* ---------------------------------------------------------------------- */
double FixTrackIn::compute_scalar()
{
  // only sum across procs when changed since last call
  return 1.0;
}
/* ---------------------------------------------------------------------- */

