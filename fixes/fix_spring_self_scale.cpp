// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Fix spring/self/scale: like fix spring/self, but the spring constant k
   may be an equal-style variable that is re-evaluated every step, so the
   spring can be scaled continuously during a thermodynamic integration
   (e.g. k = k0 * (1 - lambda)).

   Based on fix_spring_self.cpp by Naveen Michaud-Agrawal (JHU).
------------------------------------------------------------------------- */

#include "fix_spring_self_scale.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpringSelfScale::FixSpringSelfScale(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  kstr(nullptr),
  xoriginal(nullptr)
{
  if ((narg < 4) || (narg > 6))
    error->all(FLERR,"Illegal fix spring/self/scale command");

  restart_peratom = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;

  // k: either a constant or an equal-style variable referenced as v_name

  kvar = -1;
  kconst = 0.0;
  k = 0.0;

  if (utils::strmatch(arg[3],"^v_")) {
    kstr = utils::strdup(arg[3]+2);
    kvar = -2;                       // resolved in init()
  } else {
    kconst = utils::numeric(FLERR,arg[3],false,lmp);
    if (kconst < 0.0) error->all(FLERR,"Illegal fix spring/self/scale command");
    k = kconst;
  }

  xflag = yflag = zflag = 1;
  probe = 0;

  for (int i = 4; i < narg; i++) {
    if (strcmp(arg[i],"probe") == 0) {
      probe = 1;
    } else if (strcmp(arg[i],"xyz") == 0) {
      xflag = yflag = zflag = 1;
    } else if (strcmp(arg[i],"xy") == 0) {
      zflag = 0;
    } else if (strcmp(arg[i],"xz") == 0) {
      yflag = 0;
    } else if (strcmp(arg[i],"yz") == 0) {
      xflag = 0;
    } else if (strcmp(arg[i],"x") == 0) {
      yflag = zflag = 0;
    } else if (strcmp(arg[i],"y") == 0) {
      xflag = zflag = 0;
    } else if (strcmp(arg[i],"z") == 0) {
      xflag = yflag = 0;
    } else error->all(FLERR,"Illegal fix spring/self/scale command");
  }

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = nullptr;
  FixSpringSelfScale::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // xoriginal = initial unwrapped positions of atoms

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  espring = 0.0;
}

/* ---------------------------------------------------------------------- */

FixSpringSelfScale::~FixSpringSelfScale()
{
  if (copymode) return;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  // delete locally stored array

  memory->destroy(xoriginal);
  delete[] kstr;
}

/* ---------------------------------------------------------------------- */

int FixSpringSelfScale::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::init()
{
  // resolve the spring-constant variable

  if (kstr) {
    kvar = input->variable->find(kstr);
    if (kvar < 0)
      error->all(FLERR,"Variable name for fix spring/self/scale does not exist");
    if (!input->variable->equalstyle(kvar))
      error->all(FLERR,"Variable for fix spring/self/scale must be equal-style");
    k = input->variable->compute_equal(kvar);
    if (k < 0.0) error->all(FLERR,"Illegal fix spring/self/scale command");
  }

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // re-evaluate the spring constant if it is a variable

  if (kstr) {
    k = input->variable->compute_equal(kvar);
    if (k < 0.0) error->all(FLERR,"Fix spring/self/scale variable became negative");
  }

  double dx,dy,dz;
  double unwrap[3];

  espring = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xoriginal[i][0];
      dy = unwrap[1] - xoriginal[i][1];
      dz = unwrap[2] - xoriginal[i][2];
      if (!xflag) dx = 0.0;
      if (!yflag) dy = 0.0;
      if (!zflag) dz = 0.0;
      if (!probe) {                 // probe = tally energy only, no force
        f[i][0] -= k*dx;
        f[i][1] -= k*dy;
        f[i][2] -= k*dz;
      }
      espring += k * (dx*dx + dy*dy + dz*dz);
    }

  espring *= 0.5;
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringSelfScale::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched springs
------------------------------------------------------------------------- */

double FixSpringSelfScale::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSpringSelfScale::memory_usage()
{
  double bytes = (double)atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSpringSelfScale::grow_arrays(int nmax)
{
  memory->grow(xoriginal,nmax,3,"fix_spring/self/scale:xoriginal");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixSpringSelfScale::copy_arrays(int i, int j, int /*delflag*/)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixSpringSelfScale::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixSpringSelfScale::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixSpringSelfScale::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixSpringSelfScale::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixSpringSelfScale::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixSpringSelfScale::size_restart(int /*nlocal*/)
{
  return 4;
}
