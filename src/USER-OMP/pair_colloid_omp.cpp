/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_colloid_omp.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairColloidOMP::PairColloidOMP(LAMMPS *lmp) :
  PairColloid(lmp), ThrOMP(lmp, PAIR)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairColloidOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(this);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int ifrom, ito, tid;
    double **f;

    f = loop_setup_thr(atom->f, ifrom, ito, tid, inum, nall, nthreads);

    if (evflag) {
      if (eflag) {
	if (force->newton_pair) eval<1,1,1>(f, ifrom, ito, tid);
	else eval<1,1,0>(f, ifrom, ito, tid);
      } else {
	if (force->newton_pair) eval<1,0,1>(f, ifrom, ito, tid);
	else eval<1,0,0>(f, ifrom, ito, tid);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(f, ifrom, ito, tid);
      else eval<0,0,0>(f, ifrom, ito, tid);
    }

    // reduce per thread forces into global force array.
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
  } // end of omp parallel region

  // reduce per thread energy and virial, if requested.
  if (evflag) ev_reduce_thr(this);
  if (vflag_fdotr) virial_fdotr_compute();
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairColloidOMP::eval(double **f, int iifrom, int iito, int tid)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,r2inv,r6inv,forcelj,factor_lj;
  double c1,c2,fR,dUR,dUA,K[9],h[4],g[4];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;
      
      switch(form[itype][jtype]) {
      case SMALL_SMALL:
	r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	fpair = factor_lj*forcelj*r2inv;
	if (EFLAG) 
	  evdwl = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
	    offset[itype][jtype];
	break;

      case SMALL_LARGE:
	c2 = a2[itype][jtype];
	K[1] = c2*c2;
	K[2] = rsq;
	K[0] = K[1] - rsq;
	K[4] = rsq*rsq;
	K[3] = K[1] - K[2];
	K[3] *= K[3]*K[3];
	K[6] = K[3]*K[3];
	fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
	fpair = 4.0/15.0*fR*factor_lj * 
	  (2.0*(K[1]+K[2]) * (K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) * 
	   sigma6[itype][jtype]/K[6]-5.0) / K[0];
	if (EFLAG) 
	  evdwl = 2.0/9.0*fR * 
	    (1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) *
	     sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
	if (rsq <= K[1]) error->one(FLERR,"Overlapping small/large in pair colloid");
	break;

      case LARGE_LARGE:
	r = sqrt(rsq);
	c1 = a1[itype][jtype];
	c2 = a2[itype][jtype];
	K[0] = c1*c2;
	K[1] = c1+c2;
	K[2] = c1-c2;
	K[3] = K[1]+r;
	K[4] = K[1]-r;
	K[5] = K[2]+r;
	K[6] = K[2]-r;
	K[7] = 1.0/(K[3]*K[4]);
	K[8] = 1.0/(K[5]*K[6]);
	g[0] = pow(K[3],-7.0);
	g[1] = pow(K[4],-7.0);
	g[2] = pow(K[5],-7.0);
	g[3] = pow(K[6],-7.0);
	h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
	h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
	h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
	h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
	g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
	g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
	g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
	g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];
	
	fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
	evdwl = fR * (h[0]-h[1]-h[2]+h[3]);
	dUR = evdwl/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
	dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] + 
					(2.0*K[0]*K[8]-1.0)*K[8]);
	fpair = factor_lj * (dUR+dUA)/r;
	if (EFLAG)
	  evdwl += a12[itype][jtype]/6.0 * 
	    (2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) - offset[itype][jtype];
	if (r <= K[1]) error->one(FLERR,"Overlapping large/large in pair colloid");
	break;
      }
      
      if (EFLAG) evdwl *= factor_lj;
    
      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      if (NEWTON_PAIR || j < nlocal) {
	f[j][0] -= delx*fpair;
	f[j][1] -= dely*fpair;
	f[j][2] -= delz*fpair;
      }

      if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
			       evdwl,0.0,fpair,delx,dely,delz,tid);
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairColloidOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairColloid::memory_usage();

  return bytes;
}
