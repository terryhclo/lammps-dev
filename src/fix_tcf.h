/* 
 * File:   fix_tcf.h
 * Author: Mahdy
 *
 * Usage syntax:
 *
 * fix [fix_id] [atom_type] [style=tcf] [nevery=1] [starting_time_step] [tau]  &
 *     <file filename>
 *
 * e.g.
 *
 * fix tcor 2 tcf 1 1800000 50000 file TCF03.txt
 *
 *
 * Created on April 13, 2010, 1:18 AM
 */

#ifdef FIX_CLASS

FixStyle(tcf,FixTCF)

#else

#ifndef LMP_FIX_TCF_H
#define	LMP_FIX_TCF_H

#include "fix.h"

namespace LAMMPS_NS {

 class FixTCF : public Fix {

 public:

  FixTCF(class LAMMPS *, int, char **);

  ~FixTCF();

  void   init();
  int    setmask();
  void   setup(int);
  void   end_of_step();
  double compute_scalar();
  int    bigflag;

 protected:

  FILE *fpt;
  int rank,atype,tcftype,totn,share,myshare; //-,atom_type,type to find it's tcf
  int t_start,tau,t_current,mthd,cg;
  int mintag_all;
  int jstart,istart,intnpair;
  int *type;                            //pointer to atom->type
  double **rij,**xyz,*ct,*fct,*ave,*fave,dnpair;
  double *n_ave,*n_fave,*n_nave;        //new averages
  class Compute *local_pair;
  class Compute *local_dist;

     /*MPI_Datatype anxyz;
      MPI_Group tcfgrp;
      MPI_Comm  tcfworld;*/

  void options(int, char **);
  void setij(int *,int *);
  void rcalc(double *, int, int);
  void endofstep1(int, int , int);
  void printdata(int);
  int testtcf(double,int, double);

    };

}

#endif	/* _FIX_TCF_H */
#endif

