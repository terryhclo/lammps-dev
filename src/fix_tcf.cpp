#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_tcf.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "group.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "memory.h"
#include "universe.h"
#include "domain.h"
#include "mpi.h"
#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define SIZE 3
#define MAX_INT 2147483647
#define METHOD 1

using namespace LAMMPS_NS;

FixTCF::FixTCF(LAMMPS* lmp, int narg, char** arg) :
  Fix(lmp, narg, arg)
{
   if (narg < 6 || (narg-6) % 2) error->all("Illegal fix tcf command");
   if (narg>6) options(narg,arg);
   else error->all("no file name specified <file filename>");
   int mintag = MAX_INT;
   int i;
   int *mpigrpset;
   int nlocal = atom->nlocal;
   int nproc=universe->nprocs;
   int ntypes=atom->ntypes;
   bigflag=0;scalar_flag = 1;global_freq=1;
   nevery=atoi(arg[3]);                    //mendatory by FIX
   t_start=atoi(arg[4]);                   //When to start c(0)
   tau=atoi(arg[5]);                       //#time step for each c(t) point
   type = atom->type;                      //array with type numbers
   tcftype=atoi(arg[1]);                   //type of atm we wnt its tcf(is not grp)
   time_depend=1;                          //this FIX is time dependent
   totn=group->count(igroup);              //How many atoms are in this group
   intnpair=(totn*(totn-1)/2);
   dnpair=(double) ((totn*(totn-1)/2));
   share=intnpair/nproc;                   //each proc. share of total #of pairs
   rank=universe->me;
   if (share<1) share=1;              //check if too many proc compr to ttal pairs
   jstart=rank*share+1;               //relating start index to pair indices
   istart=0;                          //relating start index to pair indices
   setij(&istart,&jstart);            //relating start index to pair indices
   if (rank==nproc-1) myshare=share+intnpair%nproc;
   else myshare=share;
   if (ntypes>2)
       fprintf(screen,"Warning:This code is not yet debugged for more than 2 types\n");
   for ( i = 0; i < nlocal; i++)         //Finding the min tag of type tcf type
        if (type[i]==tcftype) mintag = MIN(mintag,atom->tag[i]);
   MPI_Allreduce(&mintag,&mintag_all,1,MPI_INT,MPI_MIN,world);
   memory->create(xyz,totn,3,"tcf:xyz");  //where to save xyz of 0 to tau
   if (mthd==0) memory->create(rij,myshare,tau/cg+1,"tcf:rij");         //rij for #pairs in last proc/memory saving method but inaccurate
   if (rank==0) fprintf(screen,"\ntau: %d MinTag: %d totl atm: %d method %d CG= %d\n",tau,mintag_all,totn,mthd,cg);
}

FixTCF::~FixTCF()
{
 /*MPI_Group_free(&tcfgrp)*/
 delete type;
 delete [] ct;
 delete [] ave;
 delete [] fct;
 memory->destroy(rij);
 memory->destroy(xyz);
 if (fpt && rank == 0) fclose(fpt);
}

void FixTCF::init()
{
   int final=update->endstep;
   if (final<t_start) return;
   if ((final-t_start)%tau!=0) error->all("remainder of TCF time to tau should be zero");
   if (tau%cg!=0) error->all("remainder of tau to CG parameter should be zero");
   if (mthd==1)  
           memory->create(rij,myshare,(final-t_start)/cg,"tcf:rij");         //rij for #pairs in last proc
   ct =(double *)calloc((final-t_start)/tau,sizeof(double));
   ave=(double *)calloc((final-t_start)/tau,sizeof(double));
   fct =(double *)calloc((final-t_start)/tau,sizeof(double));
   if (ct==NULL || ave==NULL) fprintf(screen,"Unable to init array ct or ave at fix tcf/n");
   if (rank==0){
       fprintf(screen,"\n END of TCF INIT, fpt pointer %p, col %d\n", fpt,(final-t_start)/tau);
       if (fputs("Testing...\n",fpt)==EOF) error->one("intializong file error");
       fprintf(fpt,"Time Auto Correlation Data\n");
       fflush(fpt);
   }
}

int FixTCF::setmask(){
   int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

void FixTCF::setup(int vflag)
{
  end_of_step();
}

void FixTCF::end_of_step()
{
  int i,m,t_bin;
  int ntimestep = update->ntimestep;
  int endstep=update->endstep;
  double **x = atom->x;
  int *tag=atom->tag;
  int nlocal = atom->nlocal;
  int *mask = atom->type;
  int now=(ntimestep-t_start);
  if (ntimestep<t_start || ntimestep>=endstep ) {
     return;
  }
  else{
      t_bin=(ntimestep-t_start)/tau;
      if (now%cg==0){
          now=now/cg;
          for (i=0;i<totn;i++)  {
          xyz[i][0]=0;xyz[i][1]=0;xyz[i][2]=0;
          }
          for (i=0;i<nlocal;i++){
              m=tag[i]-mintag_all;
              if (type[i]==tcftype){
                  xyz[m][0]=x[i][0];
                  xyz[m][1]=x[i][1];
                  xyz[m][2]=x[i][2];
              }
          }
      MPI_Barrier(world);//Do I need it??
      if (universe->nprocs>1)  MPI_Allreduce(&xyz[0][0],&xyz[0][0],3*totn,MPI_DOUBLE,MPI_SUM,world);
      endofstep1(t_bin,now,jstart); 
      }
      if (ntimestep==endstep-1) {
          MPI_Barrier(world);
          printdata(t_bin);
      }
  }
}

void FixTCF::endofstep1(int t_bin, int now,int dummyj){
    int i,j,k;
    int cell=now;
    int kstart=0;
    int counter=0;
    double temp=0;
    for (i=istart-1;i<totn;i++){
          for (j=dummyj-1;j<totn;j++){
              if (t_bin!=0 && mthd==0) cell=tau/cg;
              rcalc(&temp,i,j);                                  //cell is either tau now
              rij[0][cell]+=temp;
              counter++;
              if (counter>=myshare) break;
          }
          dummyj=i+3;
          if (counter>=myshare) break;
    }
    MPI_Reduce(&rij[0][cell],&rij[1][cell],1,MPI_DOUBLE,MPI_SUM,0,world);
    if (rank==0){
        rij[1][cell]/=dnpair;
        if (mthd==0) kstart=t_bin;
        for (k=kstart;k<t_bin+1;k++){
              ave[k]+=rij[1][cell];
              ct[k]+=(rij[1][cell]*rij[1][now-k*tau/cg]);
        }
    }
    else rij[1][cell]=0;
}

void FixTCF::setij( int *i,int *in)
{
    (*i)++;
    *in=*in+*i-totn;
    if (*in>0 || (*i)<totn) setij(i,in);
    if (*in<=0) *in=*in+totn;
}

void FixTCF::rcalc(double *a, int i,int j)
{
    double dx,dy,dz;
    dx=fabs(xyz[i][0]-xyz[j][0]);
    dy=fabs(xyz[i][1]-xyz[j][1]);
    dz=fabs(xyz[i][2]-xyz[j][2]);
    if (domain->xperiodic && dx>domain->xprd_half) dx=domain->xprd-dx;
    if (domain->yperiodic && dy>domain->yprd_half) dy=domain->xprd-dy;
    if (domain->zperiodic && dz>domain->zprd_half) dz=domain->xprd-dz;
    *a=sqrt(dx*dx+dy*dy+dz*dz);
}

double FixTCF::compute_scalar()
{
    return ave[0];
}

void FixTCF::options(int narg, char **arg)
{
    mthd=1;
    cg=1;
    if (rank==0) fpt=NULL;
    int iarg=6;
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix tcf command");
      if (rank == 0) {
	fpt = fopen(arg[iarg+1],"w"); //arg[iarg+1]
	if (fpt == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix tcf file %s",arg[iarg+1]);
	  error->one(str);
	}
      }
    }
    iarg+=2;
    if (iarg>=narg)  return;
    if (strcmp(arg[iarg],"method") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix tcf command");
      mthd=atoi(arg[iarg+1]);
    }
    iarg+=2;
    if (iarg>=narg) return;
    if (strcmp(arg[iarg],"CG") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix tcf command");
      cg=atoi(arg[iarg+1]);
      if (cg<1) cg=1;
    }
}

void FixTCF::printdata(int t_bin)
{
    int endstep=update->endstep;
    int ntcfstep=endstep-t_start;
    int tcg=tau/cg;
    double tave=0;
    double *nave;
    double dt=tau*update->dt/1000.0;
    if (!fpt || rank!=0) return;
    nave=(double *)calloc((endstep-t_start)/tau,sizeof(double));
    for (int i=0;i<t_bin;i++) nave[t_bin-i]=(ave[0]-ave[i+1])/(i+1)/tcg;
    nave[0]=ave[0]/(ntcfstep)*cg;
    fprintf(fpt,"time unit as unite[dt]/1000 (r.g ns)\n");
    for (int i=0;i<t_bin+1;i++) {
          fprintf(fpt,"raw c(t): %f  raw ave: %f   ",ct[i],ave[i]);
          if (mthd==1) ct[i]=ct[i]/(ntcfstep-i*tau)*cg;
          if (mthd==1) ave[i]=ave[i]/(ntcfstep-i*tau)*cg;
          if (mthd==0) ct[i]=ct[i]/tcg;
          if (mthd==0) ave[i]=ave[i]/tcg;
          tave=tave+ave[i];
          fprintf(fpt,"t: %f  nonnormal ct(t): %f c(t): %f ave(t): %f \n",i*dt,ct[i],ct[i]/ct[0],ave[i]);
    }
    if (mthd==0) tave=tave/(t_bin+1);
    else tave=ave[0];
    fprintf(fpt,"Auto Time Correlation Main Dataset \n");
    fprintf(fpt,"i  | 2nd Ave(nave) |  Time  |  C(t) \n");
    fprintf(fpt,"------------------------------ \n");
    //testtcf(tave,ntcfstep,dt);
    for(int i=0;i<t_bin+1;i++) {
        fprintf(fpt,"%d    %f ",i, nave[i]);
        if (mthd==1) fprintf(fpt, " %f  %f\n",i*dt,(ct[i]+tave*tave-tave*ave[i]-tave*nave[i])/(ct[0]-tave*tave));
        if (mthd==0) fprintf(fpt, " %f  %f\n",i*dt,(ct[i]+tave*tave-tave*ave[0]-tave*ave[i])/(ct[0]+tave*tave-2*tave*ave[0]));
    }
    fprintf(fpt, "average distance %f \n",tave);
    
    if (fflush(fpt)!=0) fprintf(screen,"Error in fflush \n");
    delete [] nave;
}

//DEbug porpous only//
int FixTCF::testtcf(double av, int ntimestep, double dt)
{
    int i,t_bin,k;
    int kstart=0;
    for (i=0;i<ntimestep;i+=cg){
        t_bin=i/tau;
        if (mthd==0) kstart=t_bin;
        for (k=kstart;k<t_bin+1;k++){
              fct[k]+=(rij[1][i/cg]-av)*(rij[1][i/cg-k*tau/cg]-av);
        }
    }
     fprintf(fpt, "----TEST TCF WITH CONVENTIONAL-------------- \n");
        fprintf(fpt,"t   nonnormal fct(t)   fc(t) \n");
        for (int i=0;i<t_bin+1;i++) {
            if (mthd==1) fct[i]=fct[i]/(ntimestep-i*tau)*cg;
            if (mthd==0) fct[i]=fct[i]/(tau/cg);
            fprintf(fpt," %f  %f  %f \n",i*dt,fct[i]*1000,fct[i]/fct[0]);
        }
    return k;
}
