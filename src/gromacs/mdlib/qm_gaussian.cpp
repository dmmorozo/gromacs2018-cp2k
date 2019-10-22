/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_QMMM_GAUSSIAN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include <libcp2k.h>

/* TODO: this should be made thread-safe */

/* Gaussian interface routines */

void write_cp2k_input_pdb(t_QMrec *qm, t_MMrec *mm)
{
    char
          periodic_system[37][4]={"X  ","H  ","He ","Li ","Be ","B  ","C  ","N  ",
			"O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ",
			"S  ","Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ",
			"Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ",
			"Ge ","As ","Se ","Br ","Kr "};
    FILE *cp2kinp = NULL;
    int
          i,j,k;

    //Write CP2K input file
    cp2kinp = fopen("cp2k.inp", "w");
    fprintf(cp2kinp,"&GLOBAL\n");
    fprintf(cp2kinp,"  PRINT_LEVEL LOW\n");
    fprintf(cp2kinp,"  PROJECT H2O\n");
    fprintf(cp2kinp,"  RUN_TYPE ENERGY_FORCE\n"); 
    fprintf(cp2kinp,"&END GLOBAL\n");
   
    fprintf(cp2kinp,"&FORCE_EVAL\n");
    fprintf(cp2kinp,"  METHOD QMMM\n"); 
    fprintf(cp2kinp,"  &DFT\n");
    fprintf(cp2kinp,"    CHARGE %d\n",qm->QMcharge);
    fprintf(cp2kinp,"    MULTIPLICITY %d\n",qm->multiplicity);
    fprintf(cp2kinp,"    BASIS_SET_FILE_NAME  BASIS_MOLOPT\n");
    fprintf(cp2kinp,"    POTENTIAL_FILE_NAME  POTENTIAL\n");
    fprintf(cp2kinp,"    &MGRID\n");
    fprintf(cp2kinp,"      CUTOFF 400\n");
    fprintf(cp2kinp,"      COMMENSURATE\n");
    fprintf(cp2kinp,"    &END MGRID\n");
    fprintf(cp2kinp,"    &SCF\n");
    fprintf(cp2kinp,"      SCF_GUESS RESTART\n");
    fprintf(cp2kinp,"      EPS_SCF 1.0E-5\n");
    fprintf(cp2kinp,"      MAX_SCF 300\n");
    fprintf(cp2kinp,"      &OT  T\n");
    fprintf(cp2kinp,"        MINIMIZER  DIIS\n");
    fprintf(cp2kinp,"        STEPSIZE     1.4999999999999999E-01\n");
    fprintf(cp2kinp,"        PRECONDITIONER FULL_ALL\n");
    fprintf(cp2kinp,"      &END OT\n");
    fprintf(cp2kinp,"    &END SCF\n");
    fprintf(cp2kinp,"    &XC\n");
    fprintf(cp2kinp,"      &XC_FUNCTIONAL PBE\n");
    fprintf(cp2kinp,"      &END XC_FUNCTIONAL\n");
    fprintf(cp2kinp,"    &END XC\n");
    fprintf(cp2kinp,"    &QS\n");
    fprintf(cp2kinp,"     METHOD GPW\n");
    fprintf(cp2kinp,"     EPS_DEFAULT 1.0E-10\n");
    fprintf(cp2kinp,"     EXTRAPOLATION ASPC\n");
    fprintf(cp2kinp,"     EXTRAPOLATION_ORDER  4\n");

    fprintf(cp2kinp,"    &END QS\n");
    fprintf(cp2kinp,"  &END DFT\n");
    
    fprintf(cp2kinp,"  &QMMM\n");
    fprintf(cp2kinp,"    &CELL\n");

    //calculate maximum size of QM system
    double dimx, dimy, dimz, maxdim = 0.0;
    for(i=0;i<qm->nrQMatoms;i++)
       for(j=i+1;j<qm->nrQMatoms;j++)
    {
        dimx=sqrt((qm->xQM[j][0]-qm->xQM[i][0])*(qm->xQM[j][0]-qm->xQM[i][0]) + (qm->xQM[j][1]-qm->xQM[i][1])*(qm->xQM[j][1]-qm->xQM[i][1]) + (qm->xQM[j][2]-qm->xQM[i][2])*(qm->xQM[j][2]-qm->xQM[i][2]));
        maxdim = (maxdim < dimx) ? dimx : maxdim;
    }

    //increase the dimensions and check if they are larger than full box
    dimx = (1.5*maxdim > qm->box[0][0]) ? qm->box[0][0] : 1.5*maxdim;
    dimy = (1.5*maxdim > qm->box[1][1]) ? qm->box[1][1] : 1.5*maxdim;
    dimz = (1.5*maxdim > qm->box[2][2]) ? qm->box[2][2] : 1.5*maxdim;
    

    fprintf(stderr,"Cell dimension of QM system= %.3lf x %.3lf x %.3lf nm\n", dimx, dimy, dimz);

    //Label for further modifications
    fprintf(cp2kinp,"#QMCELL\n"); 
 
    fprintf(cp2kinp,"      ABC %.3lf %.3lf %.3lf\n", dimx*10, dimy*10, dimz*10); 
    fprintf(cp2kinp,"      PERIODIC XYZ\n");
    fprintf(cp2kinp,"    &END CELL\n");
    fprintf(cp2kinp,"    ECOUPL GAUSS\n"); 
    fprintf(cp2kinp,"    USE_GEEP_LIB 8\n");  
    fprintf(cp2kinp,"    &PERIODIC\n");
    fprintf(cp2kinp,"      GMAX     1.0000000000000000E+00\n");
 
    fprintf(cp2kinp,"      &MULTIPOLE ON\n");
    fprintf(cp2kinp,"         RCUT     1.0000000000000007E+01\n"); 
    fprintf(cp2kinp,"         EWALD_PRECISION     9.9999999999999995E-07\n");

    fprintf(cp2kinp,"      &END\n");
    fprintf(cp2kinp,"    &END PERIODIC\n");

    //search for all QM Hydrogenes
    k = 0;
    for(i=0;i<qm->nrQMatoms;i++)
       if (qm->atomicnumberQM[i]==1) k++;
    if (k>0) {
       fprintf(cp2kinp,"    &QM_KIND H\n");
       fprintf(cp2kinp,"      MM_INDEX");
       for(i=0;i<qm->nrQMatoms;i++)
          if (qm->atomicnumberQM[i]==1) fprintf(cp2kinp," %d",i+1);
       fprintf(cp2kinp,"\n");
       fprintf(cp2kinp,"    &END QM_KIND\n");
    }

    //search for all QM Carbons
    k = 0;
    for(i=0;i<qm->nrQMatoms;i++)
       if (qm->atomicnumberQM[i]==6) k++;
    if (k>0) {
       fprintf(cp2kinp,"    &QM_KIND C\n");
       fprintf(cp2kinp,"      MM_INDEX");
       for(i=0;i<qm->nrQMatoms;i++)
          if (qm->atomicnumberQM[i]==6) fprintf(cp2kinp," %d",i+1);
       fprintf(cp2kinp,"\n");
       fprintf(cp2kinp,"    &END QM_KIND\n");
    }

    //search for all QM Nitrogenes
    k = 0;
    for(i=0;i<qm->nrQMatoms;i++)
       if (qm->atomicnumberQM[i]==7) k++;
    if (k>0) {
       fprintf(cp2kinp,"    &QM_KIND N\n");
       fprintf(cp2kinp,"      MM_INDEX");
       for(i=0;i<qm->nrQMatoms;i++)
          if (qm->atomicnumberQM[i]==7) fprintf(cp2kinp," %d",i+1);
       fprintf(cp2kinp,"\n");
       fprintf(cp2kinp,"    &END QM_KIND\n");
    }

    //search for all QM Oxygenes
    k = 0;
    for(i=0;i<qm->nrQMatoms;i++)
       if (qm->atomicnumberQM[i]==8) k++;
    if (k>0) {
       fprintf(cp2kinp,"    &QM_KIND O\n");
       fprintf(cp2kinp,"      MM_INDEX");
       for(i=0;i<qm->nrQMatoms;i++)
          if (qm->atomicnumberQM[i]==8) fprintf(cp2kinp," %d",i+1);
       fprintf(cp2kinp,"\n");
       fprintf(cp2kinp,"    &END QM_KIND\n");
    }

    fprintf(cp2kinp,"  &END QMMM\n");
    
    fprintf(cp2kinp,"  &MM\n");
    fprintf(cp2kinp,"    &FORCEFIELD\n");
    fprintf(cp2kinp,"      DO_NONBONDED FALSE\n");
    fprintf(cp2kinp,"    &END FORCEFIELD\n");
    fprintf(cp2kinp,"    &POISSON\n");
    fprintf(cp2kinp,"      &EWALD\n");
    fprintf(cp2kinp,"        EWALD_TYPE NONE\n");
    fprintf(cp2kinp,"      &END EWALD\n");
    fprintf(cp2kinp,"    &END POISSON\n");
    fprintf(cp2kinp,"  &END MM\n");
    
    fprintf(cp2kinp,"  &SUBSYS\n");

    fprintf(cp2kinp,"    &CELL\n");

    //Label for further modifiactions
    fprintf(cp2kinp,"#CELL\n");
    
    fprintf(cp2kinp,"      A %.3lf %.3lf %.3lf\n",qm->box[0][0]*10,qm->box[0][1]*10,qm->box[0][2]*10); 
    fprintf(cp2kinp,"      B %.3lf %.3lf %.3lf\n",qm->box[1][0]*10,qm->box[1][1]*10,qm->box[1][2]*10); 
    fprintf(cp2kinp,"      C %.3lf %.3lf %.3lf\n",qm->box[2][0]*10,qm->box[2][1]*10,qm->box[2][2]*10); 
    
    fprintf(cp2kinp,"      PERIODIC XYZ\n");
    fprintf(cp2kinp,"    &END CELL\n");
    fprintf(cp2kinp,"    &TOPOLOGY\n");
    fprintf(cp2kinp,"      COORD_FILE_NAME cp2k.pdb\n"); 
    fprintf(cp2kinp,"      COORD_FILE_FORMAT PDB\n");
    fprintf(cp2kinp,"      CHARGE_EXTENDED TRUE\n");    
    fprintf(cp2kinp,"      CONNECTIVITY OFF\n");
    fprintf(cp2kinp,"      &GENERATE\n");
    fprintf(cp2kinp,"         &ISOLATED_ATOMS\n");
    
    fprintf(cp2kinp,"            LIST %d..%d\n",1,qm->nrQMatoms+mm->nrMMatoms);
    
    fprintf(cp2kinp,"         &END\n");
    fprintf(cp2kinp,"      &END GENERATE\n");
    fprintf(cp2kinp,"    &END TOPOLOGY\n");
    fprintf(cp2kinp,"    &KIND H\n");
    fprintf(cp2kinp,"      ELEMENT H\n");
    fprintf(cp2kinp,"      BASIS_SET DZVP-MOLOPT-GTH\n");
    fprintf(cp2kinp,"      POTENTIAL GTH-PBE-q1\n");
    fprintf(cp2kinp,"    &END KIND\n");
    fprintf(cp2kinp,"    &KIND C\n");

    fprintf(cp2kinp,"      ELEMENT C\n");
    fprintf(cp2kinp,"      BASIS_SET DZVP-MOLOPT-GTH\n");
    fprintf(cp2kinp,"      POTENTIAL GTH-PBE-q4\n");
    fprintf(cp2kinp,"    &END KIND\n");
    fprintf(cp2kinp,"    &KIND N\n");
    fprintf(cp2kinp,"      ELEMENT N\n");
    fprintf(cp2kinp,"      BASIS_SET DZVP-MOLOPT-GTH\n");
    fprintf(cp2kinp,"      POTENTIAL GTH-PBE-q5\n");
    fprintf(cp2kinp,"    &END KIND\n");
    fprintf(cp2kinp,"    &KIND O\n");
    fprintf(cp2kinp,"      ELEMENT O\n");
    fprintf(cp2kinp,"      BASIS_SET DZVP-MOLOPT-GTH\n");
    fprintf(cp2kinp,"      POTENTIAL GTH-PBE-q6\n");
    fprintf(cp2kinp,"    &END KIND\n");
    fprintf(cp2kinp,"    &KIND X\n");
    fprintf(cp2kinp,"      ELEMENT H\n");
    fprintf(cp2kinp,"      BASIS_SET DZVP-MOLOPT-GTH\n");
    fprintf(cp2kinp,"      POTENTIAL GTH-PBE-q1\n");
    fprintf(cp2kinp,"    &END KIND\n");
    fprintf(cp2kinp,"  &END SUBSYS\n");
    fprintf(cp2kinp,"&END FORCE_EVAL\n");
    fclose(cp2kinp);

    //write PDB file
    cp2kinp = fopen("cp2k.pdb", "w");

    //write QM atoms
    for(i=0;i<qm->nrQMatoms;i++) {
       fprintf(cp2kinp,"ATOM ");
       fprintf(cp2kinp,"%6d ",i+1);
       fprintf(cp2kinp," %3s ",periodic_system[qm->atomicnumberQM[i]]);
       fprintf(cp2kinp," QM");
       fprintf(cp2kinp,"%6d     ",1);
       fprintf(cp2kinp,"%7.3lf %7.3lf %7.3lf  1.00  0.00         ",qm->xQM[i][0]*10,qm->xQM[i][1]*10,qm->xQM[i][2]*10);
       fprintf(cp2kinp," %3s ",periodic_system[qm->atomicnumberQM[i]]);
       fprintf(cp2kinp,"%.2lf\n",qm->atomicnumberQM[i]/1.0);
    }

    //write MM atoms
    for(i=0;i<mm->nrMMatoms;i++) {
       fprintf(cp2kinp,"ATOM ");
       fprintf(cp2kinp,"%6d ",qm->nrQMatoms+i+1);
       fprintf(cp2kinp," %3s ",periodic_system[0]);
       fprintf(cp2kinp," MM");
       fprintf(cp2kinp,"%6d     ",2);
       fprintf(cp2kinp,"%7.3lf %7.3lf %7.3lf  1.00  0.00         ",mm->xMM[i][0]*10,mm->xMM[i][1]*10,mm->xMM[i][2]*10);
       fprintf(cp2kinp," %3s ",periodic_system[0]);
       fprintf(cp2kinp,"%lf\n",mm->MMcharges[i]);
    }
    fprintf(cp2kinp,"END\n");

          
    fclose(cp2kinp);


}

void update_cp2k_input(t_QMrec *qm, t_MMrec *mm)
{
    char
          periodic_system[37][4]={"X  ","H  ","He ","Li ","Be ","B  ","C  ","N  ",
			          "O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ",
			          "S  ","Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ",
			          "Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ",
			          "Ge ","As ","Se ","Br ","Kr "};
    FILE *cp2kinp = NULL, *tmp = NULL;
    int  i,j,k;
    char buf[300];
    double dimx, dimy, dimz, maxdim = 0.0;

    //open CP2K input file
    if (qm->cp2k_inp)
    {
        cp2kinp = fopen(qm->cp2k_inp, "r");
    }
    else
    {
        cp2kinp = fopen("cp2k.inp", "r");
    }

    if (!cp2kinp)
    {
        gmx_fatal(FARGS, "Could not find CP2K input file, check CP2K_INPUT environmental variable\n");
    }


    //calculate maximum size of QM system
    for(i=0;i<qm->nrQMatoms;i++)
       for(j=i+1;j<qm->nrQMatoms;j++)
    {
        dimx=sqrt((qm->xQM[j][0]-qm->xQM[i][0])*(qm->xQM[j][0]-qm->xQM[i][0]) + (qm->xQM[j][1]-qm->xQM[i][1])*(qm->xQM[j][1]-qm->xQM[i][1]) + (qm->xQM[j][2]-qm->xQM[i][2])*(qm->xQM[j][2]-qm->xQM[i][2]));
        maxdim = (maxdim < dimx) ? dimx : maxdim;
    }

    //increase the dimensions and check if they are larger than full box
    dimx = (1.5*maxdim > qm->box[0][0]) ? qm->box[0][0] : 1.5*maxdim;
    dimy = (1.5*maxdim > qm->box[1][1]) ? qm->box[1][1] : 1.5*maxdim;
    dimz = (1.5*maxdim > qm->box[2][2]) ? qm->box[2][2] : 1.5*maxdim;

    //Open temporary input file
    tmp = fopen ("cp2k.tmp.inp","w");

    fgets(buf, 300, cp2kinp);
    
    while (!feof(cp2kinp))
    {
        //if the next line contains cell dimensions then update them
        if (strncmp(buf,"#CELL",5) == 0)
        {
            fputs(buf, tmp);
            //A
            fgets(buf, 300, cp2kinp);
            //B
            fgets(buf, 300, cp2kinp);
            //C
            fgets(buf, 300, cp2kinp);
            //Wrtie new cell dimensions instead
            fprintf(tmp,"      A %.3lf %.3lf %.3lf\n",qm->box[0][0]*10,qm->box[0][1]*10,qm->box[0][2]*10); 
            fprintf(tmp,"      B %.3lf %.3lf %.3lf\n",qm->box[1][0]*10,qm->box[1][1]*10,qm->box[1][2]*10); 
            fprintf(tmp,"      C %.3lf %.3lf %.3lf\n",qm->box[2][0]*10,qm->box[2][1]*10,qm->box[2][2]*10);

             fgets(buf, 300, cp2kinp);
             continue;
        }

        //if the next line contains QM cell dimensions then update them
/*        if (strncmp(buf,"#QMCELL",7) == 0)
        {
            fputs(buf, tmp);
            //ABC
            fgets(buf, 300, cp2kinp);
            //Wrtie new cell dimensions instead
            fprintf(stderr,"Dimensions of QM subsystem cell has been changed to %.3lf x %.3lf x %.3lf nm\n", dimx, dimy, dimz);
            //Label for further modifications
            fprintf(tmp,"      ABC %.3lf %.3lf %.3lf\n", dimx*10, dimy*10, dimz*10);

            fgets(buf, 300, cp2kinp);
            continue;
        }*/

        fputs(buf, tmp);

        //read the next line from cp2k input
        fgets(buf, 300, cp2kinp);
    }

    fclose(tmp);
    fclose(cp2kinp);

    //Swap old input with a new one
    if (qm->cp2k_inp)
    {
        remove(qm->cp2k_inp);
        rename("cp2k.tmp.inp",qm->cp2k_inp);
    }
    else
    {
        remove("cp2k.inp");
        rename("cp2k.tmp.inp","cp2k.inp");
    }

}


void update_cp2k_pdb(t_QMrec *qm, t_MMrec *mm)
{
    char
          periodic_system[37][4]={"X  ","H  ","He ","Li ","Be ","B  ","C  ","N  ",
			"O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ",
			"S  ","Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ",
			"Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ",
			"Ge ","As ","Se ","Br ","Kr "};
    FILE *cp2kinp = NULL;
    int  i,j,k;

    //write PDB file
    cp2kinp = fopen("cp2k.pdb", "w");

    //write QM atoms
    for(i=0;i<qm->nrQMatoms;i++) {
       fprintf(cp2kinp,"ATOM ");
       fprintf(cp2kinp,"%6d ",i+1);
       fprintf(cp2kinp," %3s ",periodic_system[qm->atomicnumberQM[i]]);
       fprintf(cp2kinp," QM");
       fprintf(cp2kinp,"%6d     ",1);
       fprintf(cp2kinp,"%7.3lf %7.3lf %7.3lf  1.00  0.00         ",qm->xQM[i][0]*10,qm->xQM[i][1]*10,qm->xQM[i][2]*10);
       fprintf(cp2kinp," %3s ",periodic_system[qm->atomicnumberQM[i]]);
       fprintf(cp2kinp,"%.2lf\n",qm->atomicnumberQM[i]/1.0);
    }

    //write MM atoms
    for(i=0;i<mm->nrMMatoms;i++) {
       fprintf(cp2kinp,"ATOM ");
       fprintf(cp2kinp,"%6d ",qm->nrQMatoms+i+1);
       fprintf(cp2kinp," %3s ",periodic_system[0]);
       fprintf(cp2kinp," MM");
       fprintf(cp2kinp,"%6d     ",2);
       fprintf(cp2kinp,"%7.3lf %7.3lf %7.3lf  1.00  0.00         ",mm->xMM[i][0]*10,mm->xMM[i][1]*10,mm->xMM[i][2]*10);
       fprintf(cp2kinp," %3s ",periodic_system[0]);
       fprintf(cp2kinp,"%lf\n",mm->MMcharges[i]);
    }
    fprintf(cp2kinp,"END\n");

          
    fclose(cp2kinp);


}


void init_gaussian(t_QMrec *qm, t_MMrec *mm)
{
    FILE *out = NULL;
    ivec
          basissets[eQMbasisNR] = {{0, 3, 0},
                                 {0, 3, 0}, /*added for double sto-3g entry in names.c*/
                                 {5, 0, 0},
                                 {5, 0, 1},
                                 {5, 0, 11},
                                 {5, 6, 0},
                                 {1, 6, 0},
                                 {1, 6, 1},
                                 {1, 6, 11},
                                 {4, 6, 0}};
    char
         *buf = NULL;
    int
          i,j,k;


    //Check if cp2k input file exists
    buf = getenv("CP2K_INPUT");

    if (buf)
    {
        qm->cp2k_inp = buf;
        fprintf(stderr, "CP2K input = %s\n",qm->cp2k_inp);

    }


    //attempt to init CP2K
    cp2k_init();

    fprintf(stderr, "CP2K initialised...\n");
}





void write_gaussian_SH_input(int step, gmx_bool swap,
                             t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
    int
        i;
    gmx_bool
        bSA;
    FILE
       *out;
    t_QMMMrec
       *QMMMrec;
    QMMMrec = fr->qr;
    bSA     = (qm->SAstep > 0);

    out = fopen("input.com", "w");
    /* write the route */
    fprintf(out, "%s", "%scr=input\n");
    fprintf(out, "%s", "%rwf=input\n");
    fprintf(out, "%s", "%int=input\n");
    fprintf(out, "%s", "%d2e=input\n");
/*  if(step)
 *   fprintf(out,"%s","%nosave\n");
 */
    fprintf(out, "%s", "%chk=input\n");
    fprintf(out, "%s%d\n", "%mem=", qm->QMmem);
    fprintf(out, "%s%3d\n", "%nprocshare=", qm->nQMcpus);

    /* use the versions of
     * l301 that computes the interaction between MM and QM atoms.
     * l510 that can punch the CI coefficients
     * l701 that can do gradients on MM atoms
     */

    /* local version */
    fprintf(out, "%s%s%s",
            "%subst l510 ",
            qm->devel_dir,
            "/l510\n");
    fprintf(out, "%s%s%s",
            "%subst l301 ",
            qm->devel_dir,
            "/l301\n");
    fprintf(out, "%s%s%s",
            "%subst l701 ",
            qm->devel_dir,
            "/l701\n");

    fprintf(out, "%s%s%s",
            "%subst l1003 ",
            qm->devel_dir,
            "/l1003\n");
    fprintf(out, "%s%s%s",
            "%subst l9999 ",
            qm->devel_dir,
            "/l9999\n");
    /* print the nonstandard route
     */
    fprintf(out, "%s",
            "#P nonstd\n 1/18=10,20=1,38=1/1;\n");
    fprintf(out, "%s",
            " 2/9=110,15=1,17=6,18=5,40=1/2;\n");
    if (mm->nrMMatoms)
    {
        fprintf(out,
                " 3/5=%d,6=%d,7=%d,25=1,32=1,43=1,94=-2/1,2,3;\n",
                qm->SHbasis[0],
                qm->SHbasis[1],
                qm->SHbasis[2]); /*basisset stuff */
    }
    else
    {
        fprintf(out,
                " 3/5=%d,6=%d,7=%d,25=1,32=1,43=0,94=-2/1,2,3;\n",
                qm->SHbasis[0],
                qm->SHbasis[1],
                qm->SHbasis[2]); /*basisset stuff */
    }
    /* development */
    if (step+1) /* fetch initial guess from check point file */
    {           /* hack, to alyays read from chk file!!!!! */
        fprintf(out, "%s%d,%s%d%s", " 4/5=1,7=6,17=",
                qm->CASelectrons,
                "18=", qm->CASorbitals, "/1,5;\n");
    }
    else /* generate the first checkpoint file */
    {
        fprintf(out, "%s%d,%s%d%s", " 4/5=0,7=6,17=",
                qm->CASelectrons,
                "18=", qm->CASorbitals, "/1,5;\n");
    }
    /* the rest of the input depends on where the system is on the PES
     */
    if (swap && bSA)             /* make a slide to the other surface */
    {
        if (qm->CASorbitals > 6) /* use direct and no full diag */
        {
            fprintf(out, " 5/5=2,16=-2,17=10000000,28=2,32=2,38=6,97=100/10;\n");
        }
        else
        {
            if (qm->cpmcscf)
            {
                fprintf(out, " 5/5=2,6=%d,17=31000200,28=2,32=2,38=6,97=100/10;\n",
                        qm->accuracy);
                if (mm->nrMMatoms > 0)
                {
                    fprintf(out, " 7/7=1,16=-2,30=1/1;\n");
                }
                fprintf(out, " 11/31=1,42=1,45=1/1;\n");
                fprintf(out, " 10/6=1,10=700006,28=2,29=1,31=1,97=100/3;\n");
                fprintf(out, " 7/30=1/16;\n 99/10=4/99;\n");
            }
            else
            {
                fprintf(out, " 5/5=2,6=%d,17=11000000,28=2,32=2,38=6,97=100/10;\n",
                        qm->accuracy);
                fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
            }
        }
    }
    else if (bSA)                /* do a "state-averaged" CAS calculation */
    {
        if (qm->CASorbitals > 6) /* no full diag */
        {
            fprintf(out, " 5/5=2,16=-2,17=10000000,28=2,32=2,38=6/10;\n");
        }
        else
        {
            if (qm->cpmcscf)
            {
                fprintf(out, " 5/5=2,6=%d,17=31000200,28=2,32=2,38=6/10;\n",
                        qm->accuracy);
                if (mm->nrMMatoms > 0)
                {
                    fprintf(out, " 7/7=1,16=-2,30=1/1;\n");
                }
                fprintf(out, " 11/31=1,42=1,45=1/1;\n");
                fprintf(out, " 10/6=1,10=700006,28=2,29=1,31=1/3;\n");
                fprintf(out, " 7/30=1/16;\n 99/10=4/99;\n");
            }
            else
            {
                fprintf(out, " 5/5=2,6=%d,17=11000000,28=2,32=2,38=6/10;\n",
                        qm->accuracy);
                fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
            }
        }
    }
    else if (swap) /* do a "swapped" CAS calculation */
    {
        if (qm->CASorbitals > 6)
        {
            fprintf(out, " 5/5=2,16=-2,17=0,28=2,32=2,38=6,97=100/10;\n");
        }
        else
        {
            fprintf(out, " 5/5=2,6=%d,17=1000000,28=2,32=2,38=6,97=100/10;\n",
                    qm->accuracy);
        }
        fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
    }
    else /* do a "normal" CAS calculation */
    {
        if (qm->CASorbitals > 6)
        {
            fprintf(out, " 5/5=2,16=-2,17=0,28=2,32=2,38=6/10;\n");
        }
        else
        {
            fprintf(out, " 5/5=2,6=%d,17=1000000,28=2,32=2,38=6/10;\n",
                    qm->accuracy);
        }
        fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
    }
    fprintf(out, "\ninput-file generated by gromacs\n\n");
    fprintf(out, "%2d%2d\n", qm->QMcharge, qm->multiplicity);
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n",
                qm->atomicnumberQM[i],
                qm->xQM[i][XX]/BOHR2NM,
                qm->xQM[i][YY]/BOHR2NM,
                qm->xQM[i][ZZ]/BOHR2NM);
    }
    /* MM point charge data */
    if (QMMMrec->QMMMscheme != eQMMMschemeoniom && mm->nrMMatoms)
    {
        fprintf(out, "\n");
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            fprintf(out, "%10.7f  %10.7f  %10.7f %8.4f\n",
                    mm->xMM[i][XX]/BOHR2NM,
                    mm->xMM[i][YY]/BOHR2NM,
                    mm->xMM[i][ZZ]/BOHR2NM,
                    mm->MMcharges[i]);
        }
    }
    if (bSA) /* put the SA coefficients at the end of the file */
    {
        fprintf(out, "\n%10.8f %10.8f\n",
                qm->SAstep*0.5/qm->SAsteps,
                1-qm->SAstep*0.5/qm->SAsteps);
        fprintf(stderr, "State Averaging level = %d/%d\n", qm->SAstep, qm->SAsteps);
    }
    fprintf(out, "\n");
    fclose(out);
}  /* write_gaussian_SH_input */

void write_gaussian_input(int step, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
    int
        i;
    t_QMMMrec
       *QMMMrec;
    FILE
       *out;

    QMMMrec = fr->qr;
    out     = fopen("input.com", "w");
    /* write the route */

    if (qm->QMmethod >= eQMmethodRHF)
    {
        fprintf(out, "%s",
                "%chk=input\n");
    }
    else
    {
        fprintf(out, "%s",
                "%chk=se\n");
    }
    if (qm->nQMcpus > 1)
    {
        fprintf(out, "%s%3d\n",
                "%nprocshare=", qm->nQMcpus);
    }
    fprintf(out, "%s%d\n",
            "%mem=", qm->QMmem);
    fprintf(out, "%s%s%s",
            "%subst l701 ", qm->devel_dir, "/l701\n");
    fprintf(out, "%s%s%s",
            "%subst l301 ", qm->devel_dir, "/l301\n");
    fprintf(out, "%s%s%s",
            "%subst l9999 ", qm->devel_dir, "/l9999\n");
    if (step)
    {
        fprintf(out, "%s",
                "#T ");
    }
    else
    {
        fprintf(out, "%s",
                "#P ");
    }
    if (qm->QMmethod == eQMmethodB3LYPLAN)
    {
        fprintf(out, " %s",
                "B3LYP/GEN Pseudo=Read");
    }
    else
    {
        fprintf(out, " %s",
                eQMmethod_names[qm->QMmethod]);

        if (qm->QMmethod >= eQMmethodRHF)
        {
            if (qm->QMmethod == eQMmethodCASSCF)
            {
                /* in case of cas, how many electrons and orbitals do we need?
                 */
                fprintf(out, "(%d,%d)",
                        qm->CASelectrons, qm->CASorbitals);
            }
            fprintf(out, "/%s",
                    eQMbasis_names[qm->QMbasis]);
        }
    }
    if (QMMMrec->QMMMscheme == eQMMMschemenormal && mm->nrMMatoms)
    {
        fprintf(out, " %s",
                "Charge ");
    }
    if (step || qm->QMmethod == eQMmethodCASSCF)
    {
        /* fetch guess from checkpoint file, always for CASSCF */
        fprintf(out, "%s", " guess=read");
    }
    fprintf(out, "\nNosymm units=bohr\n");

    fprintf(out, "FORCE Punch=(Derivatives) ");
    fprintf(out, "iop(3/33=1)\n\n");
    fprintf(out, "input-file generated by gromacs\n\n");
    fprintf(out, "%2d%2d\n", qm->QMcharge, qm->multiplicity);
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n",
                qm->atomicnumberQM[i],
                qm->xQM[i][XX]/BOHR2NM,
                qm->xQM[i][YY]/BOHR2NM,
                qm->xQM[i][ZZ]/BOHR2NM);
    }

    /* Pseudo Potential and ECP are included here if selected (MEthod suffix LAN) */
    if (qm->QMmethod == eQMmethodB3LYPLAN)
    {
        fprintf(out, "\n");
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            if (qm->atomicnumberQM[i] < 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n****\n", eQMbasis_names[qm->QMbasis]);

        for (i = 0; i < qm->nrQMatoms; i++)
        {
            if (qm->atomicnumberQM[i] > 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n****\n\n", "lanl2dz");

        for (i = 0; i < qm->nrQMatoms; i++)
        {
            if (qm->atomicnumberQM[i] > 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n", "lanl2dz");
    }



    /* MM point charge data */
    if (QMMMrec->QMMMscheme != eQMMMschemeoniom && mm->nrMMatoms)
    {
        fprintf(stderr, "nr mm atoms in gaussian.c = %d\n", mm->nrMMatoms);
        fprintf(out, "\n");
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            fprintf(out, "%10.7f  %10.7f  %10.7f %8.4f\n",
                    mm->xMM[i][XX]/BOHR2NM,
                    mm->xMM[i][YY]/BOHR2NM,
                    mm->xMM[i][ZZ]/BOHR2NM,
                    mm->MMcharges[i]);
        }
    }
    fprintf(out, "\n");


    fclose(out);

}  /* write_gaussian_input */

real read_gaussian_output(rvec QMgrad[], rvec MMgrad[], t_QMrec *qm, t_MMrec *mm)
{
    int
        i, j, atnum;
    char
        buf[300];
    real
        QMener;
    FILE
       *in;

    in = fopen("fort.7", "r");

    /* (There was additional content in the file in case
     *    of QM optimizations / transition state search,
     *    which was removed.
     */
    /* the next line is the energy and in the case of CAS, the energy
     * difference between the two states.
     */
    if (NULL == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }

#if GMX_DOUBLE
    sscanf(buf, "%lf\n", &QMener);
#else
    sscanf(buf, "%f\n", &QMener);
#endif
    /* next lines contain the gradients of the QM atoms */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        if (NULL == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#endif
    }
    /* the next lines are the gradients of the MM atoms */
    if (qm->QMmethod >= eQMmethodRHF)
    {
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            if (NULL == fgets(buf, 300, in))
            {
                gmx_fatal(FARGS, "Error reading Gaussian output");
            }
#if GMX_DOUBLE
            sscanf(buf, "%lf %lf %lf\n",
                   &MMgrad[i][XX],
                   &MMgrad[i][YY],
                   &MMgrad[i][ZZ]);
#else
            sscanf(buf, "%f %f %f\n",
                   &MMgrad[i][XX],
                   &MMgrad[i][YY],
                   &MMgrad[i][ZZ]);
#endif
        }
    }
    fclose(in);
    return(QMener);
}

real read_gaussian_SH_output(rvec QMgrad[], rvec MMgrad[], int step, t_QMrec *qm, t_MMrec *mm)
{
    int
        i;
    char
        buf[300];
    real
        QMener, DeltaE;
    FILE
       *in;

    in = fopen("fort.7", "r");
    /* first line is the energy and in the case of CAS, the energy
     * difference between the two states.
     */
    if (NULL == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }

#if GMX_DOUBLE
    sscanf(buf, "%lf %lf\n", &QMener, &DeltaE);
#else
    sscanf(buf, "%f %f\n",  &QMener, &DeltaE);
#endif

    /* switch on/off the State Averaging */

    if (DeltaE > qm->SAoff)
    {
        if (qm->SAstep > 0)
        {
            qm->SAstep--;
        }
    }
    else if (DeltaE < qm->SAon || (qm->SAstep > 0))
    {
        if (qm->SAstep < qm->SAsteps)
        {
            qm->SAstep++;
        }
    }

    /* for debugging: */
    fprintf(stderr, "Gap = %5f,SA = %3d\n", DeltaE, (qm->SAstep > 0));
    /* next lines contain the gradients of the QM atoms */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        if (NULL == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }

#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#endif
    }
    /* the next lines, are the gradients of the MM atoms */

    for (i = 0; i < mm->nrMMatoms; i++)
    {
        if (NULL == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &MMgrad[i][XX],
               &MMgrad[i][YY],
               &MMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &MMgrad[i][XX],
               &MMgrad[i][YY],
               &MMgrad[i][ZZ]);
#endif
    }

    /* the next line contains the two CI eigenvector elements */
    if (NULL == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }
    if (!step)
    {
        sscanf(buf, "%d", &qm->CIdim);
        snew(qm->CIvec1, qm->CIdim);
        snew(qm->CIvec1old, qm->CIdim);
        snew(qm->CIvec2, qm->CIdim);
        snew(qm->CIvec2old, qm->CIdim);
    }
    else
    {
        /* before reading in the new current CI vectors, copy the current
         * CI vector into the old one.
         */
        for (i = 0; i < qm->CIdim; i++)
        {
            qm->CIvec1old[i] = qm->CIvec1[i];
            qm->CIvec2old[i] = qm->CIvec2[i];
        }
    }
    /* first vector */
    for (i = 0; i < qm->CIdim; i++)
    {
        if (NULL == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf\n", &qm->CIvec1[i]);
#else
        sscanf(buf, "%f\n", &qm->CIvec1[i]);
#endif
    }
    /* second vector */
    for (i = 0; i < qm->CIdim; i++)
    {
        if (NULL == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf\n", &qm->CIvec2[i]);
#else
        sscanf(buf, "%f\n", &qm->CIvec2[i]);
#endif
    }
    fclose(in);
    return(QMener);
}

real inproduct(real *a, real *b, int n)
{
    int
        i;
    real
        dot = 0.0;

    /* computes the inner product between two vectors (a.b), both of
     * which have length n.
     */
    for (i = 0; i < n; i++)
    {
        dot += a[i]*b[i];
    }
    return(dot);
}

int hop(int step, t_QMrec *qm)
{
    int
        swap = 0;
    real
        d11 = 0.0, d12 = 0.0, d21 = 0.0, d22 = 0.0;

    /* calculates the inproduct between the current Ci vector and the
     * previous CI vector. A diabatic hop will be made if d12 and d21
     * are much bigger than d11 and d22. In that case hop returns true,
     * otherwise it returns false.
     */
    if (step) /* only go on if more than one step has been done */
    {
        d11 = inproduct(qm->CIvec1, qm->CIvec1old, qm->CIdim);
        d12 = inproduct(qm->CIvec1, qm->CIvec2old, qm->CIdim);
        d21 = inproduct(qm->CIvec2, qm->CIvec1old, qm->CIdim);
        d22 = inproduct(qm->CIvec2, qm->CIvec2old, qm->CIdim);
    }
    fprintf(stderr, "-------------------\n");
    fprintf(stderr, "d11 = %13.8f\n", d11);
    fprintf(stderr, "d12 = %13.8f\n", d12);
    fprintf(stderr, "d21 = %13.8f\n", d21);
    fprintf(stderr, "d22 = %13.8f\n", d22);
    fprintf(stderr, "-------------------\n");

    if ((fabs(d12) > 0.5) && (fabs(d21) > 0.5))
    {
        swap = 1;
    }

    return(swap);
}

void do_gaussian(int step, char *exe)
{
    char
        buf[STRLEN];

    /* make the call to the gaussian binary through system()
     * The location of the binary will be picked up from the
     * environment using getenv().
     */
    if (step) /* hack to prevent long inputfiles */
    {
        sprintf(buf, "%s < %s > %s",
                exe,
                "input.com",
                "input.log");
    }
    else
    {
        sprintf(buf, "%s < %s > %s",
                exe,
                "input.com",
                "input.log");
    }
    fprintf(stderr, "Calling '%s'\n", buf);
    if (system(buf) != 0)
    {
        gmx_fatal(FARGS, "Call to '%s' failed\n", buf);
    }
}

real call_gaussian(t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    /* normal gaussian jobs */
    static int
        step = 0;
    int
        i, j;
    real
        QMener = 0.0;
    rvec
       *QMgrad, *MMgrad;
    double
       *Crd;
    char
       *exe;
    gmx_bool
        samebox;

    snew(QMgrad, qm->nrQMatoms);
    snew(MMgrad, mm->nrMMatoms);

    if (!step) {
       //Initialize CP2K Force_env on the first step
       if (qm->cp2k_inp)
       {
            cp2k_create_force_env(&qm->force_env, qm->cp2k_inp, "cp2k.out");
       }
       else
       {
            write_cp2k_input_pdb(qm, mm);
            cp2k_create_force_env(&qm->force_env, "cp2k.inp", "cp2k.out");
       }

       fprintf(stderr,"Number of QM atoms= %d\n",qm->nrQMatoms);
       fprintf(stderr,"Number of MM atoms= %d\n",mm->nrMMatoms);
    }

    //Update PDB file
    update_cp2k_pdb(qm, mm);

    //check if old box is the same as new
    samebox = TRUE;
    for (i = 0; i<DIM; i++)
       for (j = 0; j<DIM; j++)
          if (fabs(qm->oldbox[i][j] - qm->box[i][j]) > GMX_REAL_EPS)
          {
             samebox = FALSE;
             break;
          }

    //if box is different then we should reinitialize CP2K Force environment
    if (!samebox)
    {
        fprintf(stderr,"Box changed - reinitializing CP2K\n");
        //update input
        update_cp2k_input(qm, mm);
        //destroy old force environment
        cp2k_destroy_force_env(qm->force_env);
        //Create new force environment
        if (qm->cp2k_inp)
        {
            cp2k_create_force_env(&qm->force_env, qm->cp2k_inp, "cp2k.out");
        }
        else
        {
            cp2k_create_force_env(&qm->force_env, "cp2k.inp", "cp2k.out");
        }
    }
        
    //Gather all coordinates in one array

    snew(Crd, 3*(qm->nrQMatoms + mm->nrMMatoms));
    for(i=0;i<qm->nrQMatoms;i++) {
       Crd[3*i]  =qm->xQM[i][0]/BOHR2NM;
       Crd[3*i+1]=qm->xQM[i][1]/BOHR2NM;
       Crd[3*i+2]=qm->xQM[i][2]/BOHR2NM;
    }
    for(i=0;i<mm->nrMMatoms;i++) {
       Crd[3*(qm->nrQMatoms+i)]  =mm->xMM[i][0]/BOHR2NM;
       Crd[3*(qm->nrQMatoms+i)+1]=mm->xMM[i][1]/BOHR2NM;
       Crd[3*(qm->nrQMatoms+i)+2]=mm->xMM[i][2]/BOHR2NM;
    }


    //Set new positions
    cp2k_set_positions(qm->force_env, Crd, 3*(qm->nrQMatoms + mm->nrMMatoms));
    fprintf(stderr,"cp2k_set_positions - DONE\n");

    //Run CP2K
    cp2k_calc_energy_force(qm->force_env);
    fprintf(stderr,"cp2k_calc_energy_force - DONE\n");
    
    //Get QM + QM-MM Energy
    cp2k_get_potential_energy(qm->force_env,&QMener);
    fprintf(stderr,"QMener=%.12lf\n",QMener);

    //Get Gradient
    cp2k_get_forces(qm->force_env,Crd,3*(qm->nrQMatoms + mm->nrMMatoms));
    //for (i = 0; i < qm->nrQMatoms + mm->nrMMatoms; i++)
    //   fprintf(stderr,"Grad[%d]=%.12lf %.12lf %.12lf\n",i,Crd[3*i],Crd[3*i+1],Crd[3*i+2]);

    for (i = 0; i < qm->nrQMatoms; i++)
    {
       QMgrad[i][0]=-Crd[3*i];
       QMgrad[i][1]=-Crd[3*i+1];
       QMgrad[i][2]=-Crd[3*i+2];

    }
    for(i=0;i<mm->nrMMatoms;i++) {
       MMgrad[i][0]=-Crd[3*(qm->nrQMatoms+i)];  
       MMgrad[i][1]=-Crd[3*(qm->nrQMatoms+i)+1];
       MMgrad[i][2]=-Crd[3*(qm->nrQMatoms+i)+2];
    }

    sfree(Crd);

    /* put the QMMM forces in the force array and to the fshift
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (i = 0; i < mm->nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
            fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
        }
    }
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    step++;
    return(QMener);

} /* call_gaussian */

real call_gaussian_SH(t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
    /* a gaussian call routine intended for doing diabatic surface
     * "sliding". See the manual for the theoretical background of this
     * TSH method.
     */
    static int
        step = 0;
    int
        state, i, j;
    real
        QMener = 0.0;
    static  gmx_bool
        swapped = FALSE; /* handle for identifying the current PES */
    gmx_bool
        swap = FALSE;    /* the actual swap */
    rvec
       *QMgrad, *MMgrad;
    char
       *buf;
    char
       *exe;

    snew(exe, 30);
    sprintf(exe, "%s/%s", qm->gauss_dir, qm->gauss_exe);
    /* hack to do ground state simulations */
    if (!step)
    {
        snew(buf, 20);
        buf = getenv("GMX_QM_GROUND_STATE");
        if (buf)
        {
            sscanf(buf, "%d", &state);
        }
        else
        {
            state = 2;
        }
        if (state == 1)
        {
            swapped = TRUE;
        }
    }
    /* end of hack */


    /* copy the QMMMrec pointer */
    snew(QMgrad, qm->nrQMatoms);
    snew(MMgrad, mm->nrMMatoms);
    /* at step 0 there should be no SA */
    /*  if(!step)
     * qr->bSA=FALSE;*/
    /* temporray set to step + 1, since there is a chk start */
    write_gaussian_SH_input(step, swapped, fr, qm, mm);

    do_gaussian(step, exe);
    QMener = read_gaussian_SH_output(QMgrad, MMgrad, step, qm, mm);

    /* check for a surface hop. Only possible if we were already state
     * averaging.
     */
    if (qm->SAstep > 0)
    {
        if (!swapped)
        {
            swap    = (step && hop(step, qm));
            swapped = swap;
        }
        else /* already on the other surface, so check if we go back */
        {
            swap    = (step && hop(step, qm));
            swapped = !swap; /* so swapped shoud be false again */
        }
        if (swap)            /* change surface, so do another call */
        {
            write_gaussian_SH_input(step, swapped, fr, qm, mm);
            do_gaussian(step, exe);
            QMener = read_gaussian_SH_output(QMgrad, MMgrad, step, qm, mm);
        }
    }
    /* add the QMMM forces to the gmx force array and fshift
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (i = 0; i < mm->nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
            fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
        }
    }
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    fprintf(stderr, "step %5d, SA = %5d, swap = %5d\n",
            step, (qm->SAstep > 0), swapped);
    step++;
    free(exe);
    return(QMener);

} /* call_gaussian_SH */

/* end of gaussian sub routines */

#else
int
    gmx_qmmm_gaussian_empty;
#endif
