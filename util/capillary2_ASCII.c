/*****************************************************************************
 *
 *  capillary.c
 *
 *  Compile and link with -lm for the maths library. 
 *
 *  This utility produces an output file suitable for initialising
 *  a capillary structure in Ludwig. 
 *
 *  It is assumed the periodic dimension is z, and that the system
 *  is square in the x-y directions. No 'leakage' in the x-y
 *  directions is ensured by making the last site in each direction
 *  solid.
 *
 *  The various system parameters should be set at comiple time,
 *  and are described below. The output file is always in BINARY
 *  format.
 *
 *  1. Output capaillary structure
 *  Set the required parameters and invoke with no argument
 *      ./a.out
 *   
 *  2. Profiles
 *  If the program is invoked with a single phi output file
 *  argument, e.g.,
 *      ./a.out phi-001000.001-001
 *  a scatter plot of the profile of the interface in the
 *  wetting half of the capillary will be produced. That
 *  is height vs r, the radial distance from the centre.
 *  The output file should match the capillary structure!
 *
 *  Edinburgh Soft Matter and Statistcal Physics Group and
 *  Edinburgh Parallel Computing Centre
 *  (c) 2008-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

/* This is a copy from ../src/map.h; it would be better to include
 * directly, but that incurs additional dependencies on targetDP.h */

enum map_status {MAP_FLUID, MAP_BOUNDARY, MAP_COLLOID, MAP_STATUS_MAX};


//const double rho0=0.1;
//const double kon=1e-1;
//const double koff=1e-2;
const double fcr = 1 ; // fraction of CR sites in the x-axis along the channel

const int crystalline_cell_size = 10; //added for implementation of crystalline capillaries

/* CROSS SECTION */
/* You can choose a square or circular cross section */

enum {CIRCLE, SQUARE, XWALL, YWALL, ZWALL, XWALL_OBSTACLES, XWALL_BOTTOM,
      SPECIAL_CROSS, SCC, BCC, FCC,RESERVOIR}; //changed for implementation of crystalline capillaries
const int xsection = RESERVOIR;

/*Modify the local geometry of the wall*/

int obstacle_number = 1; /* number of obstacles per wall */
int obstacle_length = 6; /* along the wall direction */
int obstacle_height = 10; /* perpendicular from wall */
int obstacle_depth  = 6; /* perpendicular to length and height */
			 /* NOTE: obstacle_depth == xmax/ymax/zmax 
				  means obstacles don't have a z-boundary */

/* BINARY FREE ENERGY PARAMETERS */
/* Set the fluid and solid free energy parameters. The fluid parameters
 * must match those used in the main calculation. See Desplat et al.
 * Comp. Phys. Comm. (2001) for details. */

typedef struct fe_symm_param_s fe_symm_param_t;

struct fe_symm_param_s {
  double a;      /* Free energy parameter A < 0 */
  double b;      /* Free energy parameter B > 0 */
  double kappa;  /* Free energy parameter */
  double c;      /* Surface free energy C */
  double h;      /* Surface free energy H */
};

/* TERNARY FREE ENERGY PARAMETERS */
/* Set the fluid and solid free energy parameters. The fluid parameters
* must match those used in the main calculation. See Semprebon et al.
* Phys. Rev. E (2016) for details. */

const double kappa1 = 0.012;
const double kappa2 = 0.05;
const double kappa3 = 0.05;
const double TERNARY_H1 = 0.000175;
const double TERNARY_H2 = -0.00175;



const double alpha = 1.000;


/* WETTING */
/* A section of capillary between z1 and z2 (inclusive) will have
 * wetting property H = H, the remainder H = 0 */

const int z1_h = 1;
const int z2_h = 36;

/* OUTPUT */
/* You can generate a file with solid/fluid status information only,
 * or one which includes the wetting parameter H or charge Q. */

/* Wetting: Please don't use STATUS_WITH_H; use STATUS_WITH_C_H
   and set C = 0, if you just want H */

enum {STATUS_ONLY, STATUS_WITH_H, STATUS_WITH_C_H, STATUS_WITH_SIGMA,STATUS_WITH_H1_H2};
const int output_type = STATUS_ONLY;

/* OUTPUT FILENAME */

const char * outputfilename = "map-000000000.001-001";

//int map_special_cross(char * map, int * nsolid);

/*****************************************************************************
 *
 *  main program
 *
 *****************************************************************************/

int main(int argc, char ** argv) {

  // PARSE INPUT
  /* SYSTEM SIZE */
/* Set the system size as desired. Clearly, this must match the system
 * set in the main input file for Ludwig. */
  int xmax, ymax, zmax; // system size
  /* SURFACE CHARGE */
  double sigma, fwall, fwallplus, fwallminus, sigmaplus, sigmaminus, w, l;
  zmax = 2;
  xmax=80;
  ymax=20;
  sigma=0.05;
  sigmaplus=0.05;
  sigmaminus=-0.05;
  w=16;
  fwall=0.5;
    int opt;/*
    while ((opt = getopt(argc, argv, "l:w:y:p:m:f:")) != -1) {
        switch (opt) {
        case 'l': xmax=atoi(optarg);
        case 'w': w=atoi(optarg);
        case 'y': ymax=atoi(optarg);
        case 'p': sigmaplus=atof(optarg);
        case 'm': sigmaminus=atof(optarg);
        case 'f': fwall=atof(optarg);
     //   default:
     //       fprintf(stderr, "Usage: %s [-l (length) -w (width) -p (sigma plus) -m (sigma minus) -f (fwall plus/minus)] \n", argv[0]);
     //       exit(EXIT_FAILURE);
        }
    }
*/
// END PARSE INPUT
   fwallminus=fwall/2; // both sections are symetric. If a uniformly charged channel is desired, we can set sigmaminus=sigmaplus and f=0.5.
   fwallplus=fwall/2;
   if (zmax < 0 || ymax < 0 || zmax < 0 || fwall < 0 || w < 0 || w > ymax-2)
   {fprintf(stderr, "ERROR: All parameters must be set. Usage: %s [-l (length) -w (width) -p (sigma plus) -m (sigma minus) -f (fwall plus/minus)] \n", argv[0]);
   exit(1);
   }

  char * map_in;
  FILE * fp_orig;
  int i, j, k, n;
  int nsolid = 0;
  int k_pic;

  double * map_h;   /* for wetting coefficient H */
  double * map_c;   /* for additional wetting coefficient C */

  double * map_sig; /* for (surface) charge */

  double rc = 0.5*(xmax-2);
  double x0 = 0.5*xmax + 0.5;
  double y0 = 0.5*ymax + 0.5;
  double x, y, r;
  double h, h1, theta;
  double f1,f2,f3,cos_theta12,cos_theta23,cos_theta31;
  double theta12,theta23,theta31;
  int L0,b0,b1,length,channel_start,channel_end,width;
  //added for implementation of crystalline capillaries
  double crystalline_cell_radius;
  double diff_x_edges, diff_y_edges, diff_z_edges, r_edges;

  int iobst;
  int obst_start[2*obstacle_number][3];
  int obst_stop[2*obstacle_number][3];
  int gap_length;

  FILE  * WriteFile;
  char  file[800];

  /* Some default values */
  fe_symm_param_t fe = {.a =    -0.0625,
                        .b =     0.0625,
                        .kappa = 0.053,
                        .c     = 0.0,
                        .h     = 0.0};


  if (output_type == STATUS_WITH_H || output_type == STATUS_WITH_C_H) {

    printf("Free energy parameters:\n");
    printf("free energy parameter kappa = %f\n", fe.kappa);
    printf("free energy parameter B     = %f\n", fe.b);
    printf("surface free energy   H     = %f\n", fe.h);
    h = fe.h*sqrt(1.0/(fe.kappa*fe.b));
    printf("dimensionless parameter h   = %f\n", h);
    h1 = 0.5*(-pow(1.0 - h, 1.5) + pow(1.0 + h, 1.5));
    printf("dimensionless parameter h1=cos(theta)   = %f\n", h1);
    theta = acos(h1);
    printf("contact angle theta         = %f radians\n", theta);
    theta = theta*180.0/(4.0*atan(1.0));
    printf("                            = %f degrees\n", theta);

  }

  if (output_type == STATUS_WITH_SIGMA) {
    printf("Surface charge sigma = %f\n", sigma);
  }
  if (output_type == STATUS_WITH_H1_H2)
  {
      printf("Ternary free energy parameters:\n");
      printf("free energy parameter kappa1 = %f\n", kappa1);
      printf("free energy parameter kappa2 = %f\n", kappa2);
      printf("free energy parameter kappa3 = %f\n", kappa3);
      printf("free energy parameter alpha = %f\n", alpha);
      f1=pow(alpha*kappa1+4*TERNARY_H1,1.5)-pow(alpha*kappa1-4*TERNARY_H1,1.5);
      f1=f1/sqrt(alpha*kappa1);
      f2=pow(alpha*kappa2+4*TERNARY_H2,1.5)-pow(alpha*kappa2-4*TERNARY_H2,1.5);
      f2=f2/sqrt(alpha*kappa2);
      
      {
	/* Constraint */
	double TERNARY_H3 = kappa3*(-(TERNARY_H1/kappa1)-(TERNARY_H2/kappa2));
	f3=pow(alpha*kappa3+4*TERNARY_H3,1.5)-pow(alpha*kappa3-4*TERNARY_H3,1.5);
      }
      f3=f3/sqrt(alpha*kappa3);
      cos_theta12=f1/(2.0*(kappa1+kappa2))-f2/(2.0*(kappa1+kappa2));
      cos_theta23=f2/(2.0*(kappa2+kappa3))-f3/(2.0*(kappa2+kappa3));
      cos_theta31=f3/(2.0*(kappa3+kappa1))-f1/(2.0*(kappa3+kappa1));
      printf("dimensionless parameters cos(theta12)   = %f\n", cos_theta12);
      printf("dimensionless parameters cos(theta23)   = %f\n", cos_theta23);
      printf("dimensionless parameters cos(theta13)   = %f\n", cos_theta31);
      theta12=acos(cos_theta12);
      theta23=acos(cos_theta23);
      theta31=acos(cos_theta31);
      printf("contact angle theta12         = %f radians\n", theta12);
      theta12=theta12*180.0/(4.0*atan(1.0));
      printf("contact angle theta12         = %f degrees\n", theta12);
      printf("contact angle theta23         = %f radians\n", theta23);
      theta23=theta23*180.0/(4.0*atan(1.0));
      printf("contact angle theta23         = %f degrees\n", theta23);
      printf("contact angle theta31         = %f radians\n", theta31);
      theta31=theta31*180.0/(4.0*atan(1.0));
      printf("contact angle theta31         = %f degrees\n", theta31);




  }

  map_in = (char *) malloc(xmax*ymax*zmax*sizeof(char));
  if (map_in == NULL) exit(-1);

  map_h = (double *) malloc(xmax*ymax*zmax*sizeof(double));
  if (map_h == NULL) exit(-1);

  map_c = (double *) malloc(xmax*ymax*zmax*sizeof(double));
  if (map_c == NULL) exit(-1);

  map_sig = (double *) malloc(xmax*ymax*zmax*sizeof(double));
  if (map_sig == NULL) exit(-1);

  k_pic = 1; /* Picture */

  /* Begin switch */
  switch (xsection) {

  case RESERVOIR:
// Remember that capillary.c has a different counting: node 0 will be node 1 in simulations
// Creates a reservoir with a trap inside.
// The channel starts at b0, and ends at b1. The trap is located at a distance set_trap from b0.
// The trap has a width of "width" solid nodes.
    L0=ymax-2;
    b0=0; // xmax*0.33333333333333;
    length=xmax*0.333333333333;
    b1=xmax; // b0+length;
    // width=(ymax-w) / 2;                    //wall width (from 0 to fluid node)
    channel_start=(ymax-w) / 2 - 1;
    channel_end=(ymax+w) / 2 ; // ymax-1-w;  //we have _width nodes of solid from 0 to width, and from ymax to ymax-width
    printf("init system %i %i %i \n",b0,b1,channel_end);
    for (i = 0; i < xmax; i++) {
        for (j = 0; j < ymax; j++) {
    for (k = 0; k < zmax; k++) {
        n = ymax*zmax*i + zmax*j + k;
        map_in[n] = MAP_FLUID;
        map_sig[n] = 0.0;
        if (output_type == STATUS_WITH_H) { map_h[n] = 0.0; }
            if (output_type == STATUS_WITH_C_H) { map_h[n] = 0.0; map_c[n] = 0.0; }
            if (output_type == STATUS_WITH_H1_H2) { map_h[n] = 0.0; map_c[n] = 0.0; }
            if ((i>=b0 && i<=b1) || (j==0) || (j==ymax-1) ){

            if ((j <= channel_start || j >= channel_end)   ){
              map_in[n] = MAP_BOUNDARY;
              if (output_type == STATUS_WITH_SIGMA) {
                //if ((j == channel_start || j == channel_end)   ){ only 1 charged layer
                if ((j <= channel_start || j >= channel_end)   ){

                if (i >= xmax*(0.25 - fwallplus / 2.0) && i < xmax*(0.25 + fwallplus / 2.0)) {
                  map_sig[n] = sigmaplus;
                  }
		            else if (i >= xmax*(0.75 - fwallminus / 2.0) && i < xmax*(0.75 + fwallminus / 2.0)) {
                  map_sig[n] = sigmaminus;
                  
		            }
		            }
		            else {map_sig[n] = 1e-13;} // small number
              }

              if (output_type == STATUS_WITH_H) { map_h[n] = fe.h; }
                  if (output_type == STATUS_WITH_C_H) {
                map_h[n] = fe.h;
                map_c[n] = fe.c;
              }
                  if (output_type == STATUS_WITH_H1_H2) {
                map_h[n] = TERNARY_H1;
                map_c[n] = TERNARY_H2;
              }
              ++nsolid;
            }
          }
            }
    }
    }
          break;

  case ZWALL:

    for (i = 0; i < xmax; i++) {
      for (j = 0; j < ymax; j++) {
	for (k = 0; k < zmax; k++) {
	  n = ymax*zmax*i + zmax*j + k;
	  map_in[n] = MAP_FLUID;
	  map_sig[n] = 0.0;
	  if (output_type == STATUS_WITH_H) { map_h[n] = 0.0; }
          if (output_type == STATUS_WITH_C_H) { map_h[n] = 0.0; map_c[n] = 0.0; }
          if (output_type == STATUS_WITH_H1_H2) { map_h[n] = 0.0; map_c[n] = 0.0; }
	  if (k == 0 || k == zmax - 1) {
	    map_in[n] = MAP_BOUNDARY;
	    if (output_type == STATUS_WITH_SIGMA) {
	      map_sig[n] = sigma;
	    }
	    if (output_type == STATUS_WITH_H) { map_h[n] = fe.h; }
            if (output_type == STATUS_WITH_C_H) {
	      map_h[n] = fe.h;
	      map_c[n] = fe.c;
	    }
	    if (output_type == STATUS_WITH_H1_H2) {
	      map_h[n] = TERNARY_H1;
	      map_c[n] = TERNARY_H2;
	    }
	    ++nsolid;
	  }
	}
      }
    }

    break;
  default:
    printf("No cross-section!\n");
    /* End switch */
  }

  /* picture */

  printf("\nCross section (%d = fluid, %d = solid)\n", MAP_FLUID, MAP_BOUNDARY);

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
	n = ymax*zmax*i + zmax*j + k_pic;
      
	if (map_in[n] == MAP_BOUNDARY) printf(" %d", MAP_BOUNDARY);
	if (map_in[n] == MAP_FLUID)    printf(" %d", MAP_FLUID);
    }
    printf("\n");
  }

//surface charge cross section
k_pic = 1;


    printf("\ncharge cross section (+ = positive, - = negative)\n");

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
	n = ymax*zmax*i + zmax*j + k_pic;
      
	if (map_sig[n] == sigmaplus) printf("+");
	else if (map_sig[n] == sigmaminus) printf("-");
  else  printf("0");

    }
    printf("\n");
  }


  if (output_type == STATUS_WITH_H)  {
    sprintf(file,"Configuration_capillary.dat");
    WriteFile=fopen(file,"w");
    fprintf(WriteFile,"#x y z n map H\n");

    for (i = 0; i < xmax; i++) {
      for (j = 0; j < ymax; j++) {
	for (k = 0; k < zmax; k++) {

	n = ymax*zmax*i + zmax*j + k;

	if (map_in[n] == MAP_BOUNDARY) { fprintf(WriteFile,"%i %i %i %i %d %f\n", i, j, k, n, MAP_BOUNDARY, map_h[n]); }
	if (map_in[n] == MAP_FLUID)    { fprintf(WriteFile,"%i %i %i %i %d %f\n", i, j, k, n, MAP_FLUID, map_h[n]); }

	}  
      }  
    }
    fclose(WriteFile);
  }


  if (output_type == STATUS_WITH_C_H)  {
    sprintf(file,"Configuration_capillary.dat");
    WriteFile=fopen(file,"w");
    fprintf(WriteFile,"#x y z n map H C\n");

    for (i = 0; i < xmax; i++) {
      for (j = 0; j < ymax; j++) {
	for (k = 0; k < zmax; k++) {

	n = ymax*zmax*i + zmax*j + k;
      
	if (map_in[n] == MAP_BOUNDARY) { fprintf(WriteFile,"%i %i %i %i %d %f %f\n", i, j, k, n, MAP_BOUNDARY, map_h[n], map_c[n]); }
	if (map_in[n] == MAP_FLUID)    { fprintf(WriteFile,"%i %i %i %i %d %f %f\n", i, j, k, n, MAP_FLUID, map_h[n], map_c[n]); }

	}  
      }  
    }
    fclose(WriteFile);
  }
    if (output_type == STATUS_WITH_H1_H2)  {
      sprintf(file,"Configuration_capillary.dat");
      WriteFile=fopen(file,"w");
      fprintf(WriteFile,"#x y z n map H1 H2\n");

      for (i = 0; i < xmax; i++) {
        for (j = 0; j < ymax; j++) {
      for (k = 0; k < zmax; k++) {

      n = ymax*zmax*i + zmax*j + k;
        
      if (map_in[n] == MAP_BOUNDARY) { fprintf(WriteFile,"%i %i %i %i %d %f %f\n", i, j, k, n, MAP_BOUNDARY, map_h[n], map_c[n]); }
      if (map_in[n] == MAP_FLUID)    { fprintf(WriteFile,"%i %i %i %i %d %f %f\n", i, j, k, n, MAP_FLUID, map_h[n], map_c[n]); }

      }
        }
      }
      fclose(WriteFile);
    }


  //printf("n = %d nsolid = %d nfluid = %d\n", xmax*ymax*zmax, nsolid,
	// xmax*ymax*zmax - nsolid);
  //changed for implementation of crystalline capillaries in order to see the volume fraction of the crystals
  printf("n = %d nsolid = %d nfluid = %d nsolid fraction: %f \n", xmax*ymax*zmax, nsolid,
	 xmax*ymax*zmax - nsolid, (double)nsolid/(xmax*ymax*zmax));
// start and end refers to first solid!
  printf("IN HEIGHT: start channel %i end channel %i \n", channel_start+1, channel_end+1);
  printf("IN LENGTH: start channel %i end channel %i \n", b0+1, b1+1);


  /* Write new data as char */

  fp_orig = fopen(outputfilename, "w");
  if (fp_orig == NULL) {
    printf("Cant open output\n");
    exit(-1);
  }

  for (i = 0; i < xmax; i++) {
    for (j = 0; j < ymax; j++) {
      for (k = 0; k < zmax; k++) {
	n = ymax*zmax*i + zmax*j + k;

	fprintf(fp_orig, "  %d", map_in[n]);
	if (output_type == STATUS_WITH_H) {
	  fprintf(fp_orig, "%f ", map_h[n]);
	}
	if (output_type == STATUS_WITH_C_H) {
	  
    fprintf(fp_orig, "%f %f ", map_c[n], map_h[n]);
	}
    if (output_type == STATUS_WITH_H1_H2) {
      fprintf(fp_orig, "%f %f ", map_c[n], map_h[n]);
    }

    if (output_type == STATUS_WITH_SIGMA) {
    if (map_sig[n] == 0) {
        fprintf(fp_orig, "  %.15e\n", map_sig[n]); // Stampa con due spazi per 0
    } else if (map_sig[n] < 0) {
        fprintf(fp_orig, " %.15e\n", map_sig[n]); // Stampa con uno spazio per numeri negativi
    } else { // Se è positivo
        fprintf(fp_orig, "  %.15e\n", map_sig[n]); // Stampa con due spazi per numeri positivi
    }
}
	
if (output_type == STATUS_ONLY) {
  fprintf(fp_orig, "\n");
}
	  


	  // TC MODIFY
	  /*
	if (i<xmax*fcr){
          fwrite(rho0, sizeof(double), 1, fp_orig);
          fwrite(kon, sizeof(double), 1, fp_orig);
          fwrite(koff, sizeof(double), 1, fp_orig);
	}
	else {
          fwrite(0.0, sizeof(double), 1, fp_orig);
	  fwrite(0.0, sizeof(double), 1, fp_orig);
	  fwrite(0.0, sizeof(double), 1, fp_orig);

	   }
	   */
	  // TC END Modify
	
      }
    }
  }

  fclose(fp_orig);

  free(map_in);
  free(map_c);
  free(map_h);
  free(map_sig);

  return 0;
}
