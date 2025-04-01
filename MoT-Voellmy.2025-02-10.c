/*******************************************************************************

  File:   MoT-Voellmy.2025-02-10.c             Dieter Issler, 2025-02-10

  Simulation of the motion of a gravity mass flow on a surface composed of
  regular quadrangles that project onto rectangles in the horizontal plane.
  In this version, a simplified variant of Fey's cell-centered Method of
  Transport is implemented, where the decomposition according to the
  eigenvectors of the Jacobians of the fluxes is omitted and only the convective
  mode is considered. The scheme is first-order in time.

  The code solves the shallow-water equations with friction, using the
  conservative formulation. The pressure distribution is hydrostatic and the
  earth pressure coefficients are 1. The Voellmy bed friction law is used in
  this version, but this is easy to modify. Entrainment can be included by
  modifying the source term for the flow height and removing the eroded mass
  from the bed.

  The program is started with the name of the command file as argument. In the
  latter, the names of the grid, initialization and output files are specified
  along with all simulation parameters. The grid file has to contain all
  pertinent geometrical information while the distribution of bed depth and
  strength, initial mass and velocity are specified in the initialization file.
  In this version, formatted ASCII input and ASCII or BinaryTerrain 1.3 output
  is used. The time steps are determined adaptively within given bounds. If
  negative flow heights occur, the CFL number is reduced and the time step
  repeated. The computation is limited to the area where the flow exceeds a
  (low) threshold and is stopped when the total flow momentum drops below a
  user-specified threshold value or the maximum simulated time is exceeded.

*******************************************************************************/


/** Compiler flags */
#define GCC
#if defined(__linux__) || defined(__APPLE__)
  #define LINUX
  #define DIRSEP        "/"
#endif
#ifdef _WIN32
  #define WINDOWS
  #define DIRSEP        "\\"
#endif

/** Code version */
#define VERSION         "2025-02-10"
#define INPUT_VERSION   "2024-09-10"

/** Include files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <locale.h>
#include <libgen.h>
#include <stdbool.h>
#include <sys/stat.h>

/** General constants and simple functions */

#define SQ(A)       ((A) * (A))
#define MAX(a,b)    ((a) > (b) ? (a) : (b))
#define MIN(a,b)    ((a) < (b) ? (a) : (b))

#define TRIES_MAX   30              /**< Max. # attempts at memory allocation */
#define TRY_WAIT    3               /**< Wait (s) between allocation attempts */

/** The following variables are declared before main(...) to make them
   accessible to other subroutines within the same file, in particular the
   ones responsible for reading the input file and for memory allocation. */

/** File names and handles */

char version[11] = VERSION;         /**< Designation of program version */
char topo_name[512];                /**< Name of avalanche path */
char run_name[512];                 /**< Name of this specific simulation */
char grid_fn[512];                  /**< Name of grid file */
char h_fn[512];                     /**< Name of release depth raster file*/
char u_fn[512];                     /**< Name of initial velocity raster file */
char v_fn[512];                     /**< Name of initial velocity raster file */
char b_fn[512];                     /**< Name of bed depth raster file */
char tauc_fn[512];                  /**< Name of bed shear strength raster file */
char mu_s_fn[512];                  /**< Name of bed frict. coeff. raster file */
char mu_fn[512];                    /**< Name of dry-friction-coefficient file */
char k_fn[512];                     /**< Name of drag-coefficient file */
char nD_fn[512];                    /**< Name of file with forest density nD */
char tD_fn[512];                    /**< Name of file with tree diameter tD */
char erod_tab_fn[512];              /**< Name of file with IsPa erosion table */
char out_fn[512];                   /**< Name of output file */
char max_fn[512];                   /**< Name of file with maximum values */
char comment[880-5*sizeof(int)];    /**< Can be used for comments */
char field_names[5][16];            /**< Array with names of output fields */
FILE *cfp;                          /**< Pointer to run set-up file */
FILE *gfp;                          /**< Pointer to grid file */


/** Variables concerning the numerics */

long   loc_dump_num;                /**< Location of # dumps in output file */
int    n_dump;                      /**< Number of time slices written */
int    nts = 0;                     /**< Time step number */
double t;                           /**< Simulated time (s) */
double dt;                          /**< Time step (s) */
double dt_min = 0.0001;             /**< Minimum admissible time step */
double dt_max = 0.2;                /**< Maximum admissible time step */
double t_max = 1000.0;              /**< Maximum time for a model run (s) */
double t_dump = 0.0;                /**< Next time for writing field values */
double dt_dump = 1.0;               /**< Write interval for results (s) */
double cfl = 0.7;                   /**< Courant-Friedrichs-Levy number */
double mov_vol;                     /**< Total volume in motion */
double h_lim = 5.0;                 /**< Max. effective flow depth for drag term,
                                       reasonable range is 5–10 m. */
double h_min = 0.05;                /**< Minimum flow height in active cells */
double u_min = 0.01;                /**< Minimum flow speed in active cells */
double mom_thr;                     /**< Stop criterion for flow momentum */
char   *fmt;                        /**< w: ESRI_ASCII_Grid, wb: BinaryTerrain */
char   write_vectors[4];            /**< Switch for writing u, v components */
char   write_max_press[4];          /**< Switch for writing maximum pressure */
char   write_press[4];              /**< Switch for writing press. time slices */
char   header[512];                 /**< Header of output raster files */
char   header_nD[512];              /**< Header of forest-permeability raster */

/** Grid-related variables */

int    m;                           /**< Number of grid nodes in W-E direction*/
int    n;                           /**< Number of grid nodes in S-N direction*/
int    i_min;                       /**< W boundary of active grid domain */
int    i_max;                       /**< E boundary of active grid domain */
int    j_min;                       /**< S boundary of active grid domain */
int    j_max;                       /**< N boundary of active grid domain */
double xllcorner;                   /**< x-coord. lower left corner of grid */
double yllcorner;                   /**< y-coord. lower left corner of grid */
double cellsize;                    /**< Projected size of square cell */
double **dx;                        /**< Oblique W-E length of a cell */
double **dy;                        /**< Oblique S-N length of a cell */
double **dA;                        /**< Oblique area of a cell */

/** Physical constants and material properties */

double g = 9.81;                    /**< Gravitational acceleration (m/s^2) */
double **gx;                        /**< x-comp. of gravitat. acceleration */
double **gy;                        /**< y-comp. of gravitat. acceleration */
double **gz0;                       /**< z-comp. of gravitat. acceleration */
double **gz;                        /**< Bed-normal gravitational acceleration
                                         incl. centrifugal acceleration */
double **IIxx;                      /**< Coefficient of 2nd fundamental form */
double **IIxy;                      /**< Coefficient of 2nd fundamental form */
double **IIyy;                      /**< Coefficient of 2nd fundamental form */
double **G_xy;                      /**< Rescaled off-diagonal metric tensor */
double rho   = 250.0;               /**< Flow density (kg/m^3) */
double rho_b = 200.0;               /**< Bed (snow cover) density */
double rho_d = 200.0;               /**< Deposit density */
double rrb = 1.25;                  /**< Ratio rho / rho_b */
double rrd = 1.25;                  /**< Ratio rho / rho_d */
double sigma = 1.0;                 /**< Factor in Grigorian–Ostroumov model */
double mu_g;                        /**< Tangent of friction angle in flow */
double mu_s0;                       /**< Constant bed friction coefficient */
double k_g;                         /**< Dimensionless turb. friction coeff. */
double kp = 1.0;                    /**< Passive earth-pressure coefficient */
double cD = 1.0;                    /**< Drag coefficient of a tree */
double nD_min = 0.001;              /**< Resid. forest density after destruction*/
double MoR = 5.0e7;                 /**< Modulus of rupture of trees (Pa) */
double decay_coeff = 0.1;           /**< Tree destruction coefficient (m/s) */
double h_drag = 0.0;                /**< Effective max. flow depth in drag term */
double k_erod;                      /**< Erosion coefficient à la RAMMS */
char   rheology[16];                /**< Type of friction law / rheology */
char   params[9];                   /**< Friction parameters can be "constant"
                                         or "variable" */
int    restart = 0;                 /**< Flag for non-zero initial velocity */
int    para = 0;                    /**< Const./variable friction param. (0/1) */
int    curve = 0;                   /**< Curvature effects (1/0) */
int    forest = 0;                  /**< Drag in forest (0/1/2) */
int    dyn_surf = 0;                /**< Dynamic surface flag (1/0) */
int    dep = 0;                     /**< Deposition flag (1/0) */
int    eromod;                      /**< Numerical code for erosion model */
int    grad;                        /**< Switch for shear strength profile */
long   utm_code;                    /**< Coded UTM zone, between -60 and 60 */
int    epsg;                        /**< EPSG code for geodetic datum */

/** Field variables, declared as pointers for dynamic memory allocation. */

double **h;                         /**< Flow height field */
double **u;                         /**< x-velocity field */
double **v;                         /**< y-velocity field */
double **s;                         /**< Speed field */
double **p_imp;                     /**< Impact pressure field */
double ***f_old;                    /**< Old values of conserved fields */
double ***f_new;                    /**< Time-stepped values of conserv. flds */
double ***src;                      /**< Area-integrated source terms */
double **d;                         /**< Field of deposition depth (m) */
double **z0;                        /**< Terrain surface elevation */
double **z;                         /**< Surface elevation incl. snow/deposit */
double **b;                         /**< Snow-cover (bed) depth (m) */
double **tau_c;                     /**< Bed shear strength (Pa) */
double **mu;                        /**< Dry-friction coefficient */
double **k;                         /**< Drag-friction coefficient */
double **mu_s;                      /**< Bed friction coefficient */
double **h_max;                     /**< Field of max. attained flow depths */
double **s_max;                     /**< Field of max. attained flow speeds */
double **p_max;                     /**< Field of max. attained impact press. */
double **u_max;                     /**< Field of max. attained u-velocity */
double **v_max;                     /**< Field of max. attained v-velocity */
double **b_min;                     /**< Minimum bed depth (due to erosion) */
double **d_max;                     /**< Max. deposition depth */
double **nD;                        /**< Field of forest opacity (n·D) */
double **tD;                        /**< Field of avg. tree diameter */
double **decay_const;               /**< Coeff. in tree fall-down rate */
float  *data;                       /**< Array holding data to be written */

/** Subroutines */

double **allocate2(int, int);       /**< Dynamically allocate 2D array */
double ***allocate3(int, int, int); /**< Dynamically allocate 3D array */
void   deallocate2(double **, int); /**< Deallocate 2D array */
void   deallocate3(double ***, int, int); /**< Deallocate 3D array */
void   allocate(void);              /**< Dynamically allocate arrays */
void   deallocate(void);            /**< Deallocate dynamic arrays*/
void   read_command_file(char *);   /**< Does what ist says! */
void   read_grid_file(char *);      /**< Reads DTM into array */
void   update_surface(double **);   /**< Adapts surface to erosion/deposition,
                                         (re-)calculates slope and curvature */
void   read_init_file(int, int);    /**< Set initial conditions from input */
int    read_raster(char *, double **, int, int, double, double, double,
                   double, int);    /**< Read data from AAIGrid file to array */
void   write_data(double, double **, double **, double **, double **,
                  double **, double **, double **, double **, double **,
                  int, int, int, int, int, char *);
                                    /**< Handles output from a time slice */
void   writeout(double **, char *, char *, int, int, int, int, char *,
                char *, char *);    /**< Writes output data to files */
void   primivar(double ***, double **, double **, double **, double **,
                double **);         /**< Computes primitive variables h,u,v,s
                                         from conserved fields h, hu, hv */
double ***source_terms(double **, double **, double **);
                                    /**< Computes source terms mass, momentum */
double find_dt(double **, double **, double **, double **, double);
                                    /**< Find new time step from CFL condition */
double update_boundaries(double ***);   /**< Determine new active region and
                                             quantity of movement */
void   create_dir(char *, char *);  /**< Create output directories as needed */


/*******************/
/*                 */
/*  Main function  */
/*                 */
/*******************/

int main(int argc, char *argv[])

{

  char   reason[80];
  int    i, j, p, q, di, dj;
  int    n_step = 0;
  int    repeat_flag;                   /**< Time step needs to be repeated */
  int    stop_code = 0;                 /**< Reason why simulation terminated */
  double aux, auy;                      /**< Auxiliary quantities */
  double U, V;                          /**< Local velocity components */
  double dAx, dAy, dAd;                 /**< Area fluxes to neighboring cells */
  double qhx, qhy, qhd;                 /**< Mass fluxes between cells */
  double qxx, qxy, qxd, qyx, qyy, qyd;  /**< Momentum fluxes to neighbor cells */
  double pWx, pEx, pSy, pNy;            /**< Earth pressure at cell boundaries */
  double F_drive_x, F_drive_y;          /**< Gravity and earth-pressure grad. */
  double F_drive_2, F_fric_2;           /**< Driving & retarding forces squared */
  double dir_cos, dir_sin;
  double mom_tot;                       /**< Approx. total avalanche momentum */
  double t_dmpp = 0.0;                  /**< Time of last write-out */


  printf("\n");
  printf("*****************************************************************\n");
  printf("*                                                               *\n");
  printf("*  MoT-Voellmy v. %10s                Dieter Issler, NGI  *\n",
         version);
  printf("*                                                               *\n");
  printf("*  Quasi-3D simulation of snow avalanches over complex terrain, *\n");
  printf("*  based on the Voellmy friction law and the cell-centered for- *\n");
  printf("*  mulation of the Method of Transport (with wave effects cur-  *\n");
  printf("*  rently neglected). Various erosion models are implemented.   *\n");
  printf("*  Curvature-induced friction, braking by/breaking of forest as *\n");
  printf("*  well as dynamic surface evolution can be simulated.          *\n");
  printf("*                                                               *\n");
  printf("*****************************************************************\n");
  printf("\n\n");

  if (argc != 2) {
    printf("   Usage:  MoT-Voellmy <input filename>\n\n");
    exit(3);
  }

  /* Set up the calculation. */

  read_command_file(argv[1]);
  read_grid_file(grid_fn);      /* Load z0 and reference raster header. */
  read_init_file(m, n);         /* Initializes all field variables, too. */
  printf("   main:  read_init_file completed.\n");

  t = 0.0;
  t_dump = -dt_dump;
  n_dump = 0;
  i_min = j_min = 0;
  i_max = m;                    /* m×n nodes, (m-1)×(n-1) cells! */
  j_max = n;
  strncpy(reason, "time limit was reached", 23);
  repeat_flag = 0;              /* repeat initially false   */

  if (dyn_surf) {
    for (i = i_min; i < i_max; i++)
      for (j = j_min; j < j_max; j++)
        z[i][j] = z0[i][j];
  }


  /* Time loop: */

  while (t < t_max) {

    printf("   main:  Step %5d,  t = %8.4f s,  %7.0f m^3,  [%d,%d]x[%d,%d]\n",
           n_step, t, mov_vol, i_min, i_max, j_min, j_max);
    if (t >= t_dump + dt_dump && t_max >= dt_dump) {
      write_data(t, h, h, b, d, s, u, v, p_imp, nD,
                 i_min, i_max, j_min, j_max, 1, fmt);
      t_dmpp = t;
      t_dump += dt_dump;
      n_dump++;
    }

    if (curve == 0)                     /* Without curvature effects */
      for (i = i_min; i < i_max; i++)
        for (j = j_min; j < j_max; j++)
          memcpy(f_old[i][j], f_new[i][j], 3*sizeof(double));
          /**< NB. f_old is needed in case the timestep needs to be repeated. */

    else if (curve == 1) {              /* With curvature effects */
      for (i = i_min; i < i_max; i++)
        for (j = j_min; j < j_max; j++) {
          /**< NB. f_old is needed in case the timestep needs to be repeated. */
          memcpy(f_old[i][j], f_new[i][j], 3*sizeof(double));

          /* Normal force corrected for curvature effects, with gz limited to
             non-negative values to prevent lift-off on convex terrain.
             Contributed by Hervé Vicari, 2023. */
          U = u[i][j]; V = v[i][j];
          gz[i][j] \
            = MAX(0.0,
                  gz0[i][j]
                    + (IIxx[i][j]*U*U + IIyy[i][j]*V*V + 2.0*IIxy[i][j]*U*V)
                      / MAX(0.0001, U*U + V*V + 2.0*G_xy[i][j]*U*V));
        }
    }

    if ((dt = find_dt(h, u, v, gz, cfl)) < dt_min) {
      strncpy(reason, "timestep fell below lower bound", 32);
      stop_code = 2;
      printf("   main:  dt set to %.5f s.\n", dt);
      break;                            /* Leave time loop to shut down. */
    }

    src = source_terms(h, u, v);

    for (i = i_min; i < i_max; i++) {

      if (repeat_flag == 1) {
        printf(".");
        for (p = i_min; p < i_max; p++) /* Restore old field values */
          for (q = j_min; q < j_max; q++)
            memcpy(f_new[p][q], f_old[p][q], 3*sizeof(double));
        repeat_flag = 0;                /* Start i-loop again with reduced   */
        i = i_min - 1;                  /* timestep and f_new reset to f_old */
        dt *= 0.8;                      /* by resetting i and `continue'.    */
        if (dt < dt_min) {
          strncpy(reason, "timestep fell below lower bound", 32);
          repeat_flag = -1;             /* Signals failure by setting flag. */
          stop_code = 2;		        /* Use as exit code at shut-down. */
          break;                        /* Break out of i-loop; flag triggers*/
                                        /* break statement below i-loop. */
        }
        continue;                       /* Start the i-loop over again; i is */
                                        /* incremented from i_min-1 to i_min */
      }

      for (j = j_min; j < j_max; j++) {

        /* Quantities used in all field components: */
        di = (u[i][j] >= 0.0 ? 1 : -1);
        dj = (v[i][j] >= 0.0 ? 1 : -1);
        aux = fabs(u[i][j]) * dt;
        auy = fabs(v[i][j]) * dt;
        dAx = aux * (dy[i][j] - auy);   /* Area flowing out in x-direction */
        dAy = auy * (dx[i][j] - aux);   /* Area flowing out in y-direction */
        dAd = aux * auy;                /* Outflow in diagonal direction */

        /* Bed depth limits erosion, flow depth limits deposition: */
        if (eromod > 0 && src[i][j][0] > 0.0) {     /* Erosion */
          /* Check erosion rate limit: */
          src[i][j][0] = MIN( src[i][j][0], b[i][j]*dA[i][j]/(rrb*dt) );
          /* Update erosion reservoir: */
          b[i][j] = MAX(0.0, b[i][j] - src[i][j][0]*rrb*dt/dA[i][j]);
          /* MAX(...) used to prevent spurious −0.0 rounding errors.
             Contributed by Hervé Vicari and Callum Tregaskis. */
        }
        else if (dep > 0 && src[i][j][0] < 0) {     /* Deposition */
          /* Check deposition rate limit: */
          src[i][j][0] = MAX( src[i][j][0], -f_old[i][j][0]/dt );
          /* Update deposit reservoir: */
          d[i][j] -= src[i][j][0] * rrd * dt / dA[i][j];
        }
        else src[i][j][0] = 0.0;

        /* Advective mass fluxes: */
        qhx = h[i][j] * dAx;
        qhy = h[i][j] * dAy;
        qhd = h[i][j] * dAd;

        /* Advective momentum fluxes: */
        qxx = qhx * u[i][j];
        qxy = qhy * u[i][j];
        qxd = qhd * u[i][j];
        qyx = qhx * v[i][j];
        qyy = qhy * v[i][j];
        qyd = qhd * v[i][j];

        f_new[i][j][0] -= (qhx + qhy + qhd - src[i][j][0]*dt);
        f_new[i][j][1] -= qxx + qxy + qxd;      /* Flowing out of cell (i,j) */
        f_new[i][j][2] -= qyx + qyy + qyd;

        if (i+di >= 0 && i+di < m) {
          f_new[i+di][j][0] += qhx;             /* Can be ahead or behind */
          f_new[i+di][j][1] += qxx;             /* (i,j) depending on di  */
          f_new[i+di][j][2] += qyx;
        }
        if (j+dj >= 0 && j+dj < n) {
          f_new[i][j+dj][0] += qhy;
          f_new[i][j+dj][1] += qxy;
          f_new[i][j+dj][2] += qyy;
        }
        if (i+di >= 0 && i+di < m && j+dj >= 0 && j+dj < n) {
          f_new[i+di][j+dj][0] += qhd;
          f_new[i+di][j+dj][1] += qxd;
          f_new[i+di][j+dj][2] += qyd;
        }

        /* Test for negative flow heights: */
        if (f_new[i][j][0] < 0.0) {
          repeat_flag = 1;
          break;                        /* Break out of j-loop, start next */
        }                               /* iteration of i-loop. */

        /* Momentum fluxes due to pressure gradients:
           Need to distinguish between empty cells (no pressure transmission),
           non-empty cells in movement (friction forces fully activated), and
           non-empty cells at rest, where the static friction force may or may
           not be fully activated. */

        /* Pressure gradient in x-direction */
        if (i > i_min && i < i_max-1) {
          pEx = 0.25 * kp * dy[i+1][j]
                * (gz[i][j]+gz[i+1][j]) * h[i][j]*h[i+1][j];
          pWx = 0.25 * kp * dy[i][j]
                * (gz[i-1][j]+gz[i][j]) * h[i-1][j]*h[i][j];
        }
        /* Von Neumann boundary conditions for outermost cells: */
        else if (i == i_min)
          pWx = pEx = 0.25 * kp * dy[i+1][j]
                      * (gz[i][j]+gz[i+1][j]) * h[i][j]*h[i+1][j];
        else if (i == i_max-1)
          pEx = pWx = 0.25 * kp * dy[i][j]
                      * (gz[i-1][j]+gz[i][j]) * h[i-1][j]*h[i][j];
        /* Include (superfluous) ELSE for sake of some compilers. */
        else
          pEx = pWx = 0.0;

        /* Pressure gradient in y-direction */
        if (j > j_min && j < j_max-1) {
          pNy = 0.25 * kp * dx[i][j+1]
                * (gz[i][j]+gz[i][j+1]) * h[i][j]*h[i][j+1];
          pSy = 0.25 * kp * dx[i][j]
                  * (gz[i][j-1]+gz[i][j]) * h[i][j-1]*h[i][j];
        }
        /* Von Neumann boundary conditions for outermost cells: */
        else if (j == j_min)
          pSy = pNy = 0.25 * kp * dx[i][j+1]
                      * (gz[i][j]+gz[i][j+1]) * h[i][j]*h[i][j+1];
        else if (j == j_max-1)
          pNy = pSy = 0.25 * kp * dx[i][j]
                      * (gz[i][j-1]+gz[i][j]) * h[i][j-1]*h[i][j];
        /* Include (superfluous) ELSE for sake of some compilers. */
        else
          pNy = pSy = 0.0;

        /* Test whether non-empty cells at rest will start moving. */

        if (s[i][j] <= u_min) {
          /* Gravity and pressure gradient combined: */
          F_drive_x = gx[i][j] * f_old[i][j][0] + pWx - pEx;
          F_drive_y = gy[i][j] * f_old[i][j][0] + pSy - pNy;
          F_drive_2 = SQ(F_drive_x) + SQ(F_drive_y)
                      + 2.0 * G_xy[i][j] * F_drive_x * F_drive_y;

          /* Is dry friction fully activated? If not, the driving and
             resisting forces cancel and there is no need to add to f_new. */
          F_fric_2 = SQ(mu[i][j] * gz[i][j] * f_old[i][j][0]);
          if (F_drive_2 > F_fric_2) {
            dir_cos = F_drive_x / sqrt(F_drive_2);
			dir_sin = F_drive_y / sqrt(F_drive_2);
			f_new[i][j][1] += (F_drive_x - dir_cos*sqrt(F_fric_2)) * dt;
            f_new[i][j][2] += (F_drive_y - dir_sin*sqrt(F_fric_2)) * dt;
          }
        }
        else {
          f_new[i][j][1] += (pWx - pEx + src[i][j][1]) * dt;
          f_new[i][j][2] += (pSy - pNy + src[i][j][2]) * dt;
        }
      }
    }

    if (repeat_flag == -1) break;       /* Break out of time loop. */

    /* If the momentum vector in a cell reverses direction, arrest
       the cell unless the new direction is downhill: */
    for (i = i_min; i < i_max; i++)
      for (j = j_min; j < j_max; j++)
        if (f_old[i][j][1]*f_new[i][j][1]
                  + f_old[i][j][2]*f_new[i][j][2] < 0.0
            && f_new[i][j][1]*gx[i][j] + f_new[i][j][2]*gy[i][j] < 0.0) {
          if (dep == 1) {
            d[i][j] += f_new[i][j][0] / dA[i][j];
            f_new[i][j][0] = 0.0;
          }
          f_new[i][j][1] = 0.0;
          f_new[i][j][2] = 0.0;
        }

    /* Update surface elevation for dynamic bed computation: */
    if (dyn_surf) {
      for (i = i_min; i <= i_max; i++)
        for (j = j_min; j <= j_max; j++)
          /* Add bed layer and deposited layer to surface. */
          z[i][j] = z0[i][j] + (b[i][j] + d[i][j]) * g / gz0[i][j];
      update_surface(z);
    }

    /* Update boundaries and test if avalanche still moves: */
    primivar(f_new, h, u, v, s, p_imp);
    mom_tot = update_boundaries(f_new);
    if (mom_tot < mom_thr && n_step > 10) {
      strncpy(reason, "avalanche has stopped or left the domain", 43);
      stop_code = 1;
      break;                            /* Break out of time loop. */
    }

    /* Update time: */
    t += dt;
    n_step++;
  }
  printf("   main:  Finished time loop.\n");

  /* End of time loop */

  /* Write out last time step only if there is new data! */
  if (t > t_dmpp && t_max >= dt_dump)
    write_data(t, d, h, b, d, s, u, v, p_imp, nD, 0, m, 0, n,
               1, fmt);

  /* Write maximum fields over entire simulation (incl. deposit depth). */
  write_data(t, d, h_max, b_min, d_max, s_max, u_max, v_max, p_max, nD,
             0, m, 0, n, 2, fmt);

  deallocate();
  printf("\n   Simulation terminated because %s.\n\n", reason);

  exit(stop_code);

}

/**********************/
/*  End of main(...)  */
/**********************/


/*******************/
/*                 */
/*  primivar(...)  */
/*                 */
/*******************/

/** Computes the primitive variables h, u, v from the conservative quantities
   h dA, h u dA, h v dA. Also, the speed s is computed for a non-orthogonal
   coordinate system. */

void primivar(double ***f, double **h, double **u, double **v,
              double **s, double **p)

{
  int    i, j;
  double aux1, aux2;

  for (i = i_min; i < i_max; i++) {
    for (j = j_min; j < j_max; j++) {
      aux1 = 1.0 / dA[i][j];
      aux2 = (f[i][j][0] > 0.0 ? 1.0 / MAX(f[i][j][0], h_min*dA[i][j]) : 0.0);
      h[i][j] = f[i][j][0] * aux1;
      u[i][j] = f[i][j][1] * aux2;
      v[i][j] = f[i][j][2] * aux2;
      p[i][j] = SQ(u[i][j]) + SQ(v[i][j]) + 2.0*G_xy[i][j]*u[i][j]*v[i][j];
      s[i][j] = sqrt(p[i][j]);
      p[i][j] *= (0.001*rho);           /* Pressures in kPa */
    }
  }
}

/**************************/
/*  End of primivar(...)  */
/**************************/


/***********************/
/*                     */
/*  source_terms(...)  */
/*                     */
/***********************/

/** Computes source terms in the equations for the conservative fields.
   For the erosion rate in the TJEM model, note that tau_c is scaled with
   the flow density in read_init_file(). */

double ***source_terms(double **h, double **u, double **v)

{
  int    i, j, variant;
  double speed, dir_cos, dir_sin, cos_th;
  double tau_b, tau_c_loc;              /* Bed shear stress over density,
                                           local bed shear strength */
  double mu_loc, k_loc;                 /* Including forest effects */
  double hs = 0.0;                      /* Snow cover depth for torque */
  double bend_mom;                      /* Bending moment on tree */
  double dp;                            /* Pressure à la Grigorian–Ostroumov */
  double U, V, gxy;                     /* Local velocity, off-diag. metric */
  double dbdx, dbdy;                    /* Change of snow depth in x, y-dir. */
  double calpha, salpha, talpha;        /* Cosine, sine, tangent of snow surface
                                           slope angle in flow direction */

  variant = 2*para + forest;            /* 0: constant, no forest
                                           1: constant, with forest
                                           2: variable, no forest
                                           3: variable, with forest */

  for (i = i_min; i < i_max; i++) {
    for (j = j_min; j < j_max; j++) {

      /* Local values that will come in handy: */
      speed = s[i][j];
      U = u[i][j];
      V = v[i][j];
      gxy = G_xy[i][j];
      cos_th = SQ(cellsize) / dA[i][j];

      /* Set friction parameters according to chosen variant: */
      switch(variant) {
        case 0 :
          mu_loc = mu_g;
          k_loc = k_g;
          break;
        case 1 :
          mu_loc = mu_g + 1.25 * cos_th * nD[i][j]*h[i][j];
          k_loc = k_g + 0.5*cD*cos_th*nD[i][j]*h[i][j];
          break;
        case 2 :
          mu_loc = mu[i][j];
          k_loc = k[i][j];
          break;
        case 3 :
          mu_loc = mu[i][j] + 1.25 * cos_th * nD[i][j]*h[i][j];
          k_loc = k[i][j] + 0.5*cD*cos_th*nD[i][j]*h[i][j];
          break;
        default:                /* Cannot be reached, for compiler's sake. */
          printf("\nIllegal value %d of \'variant\' --- STOP!\n\n", variant);
          exit(21);
      }

      /* Option to increase the drag term by a factor so that it does not
         vanish if the flow depth becomes excessive in channelized areas. No
         changes as h → 0, but as h → ∞, the drag deceleration is like for
         h = h_drag. */
      if (h_drag > 0.0)
        k_loc /= (1.0 - exp(-h_drag / MAX(h[i][j], h_min)));

      /* Maximum shear stress at the bottom of the flow and shear strength
         of bed including Coulombic contribution. This must be computed before
         adding the braking effect of forest: */
      tau_b = mu_loc*gz[i][j]*h[i][j] + k_loc*SQ(speed);

      /* Erosion term: */
      switch(eromod) {
        case 0 :                        /* No erosion */
          src[i][j][0] = 0.0;
          break;
        case 1 :                        /* RAMMS erosion model */
          if (h[i][j] > h_min && speed > 1.0)
            src[i][j][0] = k_erod * speed * dA[i][j];
          else
            src[i][j][0] = 0.0;
          break;
        case 2 :                        /* Tangential-jump erosion model */
          /* grad = 0 if snow-cover strength assumed constant with depth,
             grad = 1 if vertical strength gradient is spatially constant,
             grad = 2 if vertical strength gradient is read from file. */
          tau_c_loc = (grad < 2 ? tau_c[i][j] + mu_s0*gz[i][j]*h[i][j] \
                                : tau_c[i][j] + mu_s[i][j]*gz[i][j]*h[i][j]);
          /* Erosion rate prop. to the excess of rheological stress over bed
             shear strength: */
          src[i][j][0] = (speed > 10.0*u_min && h[i][j] > 10.0*h_min ? \
                          MAX(0.0, tau_b - tau_c_loc) * dA[i][j] / speed : \
                          0.0);
          /* In the TJEM, the bed shear stress is limited to τ_c if erodible
             bed material is present: */
          if (src[i][j][0] > 0.0 && b[i][j] > 0.0)
            tau_b = MAX(tau_c_loc, tau_b);
          break;
        case 3 :                        /* com1DFA erosion model (AvaFrame) */
          /* τ_c is here taken to represent the specific erosion energy e_b
             in the com1DFA entrainment module. It has the same dimensions
             m²/s² as the specific bed shear stress μ·g_z·h + k·u². There are
             no indicative values quoted in the com1DFA manual, but one may
             assume that e_b is roughly 100 times larger than typical values
             of μ·g_z·h + k·u², i.e., in the range 300–3000 m²/s². */
          if (h[i][j] > h_min && speed > 1.0)
            src[i][j][0] = speed * dA[i][j] / tau_c[i][j] \
                           * (mu[i][j]*gz[i][j]*h[i][j] + k[i][j]*SQ(speed));
          else
            src[i][j][0] = 0.0;
          break;
        case 4 :                        /* Grigorian–Ostroumov  (GOEM) */
          /* Gradient of snow surface relative to terrain: */
          dbdx = (i > 0 && i < m-1 ? 0.5 * (b[i+1][j]-b[i-1][j]) / dx[i][j] \
                                   : (i == 0 ? (b[1][j]-b[0][j]) / dx[0][j] \
                                             : (b[m-1][j]-b[m-2][j])/dx[m-2][j]
                                     )
                 );
          dbdy = (j > 0 && j < n-1 ? 0.5 * (b[i][j+1]-b[i][j-1]) / dy[i][j] \
                                   : (j == 0 ? (b[i][1]-b[i][0]) / dy[i][0] \
                                             : (b[i][n-1]-b[i][n-2])/dy[i][n-2]
                                     )
                 );
          /* Slope angle of snow surface rel. to terrain in flow direction: */
          talpha = ((U + V*gxy)*dbdx + (V + U*gxy)*dbdy) / MAX(0.01, speed);
          calpha = 1.0 / sqrt(1.0 + SQ(talpha));
          salpha = talpha * calpha;
          /* Excess pressure dp and strength τ_c are scaled by ρ! Include
             depth-dependent bed strenfromgth as in TJEM */
          tau_c_loc = (grad < 2 ? tau_c[i][j] + mu_s0*gz[i][j]*h[i][j] \
                                : tau_c[i][j] + mu_s[i][j]*gz[i][j]*h[i][j]);
          dp = MAX(0.0, gz[i][j]*h[i][j]*calpha + k_erod*SQ(speed)*salpha \
                        - tau_c[i][j]);
          src[i][j][0] = sigma * sqrt(dp) * dA[i][j] * calpha;
          break;
        default :                       /* To satisfy purists... */
          printf("\n   Erosion model #%d not implemented. STOP!\n\n", eromod);
          exit(29);
      }

      /* Momentum sources (gravity and friction): */
      if (speed > u_min) {              /* Friction opposing flow direction */
        dir_cos = u[i][j] / speed;
        dir_sin = v[i][j] / speed;
        src[i][j][1] = (gx[i][j]*h[i][j] - dir_cos*tau_b) * dA[i][j];
        src[i][j][2] = (gy[i][j]*h[i][j] - dir_sin*tau_b) * dA[i][j];
      }

      /* What is the fate of the forest?
         forest=0: no forest; forest=1: braking effect, can be destroyed
         If forest density is below residual value, do nothing. Otherwise: */
      if ((forest == 1) && (nD[i][j] > nD_min)) {
        /* If no erosion, assume 1.0 m snow depth for moment calculation: */
        hs = (eromod > 0 ? b[i][j] : 1.0);
        /* Check whether the forest in the cell is still intact: */
        if (decay_const[i][j] == 0.0) {         /* No destruction yet */
          /* Compare bending moment on average tree to its strength: */
          bend_mom = 0.25 * cD * rho * (SQ(speed) + 5.0*g*h[i][j]*cos_th)
                     * tD[i][j] * h[i][j] * (h[i][j] + 2.0*hs);
          if (bend_mom > MoR * tD[i][j]*tD[i][j]*tD[i][j])
            /* Tree breaks or is uprooted. Estimate time for it to fall and
               lose braking effect on avalanche: */
            decay_const[i][j] = decay_coeff / tD[i][j];
        }
        else {                                  /* Destruction ongoing */
          /* Approximate an exponential decay by first two terms: */
          nD[i][j] *= MAX(0.0, 1.0 - decay_const[i][j]*dt);
        }
      }
    }
  }

  return(src);
}

/******************************/
/*  End of source_terms(...)  */
/******************************/


/******************/
/*                */
/*  find_dt(...)  */
/*                */
/******************/

/** The time step is computed on the basis of the maximum velocity of the
   forward "acoustic" wave. Note that source terms have the potential to
   "empty" a cell more quickly; the main routine checks for this and repeats
   a time step with reduced dt if necessary. */

double find_dt(double **h, double **u, double **v, double **gz, double cfl)

{
  int i, j;
  double dt = 1000.0, aux;

  for (i = i_min; i < i_max; i++) {
    for (j = j_min; j < j_max; j++) {
      aux = MAX(sqrt(SQ(u[i][j])+SQ(v[i][j])) + sqrt(gz[i][j]*h[i][j]), u_min);
      dt = MIN(cfl * MIN(dx[i][j], dy[i][j]) / aux,  dt);
    }
  }
  dt = MIN(dt, dt_max);

  return(dt);
}

/*************************/
/*  End of find_dt(...)  */
/*************************/


/****************************/
/*                          */
/*  update_boundaries(...)  */
/*                          */
/****************************/

/** This function is called to (i) determine whether the flow has effectively
   stopped, (ii) update the fields of maximum values, and (iii) determine the
   "active" region of the computational domain, i.e., the smallest rectangle
   outside of which the flow height or the speed are below user-specified
   thresholds. (One row of cells is added in every direction to prevent
   spurious effects.) */

double update_boundaries(double ***f)

{
  int    i, j, west = m, east = 0, south = n, north = 0;
  double mom = 0.0, speed, vol_min, tot_vol = 0.0;

  mov_vol = 0.0;
  for (i = i_min; i < i_max; i++) {
    for (j = j_min; j < j_max; j++) {
      vol_min = h_min * dA[i][j];
      speed = s[i][j];

      /* Boundaries of active domain */
      if (f[i][j][0] > vol_min && speed > u_min) {
        west  = MIN(west, i-1);
        east  = MAX(east, i+1);
        south = MIN(south, j-1);
        north = MAX(north, j+1);
        mov_vol += f[i][j][0];
      }

      /* Maximum fields, written to file(s) at end of run. Note u_max, v_max
         are components of max. speed and not actual max(u) and max(v). */
      h_max[i][j]   = MAX(h_max[i][j], h[i][j]);
      if (speed > s_max[i][j]) {
        s_max[i][j] = speed;
        u_max[i][j] = u[i][j];
        v_max[i][j] = v[i][j];
        p_max[i][j] = 0.001 * rho * SQ(speed);
      }
      mom += speed * f[i][j][0];

      if (eromod > 0)                   /* Update erodible snow depth */
        b_min[i][j] = MIN(b[i][j], b_min[i][j]);
      if (dep > 0)
        d_max[i][j] = MAX(d[i][j], d_max[i][j]);
    }
  }

  i_min = MAX(0, west);
  i_max = MIN(m, east+1);             /* m x n nodes, (m-1) x (n-1) cells */
  j_min = MAX(0, south);
  j_max = MIN(n, north+1);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      tot_vol += f[i][j][0];
  printf("      update_boundaries:  V_tot = %9.1f m³,  J_tot = %9.0f kg m/s\n",
         tot_vol, rho*mom);

  return(mom);
}

/***********************************/
/*  End of update_boundaries(...)  */
/***********************************/


/****************************/
/*                          */
/*  read_command_file(...)  */
/*                          */
/****************************/

/** The command file given on the command line is parsed for the values of
   variables and some basic consistency checks are performed. The parameters
   have to appear in a specific order and must be preceded by the correct
   keyword, but the amount of white space is free. */

void read_command_file(char *ifn)

{
  FILE   *ifp;                      /* Pointer to input file handle */
  char   file_version[20], tstr[20];
  char   out_path[512];             /* Root directory for the run */
  char   utm_str[5];                /* String with number of UTM zone */
  char   *tail;                     /* Ptr. to hemisphere attr. of UTM_str */
  char   output_format[32];         /* String designating output format
                                       ESRI_ASCII_Grid or Binary_Terrain */
  char   erosion[20];               /* Type of entrainment law: "none", "RAMMS",
                                       "TJEM", "AvaFrame", "GOEM". */
  char   gradient[16];              /* Local/global bed friction coefficient */
  char   dep_flag[16];              /* Flag for deposition (no, yes)*/
  char   dyn_surf_switch[16];       /* Flag for dynamical surface (no, yes) */
  char   curveff[4];                /* Curvature effects (yes/no) */
  char   foresteff[8];              /* Drag in forest (yes/no) */
  char   line[512];                 /* Hold a line of input file */
  int    ifv = 0;                   /* Input file version, 0 if "2024-09-10",
                                       1 if "2021-10-25", 2 if "2020-06-23" */
  int    n_items;                   /* # data lines in input file */
  int    lest = 0;                  /* Counts # variables read */
  int    dummy = 0;                 /* Counts # irrelevant chars read. */

  if ((ifp = fopen(ifn, "r")) == NULL) {
    printf("   Failed to open %s. STOP!\n\n", ifn);
    exit(10);
  }

  /* Read parameter values.
     This is a less safe scanning routine than in SL-1D, but it should do
     for this purpose. */

  if (fgets(line, 512, ifp) == NULL) {
    printf("\n   Failed to read first line of %s. STOP!\n\n", ifn);
    exit(10);
  }
  if (!strncmp(line, "# Run information", 17)) {
    dummy += fscanf(ifp, "#\n");            /* Read over "#\n" and get line */
    if (fgets(line, 512, ifp) == NULL) {    /* with version info. */
      printf("\n   Failed to read a line of %s. STOP!\n\n", ifn);
      exit(10);
    }
  }
  /* Check whether input file format is newest version (ifv = 0): */
  sscanf(line, "MoT-Voellmy input file version %[^\r\n]\n", file_version);
  if (strncmp(file_version, INPUT_VERSION, 10)) {
    /* Accept file formats 2021-10-25 (ifv = 1) and 2020-06-23 (ifv = 2) but
       reject earlier ones: */
    if (!strncmp(file_version, "2021-10-25", 10))
      ifv = 1;
    else if (!strncmp(file_version, "2020-06-23", 10))
      ifv = 2;
    else {
      printf("   Input file format version %s not supported. STOP!\n\n",
             file_version);
      exit(11);
    }
  }
  lest += fscanf(ifp, "Area of Interest %[^\r\n]\n", topo_name);
  printf("%2d  topo_name       = %s\n", lest, topo_name);
  lest += fscanf(ifp, "UTM zone %[^\r\n]\n", utm_str);
  printf("%2d  utm_str         = %s\n", lest, utm_str);
  lest += fscanf(ifp, "EPSG geodetic datum code %d\n", &epsg);
  printf("%2d  epsg            = %d\n", lest, epsg);
  lest += fscanf(ifp, "Run name %[^\r\n]\n#\n", run_name);
  printf("%2d  run_name        = %s\n", lest, run_name);
  if (ifv < 2)
    dummy += fscanf(ifp, "# File names\n#\n");
  lest += fscanf(ifp, "Grid filename %[^\r\n]\n", grid_fn);
  printf("%2d  grid_fn         = %s\n", lest, grid_fn);
  lest += fscanf(ifp, "Release depth filename %[^\r\n]\n", h_fn);
  printf("%2d  h_fn            = %s\n", lest, h_fn);
  lest += fscanf(ifp, "Bed depth filename %[^\r\n]\n", b_fn);
  printf("%2d  b_fn            = %s\n", lest, b_fn);
  lest += fscanf(ifp, "Bed shear strength filename %[^\r\n]\n", tauc_fn);
  printf("%2d  tauc_fn         = %s\n", lest, tauc_fn);
  if (ifv < 2) {                        /* New: forest destruction is option */
    lest += fscanf(ifp, "Forest density filename %[^\r\n]\n", nD_fn);
    printf("%2d  nD_fn           = %s\n", lest, nD_fn);
    lest += fscanf(ifp, "Tree diameter filename %[^\r\n]\n", tD_fn);
    printf("%2d  tD_fn           = %s\n", lest, tD_fn);
  }
  else {                                /* Old: no forest destruction option */
    lest += fscanf(ifp, "Forest data filename %[^\r\n]\n", nD_fn);
    printf("%2d  nD_fn           = %s\n", lest, nD_fn);
  }
  lest += fscanf(ifp, "Start velocity u filename %[^\r\n]\n", u_fn);
  printf("%2d  u_fn            = %s\n", lest, u_fn);
  lest += fscanf(ifp, "Start velocity v filename %[^\r\n]\n", v_fn);
  printf("%2d  v_fn            = %s\n", lest, v_fn);

  lest += fscanf(ifp, "Output filename root %[^\r\n]\n", out_fn);
  printf("%2d  out_fn          = %s\n", lest, out_fn);
  lest += fscanf(ifp, "Output format %[^\r\n]\n#\n", output_format);
  if (!strncmp(output_format, "ESRI_ASCII_Grid", 16))
    fmt = "w";
  else if (!strncmp(output_format, "Binary_Terrain", 15))
    fmt = "wb";
  else {
    printf("   read_command_file:  Output format %s not supported. STOP!\n\n",
           output_format);
    exit(12);
  }
  printf("%2d  fmt             = %s\n", lest, fmt);

  if (ifv < 2)
    dummy += fscanf(ifp, "# Physical parameters\n#\n");
  /* Read the value of the gravitational constant as a string, check for the
     decimal sign (period or comma) and set the locale to "C" for the former
     and to "nb-NO" for the latter. */
  lest += fscanf(ifp, "Gravitational acceleration (m/s^2) %[^\r\n]\n", tstr);
  if (strchr(tstr, '.') != NULL) {
    setlocale(LC_NUMERIC, "C");         /* Locale is guaranteed to exist. */
    printf("%2d  read_command_file:  LC_NUMERIC set to C.\n", lest);
  }
  #ifdef WINDOWS
  else if (strchr(tstr, ',') != NULL) {
    if (setlocale(LC_NUMERIC, "no_NO.UTF-8") != NULL)
      printf("   read_command_file:  LC_NUMERIC set to no_NO.UTF-8.\n");
    else {
      printf("   read_command_file:  Cannot switch to nb-NO. STOP!\n\n");
      exit(14);
    }
  }
  #endif /* defined WINDOWS */
  #ifdef LINUX
  else if (strchr(tstr, ',') != NULL) {
    if (setlocale(LC_NUMERIC, "nb_NO.UTF-8") != NULL)
      printf("   read_command_file:  LC_NUMERIC set to nb_NO.\n");
    else {
      printf("   read_command_file:  Cannot switch to nb_NO.UTF-8. STOP!\n\n");
      exit(14);
    }
  }
  #endif /* defined LINUX */
  else {
    printf("   read_command_file:  No decimal sign in value of g. STOP!\n\n");
    printf("                       g = %s\n", tstr);
    exit(15);
  }
  sscanf(tstr, "%lf", &g);
  printf("%2d  g               = %.2f\n", lest, g);

  /* Continue to read normally. */
  if (ifv == 0) {                       /* Changes in v.2024-09-10! */
    lest += fscanf(ifp, "Flow density (kg/m^3) %lf\n", &rho);
    printf("%2d  rho             = %.3f\n", lest, rho);
    lest += fscanf(ifp, "Bed density (kg/m^3) %lf\n", &rho_b);
    printf("%2d  rho_b           = %.3f\n", lest, rho_b);
    lest += fscanf(ifp, "Deposit density (kg/m^3) %lf\n", &rho_d);
    printf("%2d  rho_d           = %.3f\n", lest, rho_d);
  }
  else {                                /* Versions 2020-06-23, 2021-10-25 */
    lest += fscanf(ifp, "Density (kg/m^3) %lf\n", &rho);
    printf("%2d  rho             = %.3f\n", lest, rho);
    rho_d = 1.6*rho;
    printf("    rho_d not specified, set to %.3f kg/m³.\n", rho_d);
  }
  lest += fscanf(ifp, "Rheology %[^\r\n]\n", rheology);
  printf("%2d  rheology        = %s\n", lest, rheology);
  lest += fscanf(ifp, "Parameters %[^\r\n]\n", params);
  printf("%2d  params          = %s\n", lest, params);
  if (!strncmp(params, "constant", 9)) {
    para = 0;
    lest += fscanf(ifp, "Dry-friction coefficient (-) %lf\n", &mu_g);
    lest += fscanf(ifp, "Turbulent drag coefficient (-) %lf\n", &k_g);
    printf("%2d  mu              = %5.3f\n", lest-1, mu_g);
    printf("%2d  k               = %6.4f\n", lest, k_g);
  }
  else {
    para = 1;
    lest += fscanf(ifp, "Dry-friction coefficient (-) %[^\r\n]\n", mu_fn);
    lest += fscanf(ifp, "Turbulent drag coefficient (-) %[^\r\n]\n#\n", k_fn);
    printf("%2d  mu_fn           = %s\n", lest-1, mu_fn);
    printf("%2d  k_fn            = %s\n", lest, k_fn);
  }

  if (ifv < 2) {                        /* Option of modified drag term */
    lest += fscanf(ifp, "Effective drag height (m) %lf\n", &h_drag);
    printf("%2d  h_drag          = %.1f\n", lest, h_drag);
  }
  else                                  /* Standard drag term */
    h_drag = 0.0;

  lest += fscanf(ifp, "Centrifugal effects %s\n", curveff);
  printf("%2d  curveff         = %s\n", lest, curveff);
  if (!strncmp(curveff, "yes", 4) || !strncmp(curveff, "Yes", 4))
    curve = 1;

  if (ifv < 2) {
    lest += fscanf(ifp, "Forest effects %[^\r\n]\n", foresteff);
    printf("%2d  foresteff       = %s\n", lest, foresteff);
  }

  lest += fscanf(ifp, "Passive earth-pressure coeff. (-) %lf\n", &kp);
  printf("%2d  kp              = %.2f\n", lest, kp);

  if (ifv < 2) {
    lest += fscanf(ifp, "#\nForest effects %[^\r\n]\n", foresteff);
    printf("%2d  foresteff       = %s\n", lest, foresteff);
  }
  if (!strncmp(foresteff, "no", 3) || !strncmp(foresteff, "No", 3))
    forest = 0;
  else if (!strncmp(foresteff, "yes", 4) || !strncmp(foresteff, "Yes", 4))
    forest = 1;
  else {
    printf("   read_command_file:  Invalid value of Forest effects -- %s.\n\n",
           foresteff);
    exit(15);
  }
  lest += fscanf(ifp, "Tree drag coefficient (-) %lf\n", &cD);
  printf("%2d  cD              = %4.2f\n", lest, cD);
  if (ifv < 2) {
    lest += fscanf(ifp, "Modulus of rupture (MPa) %lf\n", &MoR);
    printf("%2d  MoR             = %.1f\n", lest, MoR);
    lest += fscanf(ifp, "Forest decay coefficient (m/s) %lf\n", &decay_coeff);
    printf("%2d  tree_fail       = %4.2f\n", lest, decay_coeff);
  }

  if (ifv < 2)
    dummy += fscanf(ifp, "#\n");
  lest += fscanf(ifp, "Entrainment %[^\r\n]\n", erosion);
  printf("%2d  erosion         = %s\n", lest, erosion);
  if (!strncmp(erosion, "none", 5))
    eromod = 0;
  else if (!strncmp(erosion, "RAMMS", 6))
    eromod = 1;
  else if (!strncmp(erosion, "TJEM", 5) || !strncmp(erosion, "IsJo", 5))
    eromod = 2;
  else if (!strncmp(erosion, "AVAFRAME", 9) || !strncmp(erosion, "AvaFrame", 9))
    eromod = 3;
  else if (!strncmp(erosion, "GOEM", 5))
    eromod = 4;
  else {
    printf("   Entrainment model \"%s\" not implemented.\n", erosion);
    printf("   Calculation is carried out without entrainment.\n\n");
    eromod = 0;
  }
  printf("%2d  eromod          = %1d\n", lest, eromod);
  lest += fscanf(ifp, "Erosion coefficient (-) %lf\n", &k_erod);
  if (eromod == 0 || eromod == 2)
    k_erod = 0.0;                       /* No need for erosion coefficient */
  else if ((eromod == 1 || eromod == 3 || eromod == 4) && k_erod <= 0.0) {
    printf("   Warning:  You need k_erod > 0 to obtain erosion!\n");
    k_erod = 0.0;
  }
  printf("%2d  k_erod          = %5.3f\n", lest, k_erod);

  lest += fscanf(ifp, "Bed strength profile %[^\r\n]\n", gradient);
  printf("%2d  gradient        = %s\n", lest, gradient);
  if (!strncmp(gradient, "constant", 9))
    grad = 0;
  else if (!strncmp(gradient, "global", 7))
    grad = 1;
  else if (!strncmp(gradient, "local", 6))
    grad = 2;
  else {
    printf("   Bad property %s of bed strength profile. STOP!\n\n", gradient);
    exit(16);
  }
  if (grad < 2) {
    lest += fscanf(ifp, "Bed friction coefficient (-) %lf\n", &mu_s0);
    if (grad == 0 && mu_s0 != 0.0) {
      mu_s0 = 0.0;
      printf("%2d  mu_s0 set to 0.\n", lest);
    }
    else
      printf("%2d  mu_s0           = %5.3f\n", lest, mu_s0);
  }
  else {
    lest += fscanf(ifp, "Bed friction coefficient (-) %[^\r\n]\n", mu_s_fn);
    printf("%2d  mu_s_fn = %s\n", lest, mu_s_fn);
  }
  if (ifv > 0) {                        /* Versions 2020-06-23, 2021-10-25 */
    lest += fscanf(ifp, "Bed density (kg/m^3) %lf\n", &rho_b);
    printf("%2d  rho_b           = %.3f\n", lest, rho_b);
  }

  lest += fscanf(ifp, "Deposition %[^\r\n]\n", dep_flag);
  if (!strncmp(dep_flag, "no", 3))
    dep = 0;
  else if (!strncmp(dep_flag, "yes", 4)) {
    dep = 1;
    printf("   Warning:  Deposition not implemented in this version.\n");
  }
  else {
    printf("   Invalid value of dep_flag -- %s, dep set to 0.\n", dep_flag);
    dep = 0;
  }
  printf("%2d  dep             = %d\n", lest, dep);

  lest += fscanf(ifp, "Evolving geometry %[^\r\n]\n", dyn_surf_switch);
  if (!strncmp(dyn_surf_switch, "no", 3))
    dyn_surf = 0;
  else if (!strncmp(dyn_surf_switch, "yes", 4))
    dyn_surf = 1;
  if (dyn_surf == 1 && eromod == 0) {
    printf("%2d  No erosion, thus dyn_surf set to 0.\n", lest);
    dyn_surf = 0;
  }
  else
    printf("%2d  dyn_surf        = %d\n", lest, dyn_surf);

  if (ifv < 2)
    dummy += fscanf(ifp, "#\n# Numerical parameters\n");
  lest += fscanf(ifp, "#\nSimulation time (s) %lf\n", &t_max);
  printf("%2d  t_max           = %.2f\n", lest, t_max);
  lest += fscanf(ifp, "Minimum time step (s) %lf\n", &dt_min);
  printf("%2d  dt_min          = %6.4f\n", lest, dt_min);
  lest += fscanf(ifp, "Maximum time step (s) %lf\n", &dt_max);
  printf("%2d  dt_max          = %.4f\n", lest, dt_max);
  lest += fscanf(ifp, "Output interval (s) %lf\n", &dt_dump);
  printf("%2d  dt_dump         = %.2f\n", lest, dt_dump);
  lest += fscanf(ifp, "Write velocity vectors %[^\r\n]\n", write_vectors);
  printf("%2d  vectors         = %s\n", lest, write_vectors);
  lest += fscanf(ifp, "Write maximum pressure %[^\r\n]\n", write_max_press);
  printf("%2d  write_max_press = %s\n", lest, write_max_press);
  lest += fscanf(ifp, "Write instant. pressure %[^\r\n]\n", write_press);
  printf("%2d  write_press     = %s\n", lest, write_press);
  lest += fscanf(ifp, "Minimum flow depth (m) %lf\n", &h_min);
  printf("%2d  h_min           = %.4f\n", lest, h_min);
  lest += fscanf(ifp, "Minimum speed (m/s) %lf\n", &u_min);
  printf("%2d  u_min           = %.4f\n", lest, u_min);
  lest += fscanf(ifp, "Momentum threshold (kg m/s) %lf\n", &mom_thr);
  printf("%2d  mom_thr         = %.1f\n", lest, mom_thr);
  lest += fscanf(ifp, "Initial CFL number (-) %lf\n", &cfl);
  printf("%2d  CFL             = %5.3f\n", lest, cfl);

  n_items = (ifv == 0 ? 46 : (ifv == 1 ? 45 : 41));
  fclose(ifp);
  if (lest != n_items) {
     printf("   %d items read instead of %d. STOP!\n\n", lest, n_items);
     exit(13);
  }

  /* Check values for consistency and compute some derived quantities. */

  if (strlen(u_fn) > 0 || strlen(v_fn) > 0)
    restart = 1;                        /* File with start velocity present */

  utm_code = strtol(utm_str, &tail, 16);
  if (utm_code == 0)
    printf("   Invalid UTM zone. Set to 0.\n");
  else if (tail[0] != 'N' && tail[0] != 'n' \
           && tail[0] != 'S' && tail[0] != 's') {
    printf("   Invalid UTM hemisphere %s. Set UTM zone to 0.\n", tail);
    utm_code = 0;
  }
  else if (tail[0] == 'S' || tail[0] == 's') {
    utm_code *= (-1);
    printf("   UTM zone on southern hemisphere.\n");
  }
  if (epsg < 0 || epsg > 32767) {
    printf("   EPSG datum code outside allowed range, set to 0.\n");
    epsg = 0;
  }

  if (strlen(out_fn) == 0) {
    printf("   No output filename specified, write to ./result.\n");
    strncpy(out_fn, "result", 7);
  }

  /* Check densities: */
  if (rho <= 0.0) {
    printf("   Flow density has unphysical value %.1f kg/m^3. STOP!\n\n", rho);
    exit(17);
  }
  if (rho_b <= 0.0 || rho_d <= 0.0) {
    printf("   Bed and deposit densities > 0 are required. STOP!\n\n");
    exit(18);
  }
  rrb = rho / rho_b;
  rrd = rho / rho_d;
  /* The Grigorian–Ostroumov erosion model requires rrb, sigma > 1: */
  if (eromod == 4 && rrb <= 1.0) {
    printf("\n   Grigorian-Ostroumov model requires rho > rho_b. STOP!\n\n");
    exit(28);
  }
  sigma = 1.0 / sqrt(rrb-1.0);

  if (g < 0.0) {
    printf("   Gravity constant has unphysical value %.2f m/s^2. STOP!\n\n", g);
    exit(19);
  }

  if (strcmp(rheology, "Voellmy")) {
    printf("   Rheology \"%s\" not supported in this version. STOP!\n\n",
           rheology);
    exit(20);
  }
  if (strcmp(params, "constant") && strcmp(params, "variable")) {
    printf("   Parameters must be either \"constant\" or \"variable\".\n");
    printf("   Set to \"%s\". STOP!\n\n", params);
    exit(21);
  }
  if (!strcmp(params, "constant") && (mu_g < 0.0 || k_g < 0.0)) {
    printf("   mu >= 0 and k >= 0 required. STOP!\n\n");
    exit(22);
  }
  if (!strncmp(params, "variable", 9) && (strstr(mu_fn, "_mu.asc") == NULL)) {
    printf("   ! Dry-friction coefficient file does not end in '_mu.asc'.\n");
    printf("        %s\n", mu_fn);
  }
  if (!strncmp(params, "variable", 9) && (strstr(k_fn, "_k.asc") == NULL)) {
    printf("   ! Turb.-friction parameter file does not end in '_k.asc'.\n");
    printf("             %s\n", k_fn);
  }
  if (!strncmp(gradient, "local", 6) && (strstr(mu_s_fn, "_mu_s.asc")==NULL)) {
    printf("   ! Bed-friction coefficient file does not end in '_mu_s.asc'.\n");
    printf("        %s\n", mu_fn);
  }
  if (h_drag < 0.0)                     /* Effective flow depth in drag term */
    h_drag = 0.0;                       /* must be non-negative. */

  if (t_max < 0.0) {
    printf("   t_max must be >= 0.0, input as %.1f s. STOP!\n\n", t_max);
    exit(23);
  }
  if (dt_min <= 0.0 || dt_min > dt_max) {
    printf("   0.0 < dt_min <= dt_max required. STOP!\n\n");
    exit(24);
  }

  if (strncmp(write_vectors, "yes", 4) && strncmp(write_vectors, "no", 3)) {
    printf("   Value of \'Write velocity vectors\' must be \'yes\'");
    printf(" or \'no\'. STOP!\n\n");
    exit(25);
  }

  if (strncmp(write_max_press, "yes", 4) && strncmp(write_max_press, "no", 3)) {
    printf("   Value of \'Write maximum pressure\' must be \'yes\'");
    printf(" or \'no\'. STOP!\n\n");
    exit(25);
  }

  if (strncmp(write_press, "yes", 4) && strncmp(write_press, "no", 3)) {
    printf("   Value of \'Write pressure\' must be \'yes\'");
    printf(" or \'no\'. STOP!\n\n");
    exit(25);
  }

  if (h_min <= 0.0 || u_min <= 0.0 || mom_thr <= 0.0) {
    printf("   h_min, u_min, mom_thr > 0 required. STOP!\n\n");
    exit(26);
  }
  mom_thr /= rho;

  strcpy(max_fn, out_fn);
  if ((strrchr(max_fn, '.') == NULL && strlen(max_fn) > 252)
      || (strrchr(max_fn, '.') != NULL
          && (int) (strrchr(max_fn, '\0') - max_fn) > 251)) {
    printf("   Cannot create filename for max. values. STOP!\n\n");
    exit(27);
  }
  else if (strrchr(max_fn, '.') == NULL)
    strcat(max_fn, ".max");
  else
    strcpy(strrchr(max_fn, '.'), ".max");

  /* Create needed directories that do not exist yet. */

  /* Extract the name of the root folder for simulation results and create it
     if it does not exist yet: */
  strncpy(out_path, out_fn, 512);
  strncpy(out_path, dirname(out_path), 511);
  create_dir(out_path, "");

  /* If time slices are to be written, further subfolders are needed: */
  if (dt_dump < t_max) {
    create_dir(out_path, "h");          /* Folder for flow depth */
    create_dir(out_path, "s");          /* Folder for speed */
    if (!strncmp(write_vectors, "yes", 4)) {
      create_dir(out_path, "u");        /* Folders for velocity components */
      create_dir(out_path, "v");
    }
    if (!strncmp(write_press, "yes", 4))
      create_dir(out_path, "p");        /* Folder for stagnation pressure */
    if (dep > 0)
      create_dir(out_path, "d");        /* Folder for deposit depth */
    if (eromod > 0)
      create_dir(out_path, "b");        /* Folder for erodible snow depth */
    if (forest > 0)
      create_dir(out_path, "n");        /* Folder for stand density */
  }
  printf("\n   read_command_file:  Completed.\n");

}

/******************************/
/* End of read_command_file() */
/******************************/


/*************************/
/*                       */
/*  read_grid_file(...)  */
/*                       */
/*************************/

/** Reads the grid file and computes the needed geometrical information. */

void read_grid_file(char *grid_fn)

{
  FILE   *gfp;
  char   nodeline[256], xll[10], yll[10];
  int    lest = 0;
  short  shorty;
  float  floaty;
  double doubl, NaN;

  if ((gfp = fopen(grid_fn, "r")) == NULL) {
    printf("   read_grid_file:  Failed to open %s. STOP!\n\n", grid_fn);
    exit(30);
  }

  /* Read the header, which will be used in the output files: */
  lest += fscanf(gfp, "ncols %d\nnrows %d\n%s %lf\n%s %lf\n",
                 &m, &n, xll, &xllcorner, yll, &yllcorner);
  lest += fscanf(gfp, "cellsize %lf\nNODATA_value %lf\n", &cellsize, &NaN);
  if (lest != 8) {
    printf("\n   read_grid_file:  Incorrect grid file header. STOP!\n\n");
    exit(31);
  }

  if (!strcmp(xll, "xllcenter"))            /* If necessary, convert cell */
    xllcorner -= (0.5*cellsize);            /* center to cell corner      */
  if (!strcmp(yll, "yllcenter"))            /* coordinates */
    yllcorner -= (0.5*cellsize);

  if (!strncmp(fmt, "wb", 2)) {         /* BinaryTerrain */
    strncpy(header, "binterr1.3", 11);
    memcpy(header+10, &m, 4);               /* bth.xdim */
    memcpy(header+14, &n, 4);               /* bth.ydim */
    shorty = (short) 4;
    memcpy(header+18, &shorty, 2);          /* bth.datasize */
    shorty = (short) 1;
    memcpy(header+20, &shorty, 2);          /* bth.fp_flag  */
    memcpy(header+22, &shorty, 2);          /* bth.horiz_unit */
    shorty = (short) utm_code;
    memcpy(header+24, &shorty, 2);          /* bth.utm */
    shorty = (short) epsg;
    memcpy(header+26, &shorty, 2);          /* bth.epsg */
    memcpy(header+28, &xllcorner, 8);       /* bth.W_ext */
    doubl = (double) (xllcorner + m*cellsize);
    memcpy(header+36, &doubl, 8);           /* bth.E_ext */
    memcpy(header+44, &yllcorner, 8);       /* bth.S_ext */
    doubl = (double) (yllcorner + n*cellsize);
    memcpy(header+52, &doubl, 8);           /* bth.N_ext */
    shorty = (short) 0;
    memcpy(header+60, &shorty, 2);          /* bth.prj_flag */
    floaty = (float) 1.0;
    memcpy(header+62, &floaty, 4);          /* bth.scale */
    strncpy(header+66, "MoT-Voellmy "VERSION, 24);
    strncpy(header+150, " s", 3);           /* Time units */

    if (forest > 0)                     /* Forest needs special handling */
      memcpy(header_nD, header, 512);   /* on output */
  }
  else {                                /* ESRI ASCII Grid */
    sprintf(header, "ncols        %d\nnrows        %d\nxllcorner    %f\n",
            m, n, xllcorner);
    sprintf(nodeline, "yllcorner    %f\ncellsize     %f\nNODATA_value -9999\n",
            yllcorner, cellsize);
    strcat(header, nodeline);

    if (forest > 0)                     /* Forest needs special handling */
      strcat(header_nD, header);        /* on output */
  }

  fclose(gfp);

  /* Allocate all dynamic arrays and read the z-coordinates: */
  allocate();

  read_raster(grid_fn, z0, m, n, xllcorner, yllcorner, cellsize, -9998.9, 0);

  update_surface(z0);

  printf("   read_grid_file:     Completed.\n");

}

/***************************/
/* End of read_grid_file() */
/***************************/


/*************************/
/*                       */
/*  update_surface(...)  */
/*                       */
/*************************/

/** Compute slope and curvature components for given surface. */

void update_surface(double **Z)

{
  int i, j;
  double aux, auy, auz, dZdX, dZdY;
  double d2ZdX2, d2ZdY2, d2ZdXY;
  double cs2 = SQ(cellsize);
  double gsq;							/* For testing */

  /* Calculate the quantities that are used in the simulation: */
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      /* Slope angles and cell sizes */
      if (i == 0)
        dZdX = (Z[1][j] - Z[0][j]) / cellsize;    /* ∂Z/∂X */
      else if (i == m-1)
        dZdX = (Z[m-1][j] - Z[m-2][j]) / cellsize;
      else
        dZdX = 0.5 * (Z[i+1][j]-Z[i-1][j]) / cellsize;
      if (j == 0)
        dZdY = (Z[i][1] - Z[i][0]) / cellsize;    /* ∂Z/∂Y */
      else if (j == n-1)
        dZdY = (Z[i][n-1] - Z[i][n-2]) / cellsize;
      else
        dZdY = 0.5 * (Z[i][j+1] - Z[i][j-1]) / cellsize;
      aux = sqrt(1.0 + SQ(dZdX));
      auy = sqrt(1.0 + SQ(dZdY));
      auz = 1.0 + SQ(dZdX) + SQ(dZdY);
      dx[i][j] = cellsize * aux;
      dy[i][j] = cellsize * auy;
      gx[i][j] = -g * dZdX * aux / auz;
      gy[i][j] = -g * dZdY * auy / auz;
      auz = sqrt(auz);
      gz0[i][j] = g / auz;      /* Count normal component of gravity positive */
      gz[i][j] = gz0[i][j];     /* May be modified for centripetal acceler. */
      dA[i][j] = cs2 * auz;

      /* Off-diagonal element of metric tensor (normalized for use with
         physical instead of contravariant velocity components) – thus this
         is actually the cosine between x- and y-directions and not Gxy. */
      G_xy[i][j] = dZdX * dZdY / (aux*auy);

      /* Test whether this gives back the correct magnitude of g: */
      if (fabs((gsq = SQ(g/auz) + SQ(gx[i][j]) + SQ(gy[i][j])
                       + 2.0*G_xy[i][j]*gx[i][j]*gy[i][j]) - SQ(g)) > 0.0001)
        printf("   %3d, %3d:  |g| = %5.3f m/s²\n", i, j, sqrt(gsq)/g);

      /* Curvature tensor */
      if ((i == 0) || (i == m-1))
        d2ZdX2 = 0.0;                   /* ∂²Z/∂X² */
      else
        d2ZdX2 = (Z[i+1][j] + Z[i-1][j] - 2.0 * Z[i][j]) / SQ(cellsize);
      if ((j == 0) || (j == n-1))
        d2ZdY2 = 0;                     /* ∂²Z/∂Y² */
      else
        d2ZdY2 = (Z[i][j+1] + Z[i][j-1] - 2.0 * Z[i][j]) / SQ(cellsize);
      if ((i == 0) || (i == m-1) || (j == 0) || (j == n-1))
        d2ZdXY = 0.0;                   /* ∂²Z/∂X∂Y */
      else
        d2ZdXY = (Z[i+1][j+1] + Z[i-1][j-1] - Z[i+1][j-1] - Z[i-1][j+1])
                 / (4.0 * SQ(cellsize));
      IIxx[i][j] = d2ZdX2 / auz;
      IIyy[i][j] = d2ZdY2 / auz;
      IIxy[i][j] = d2ZdXY / auz;
    }
}

/***************************/
/* End of update_surface() */
/***************************/


/*************************/
/*                       */
/*  read_init_file(...)  */
/*                       */
/*************************/

/** Reads the initialization file(s) containing all the start values. */

void read_init_file(int m, int n)

{
  int    i, j, ec;
  double cs2 = SQ(cellsize), g_inv = 1.0 / g;


  /* Release area and release depth (compulsory file) */
  ec = read_raster(h_fn, h, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
  if (ec == 1) {
    printf("   read_init_file:     No file for release depth. STOP!\n");
    exit(40);
  }

  /* Initial velocities in x and y-direction (assume 0 if file not present). */
  ec =  read_raster(u_fn, u, m, n, xllcorner, yllcorner, cellsize, -9999, 1);
  if (ec > 0)
    printf("   read_init_file:     Could not read initial u velocity.\n");
  else if (ec < 0) {
    printf("   read_init_file:     Value out of bound in %s. STOP!\n", u_fn);
    exit(41);
  }
  ec += read_raster(v_fn, v, m, n, xllcorner, yllcorner, cellsize, -9999, 1);
  if (ec > 0)
    printf("   read_init_file:     Could not read initial v velocity.\n");
  else if (ec < 0) {
    printf("   read_init_file:     Value out of bound in %s. STOP!\n", v_fn);
    exit(41);
  }

  /* Erodible snow depth and erodibility */
  if (eromod > 0) {
    ec = read_raster(b_fn, b, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
    if (ec == 1) {
      printf("   read_init_file:     No file for erodible snow depth. STOP!\n");
      exit(42);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Value out of bound in %s. STOP!\n", b_fn);
      exit(41);
    }
    for (i = 0; i < m; i++)             /* No erodible snow in release area */
      for (j = 0; j < n; j++)
        if (h[i][j] > 0.0)
          b[i][j] = 0.0;
  }

  if (dep > 0) {
    for (i = i_min; i < i_max; i++)
      for (j = j_min; j < j_max; j++)
        d[i][j] = 0.0;          /* Initially, deposition is zero but previous
                                   avalanches could be read from file. */
  }

  if (eromod > 1) {             /* Read τ_c for all erosion models except RAMMS.
                                   AVAFRAME: this is specific erosion energy. */
    printf("   read_init_file:     About to read tau_c file...  ");
    ec = read_raster(tauc_fn, tau_c, m, n, xllcorner, yllcorner, \
                     cellsize, 0.0, 1);
    printf("done.\n");
    if (ec == 1) {
      printf("   read_init_file:     No file for bed shear strength. STOP!\n");
      exit(43);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Value out of bound in %s. STOP!\n",
             tauc_fn);
      exit(41);
    }
    for (i = 0; i < m; i++)             /* Scale tau_c by flow density and */
      for (j = 0; j < n; j++)           /* prevent too small values */
        tau_c[i][j] = MAX(tau_c[i][j]/rho, 0.1);    /* Units m²/s² */
    if (grad == 2) {                    /* Local bed friction angle from file */
      printf("   read_init_file:     About to read μ_s file...  ");
      ec = read_raster(mu_s_fn, mu_s, m, n, xllcorner, yllcorner, cellsize, \
                       0.0, 1);
      if (ec == 1) {
        printf("   read_init_file:     No file for bed friction coeff. STOP!\n");
        exit(44);
      }
      else if (ec < 0) {
        printf("   read_init_file:     Value out of bound in %s. STOP!\n", u_fn);
        exit(41);
      }
      printf("done.\n");
    }
    else {                             /* Constant friction angle of bed */
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
          mu_s[i][j] = mu_s0;
    }
  }

  /* Embed the release area in the snowcover if avalanche starts from rest: */
  if (!restart) {
    if (eromod == 0) {                  /* Subtract h0 from z0 */
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
          z0[i][j] -= (h[i][j] * gz[i][j] * g_inv);
    }
    else {                              /* Add b0 to z0, subtract h0 */
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
          z0[i][j] += ((b[i][j] - h[i][j]) * gz[i][j] * g_inv);
    }
    update_surface(z0);                 /* Compute new slope/curvature */
  }

  /* Friction parameters */
  if (!strcmp(params, "variable")) {    /* Variable coefficients from file */
    ec = read_raster(mu_fn, mu, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
    if (ec > 0) {
      printf("   read_init_file:     Missing file for mu. STOP!\n");
      exit(45);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Out-of-bound value in %s. STOP!\n", mu_fn);
      exit(41);
    }
    ec += read_raster(k_fn, k, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
    if (ec > 0) {
      printf("   read_init_file:     Missing file for k. STOP!\n");
      exit(46);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Out-of-bound value in %s. STOP!\n", k_fn);
      exit(41);
    }
  }
  else {                                /* Constant friction parameters */
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) {
        mu[i][j] = mu_g;
        k[i][j] = k_g;
      }
  }

  /* Forest parameters: tree density times diameter, diameter.
     Note that the input file refers to horizontally projected area. */
  if (forest > 0) {
    ec = read_raster(nD_fn, nD, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
    if (ec > 0) {
      printf("   read_init_file:     Missing file for nD. STOP!\n");
      exit(47);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Out-of-bound value in %s. STOP!\n", nD_fn);
      exit(41);
    }
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        nD[i][j] *= (cs2 / dA[i][j]);
    ec = read_raster(tD_fn, tD, m, n, xllcorner, yllcorner, cellsize, 0.0, 1);
    if (ec > 0) {
      printf("   read_init_file:     Missing file for tD. STOP!\n");
      exit(47);
    }
    else if (ec < 0) {
      printf("   read_init_file:     Out-of-bound value in %s. STOP!\n", tD_fn);
      exit(41);
    }
  }

  /* Impenetrable areas */
  /* Not implemented yet! */

  /* Initialize the conservative field variables: */
  mov_vol = 0.0;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      f_new[i][j][0] = h[i][j] * dA[i][j];
      f_new[i][j][1] = f_new[i][j][0] * u[i][j];
      f_new[i][j][2] = f_new[i][j][0] * v[i][j];
      if (eromod > 0)
        b_min[i][j] = b[i][j];
      mov_vol += f_new[i][j][0];
      /* Initialize erosion rate to 0 here so that it need not be computed
         again in each timestep when running without erosion. */
      src[i][j][0] = 0.0;
      if (dep > 0)
        d_max[i][j] = 0.0;
    }
  }

  printf("   read_init_file:     Completed.\n\n");
}

/***************************/
/* End of read_init_file() */
/***************************/


/********************/
/*                  */
/* read_raster(...) */
/*                  */
/********************/

/** Opens an AAIGrid raster file and reads its content into an array, doing
    some consistency checks and eliminating values below a threshold. */

int read_raster(char *raster_fn, double** X, int m, int n, double xll,
                double yll, double cs, double min_val, int pass)
{
  FILE   *ifp;
  char   a[10], b[10];                  /* Check xllcorner or xllcenter? */
  int    mr, nr, i, j;
  double xll_read, yll_read, cs_read, nan, fval;

  /* Open raster file for reading. */
  if ((ifp = fopen(raster_fn, "r")) == NULL) {
    printf("   read_raster:        Could not open file %s.\n", raster_fn);
    return 1;
  }

  /* Read file header and check values for consistency. */
  if (fscanf(ifp, "ncols %d nrows %d %s %lf %s %lf \
             cellsize %lf NODATA_value %lf\n",
             &mr, &nr, a, &xll_read, b, &yll_read, &cs_read, &nan) != 8) {
    printf("   Error reading header of file %s. STOP!\n\n", raster_fn);
    exit(50);
  }
  if (pass > 0 && (mr != m || nr != n || fabs(cs_read - cs) > 0.0001
                   || fabs(xll_read - xll) > 0.001
                   || fabs(yll_read - yll) > 0.001)) {
    printf("   read_raster:  Contradiction in header data of %s. STOP!\n",
           raster_fn);
    printf("      m = %d, mr = %d;  n = %d, nr = %d;  cs = %.3f, csr = %.3f\n",
           m, mr, n, nr, cs, cs_read);
    exit(51);
  }

  /* Read data one by one. */
  for (j = n-1; j >= 0; j--) {
    for (i = 0; i < m; i++) {
      /* lest = fscanf(ifp, "%lf", &fval);
      if (lest == 0) { */
      if (fscanf(ifp, "%lf", &fval) != 1) {
        printf("   Error reading data from file %s at (%d,%d). STOP!\n\n",
               raster_fn, i, j);
        exit(52);
      }
      if (fval >= min_val)
        X[i][j] = fval;
      else {
        printf("   read_raster:  Reading %s.\n", raster_fn);
        printf("                 Value at (%d,%d) is %.5f < %.5f. STOP!\n",
               i, j, fval, min_val);
        exit(53);
      }
    }
  }
  fclose(ifp);
  return 0;
}

/***************************/
/* End of read_raster(...) */
/***************************/


/*********************/
/*                   */
/*  write_data(...)  */
/*                   */
/*********************/

/* This function prepares the header for the output files and calls
   writeout() repeatedly to produce the result files for timeslices and
   for the maximum values at the end of the simulation.
   To facilitate writing the deposit depth at the end (pass == 2), a separate
   argument **h_dep is required, which can be a dummy field for pass == 1
   since it will not be used. */

void write_data(double tid, double **h_dep, double **h, double **b, double **d,
              double **s, double **u, double **v, double **p_imp, double **nD,
              int imin, int imax, int jmin, int jmax, int pass, char *fmt)

{
  int    i, j, di, dj;
  float  tempus;
  double westend, eastend, southend, northend;
  char   line[256], suf[10];    /* Temporary variable for ASCII format */
  time_t now;                   /* Date and time of run */

  di = imax - imin;             /* Number of cells in x-direction */
  dj = jmax - jmin;             /* Number of cells in y-direction */
  if (di < 1 || dj < 1) {
    printf("   write_data:  Nothing to print.\n");
    return;
  }
  westend  = xllcorner + imin * cellsize;
  eastend  = xllcorner + (imin+di) * cellsize;
  southend = yllcorner + jmin * cellsize;
  northend = yllcorner + (jmin+dj) * cellsize;

  tempus = (float) tid;

  /* Prepare the header */

  if (!strncmp(fmt, "wb", 2)) {         /* BinaryTerrain format */
    memcpy(header+146, &tempus, sizeof(float));
    memcpy(header+10, &di, 4);              /* bth.xdim */
    memcpy(header+14, &dj, 4);              /* bth.ydim */
    memcpy(header+28, &westend, 8);         /* bth.W_ext */
    memcpy(header+36, &eastend, 8);         /* bth.E_ext */
    memcpy(header+44, &southend, 8);        /* bth.S_ext */
    memcpy(header+52, &northend, 8);        /* bth.N_ext */
    now = time(NULL);
    strftime(header+90, 24, "%Y-%m-%d %H:%M:%S %z", localtime(&now));
  }
  else {                                /* ESRI ASCII Grid format */
    sprintf(header, "ncols        %d\nnrows        %d\nxllcorner    %.1f\n",
            di, dj, westend);
    sprintf(line, "yllcorner    %.1f\ncellsize     %.2f\n",
            southend, cellsize);
    strncat(header, line, 256);
    sprintf(line, "NODATA_value -9999\n");
    strncat(header, line, 256);
  }

  /* Call writeout repeatedly to write the files */

  if (pass == 1) {                  /* Write timeslice of the fields */
    printf("   write_data:  Output %04d at time %7.2f...   ", n_dump, tempus);
    /* Flow depth */
    sprintf(suf, "_h_%04d", n_dump);
    writeout(h, suf, fmt, imin, imax, jmin, jmax, header,
             "h -- Flow depth (m)            ", "5.2");
    /* Speed */
    sprintf(suf, "_s_%04d", n_dump);
    writeout(s, suf, fmt, imin, imax, jmin, jmax, header,
             "s -- Flow speed (m/s)          ", "6.2");
    /* Write snow cover depth only if erosion was specified */
    if (eromod > 0) {
      sprintf(suf, "_b_%04d", n_dump);
      writeout(b, suf, fmt, imin, imax, jmin, jmax, header,
               "b -- Erodible snow depth (m)   ", "5.3");
    }
    /* Write deposition only if it is activated */
    if (dep > 0) {
      sprintf(suf, "_d_%04d", n_dump);
      writeout(d, suf, fmt, 0, m, 0, n, header,
               "d -- Deposit depth (m)         ", "6.3");
    }
    /* Write files for u, v only if requested */
    if (!strncmp(write_vectors, "yes", 4)) {
      sprintf(suf, "_u_%04d", n_dump);
      writeout(u, suf, fmt, imin, imax, jmin, jmax, header,
               "u -- x-velocity (m/s)          ", "7.2");
      sprintf(suf, "_v_%04d", n_dump);
      writeout(v, suf, fmt, imin, imax, jmin, jmax, header,
               "v -- y-velocity (m/s)          ", "7.2");
    }
    /* Write pressure only if requested */
    if (!strncmp(write_press, "yes", 4)) {
      sprintf(suf, "_p_%04d", n_dump);
      writeout(p_imp, suf, fmt, imin, imax, jmin, jmax, header,
               "p -- impact pressure (kPa)     ", "7.2");
    }
    /* Write forest density nD only if forest effects are included.
       Need to write nD over entire DEM area to see remaining forest. */
    if (forest > 0) {
      sprintf(suf, "_n_%04d", n_dump);
      writeout(nD, suf, fmt, 0, m, 0, n, header_nD,
               "nD -- braking effect (1/m)     ", "7.4");
    }
  }

  else if (pass == 2) {             /* Write fields of maximum values */
    printf("   write_data:  Write maximum values of fields...");
    if (!strncmp(fmt, "wb", 2)) {           /* BinaryTerrain format */
      /* Set simulation time to INFINITY: */
      tempus = (float) INFINITY;
      memcpy(header+146, &tempus, sizeof(float));
    }
    /* Deposit depth */
    if (dep == 0) {
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
          h_dep[i][j] = rrd * h[i][j];
    }
    writeout(h_dep, "_h_dep", fmt, 0, m, 0, n, header,
             "h_dep -- Deposit depth (m)     ", "5.2");
    /* Maximum flow depth */
    writeout(h_max, "_h_max", fmt, 0, m, 0, n, header,
             "h_max -- Max. flow depth (m)   ", "5.2");
    /* Maximum speed */
    writeout(s_max, "_s_max", fmt, 0, m, 0, n, header,
             "s_max -- Max. speed (m/s)      ", "6.2");
    /* Write min. snow cover depth only if erosion was specified: */
    if (eromod > 0)
      writeout(b_min, "_b_min", fmt, 0, m, 0, n, header,
               "b -- Min. snowpack depth (m)   ", "5.3");
    /* Write deposition only if it is activated */
    if (dep > 0)
      writeout(d_max, "_d_max", fmt, 0, m, 0, n, header,
               "d_max -- Max. deposit (m)      ", "6.3");
    /* Write files for u, v only if requested: */
    if (!strncmp(write_vectors, "yes", 4)) {
      writeout(u_max, "_u_max", fmt, 0, m, 0, n, header,
               "u_max -- Max. x-velocity (m/s) ", "7.2");
      writeout(v_max, "_v_max", fmt, 0, m, 0, n, header,
               "v_max -- Max. y-velocity (m/s) ", "7.2");
    }
    /* Write maximum pressure only if requested: */
    if (!strncmp(write_max_press, "yes", 4))
      writeout(p_max, "_p_max", fmt, 0, m, 0, n, header,
               "p_max -- Max. pressure (kPa)   ", "7.2");
    /* Write forest density nD only if forest can be destroyed */
    if (forest > 0)
      writeout(nD, "_nD_min", fmt, 0, m, 0, n, header_nD,
               "nD_min -- braking effect (1/m) ", "7.4");
  }

  printf(" done.\n");
}

/***********************/
/* End of write_data() */
/***********************/


/*******************/
/*                 */
/*  writeout(...)  */
/*                 */
/*******************/

/** Called by write_data repeatedly to create output files in AAIGrid or
    Binaryterrain 1.3 format either for time slice or max. values at the
    end of a run. */

void writeout(double **F, char* suffix, char* fmt, int imin, int imax,
              int jmin, int jmax, char* header, char* descr, char* ascfmt)

{
  int    i, j, l;
  size_t nitems, length;
  char   fn[1024], bn[511], dn[511], *addr, fmtstr[10], fmtstrn[11];
  FILE * ofp;                           /* Pointer to output file handle */

  /* For time slices, reconstruct the directory where to write the file:
     All files for maximum (minimum) values go into the folder contained in
     out_fn while time slices of fields are written to the subfolder reserved
     for the specific field. */

  strncpy(fn, out_fn, 510);
  if (strncmp(suffix+3, "m", 1) && strncmp(suffix+4, "m", 1)
      && strncmp(suffix+3, "dep", 3)) {
    strncpy(bn, basename(fn), 510);     /* Time-slice files (no 'max'/'min') */
    strncpy(dn, dirname(fn), 510);      /* go into specific subfolders. */
    sprintf(fn, "%s%s%c%s%s", dn, DIRSEP, suffix[1], DIRSEP, bn);
  }
  strncat(fn, suffix, 7);               /* Common to time-slice and max files */
  if (!strncmp(fmt, "wb", 2))           /* BinaryTerrain format */
    strncat(fn, ".bt", 4);
  else                                  /* ESRI ASCII Grid format */
    strncat(fn, ".asc", 5);

  if ((ofp = fopen(fn, fmt)) == NULL) {
    printf("\n   writeout:  Failed to open output file %s. STOP!\n\n", fn);
    exit(60);
  }

  /* Construct the file headers and afterwards write the data: */

  if (!strncmp(fmt, "wb", 2)) {           /* BinaryTerrain format */
    strncpy(header+114, descr, 32);
    /* When inserting basename into the header, catch file names w/out '/'!
       Also, if basename is longer than 103 chars and does not fit into the
       header, replace the basename by 'TRUNCATED'. */
    length = ((addr = strrchr(fn, DIRSEP[0])) == NULL ?
               strlen(fn) :                 /* No directory separator in fn. */
               strlen(addr+1));             /* There is dir. separator in fn. */
    if (length < 104)
      strncpy(header+152, addr+1, length+1);
    else
      strcpy(header+152, "TRUNCATED");
    if (fwrite(header, 1, 256, ofp) != 256) {
      printf("\n   writeout:  Could not write file header. STOP!\n\n");
      exit(61);
    }
    for (i = imin, l = 0; i < imax; i++)
      for (j = jmin; j < jmax; j++, l++)
        *(data+l) = (float) F[i][j];
    nitems = (size_t) ((imax-imin) * (jmax-jmin));
    if (fwrite(data, sizeof(float), nitems, ofp) != nitems) {
      printf("\n   writeout:  Failed to write data to file. STOP!\n\n");
      exit(62);
    }
  }

  else {                                /* ESRI ASCII Grid format */
    sprintf(fmtstr,  "\"%%%s \"",   ascfmt);
    sprintf(fmtstrn, "\"%%%s\\n\"", ascfmt);
    if (fprintf(ofp, "%s", header) < 0) {
      printf("\n   writeout:  Could not write file header. STOP!\n\n");
      exit(61);
    }
    for (j = jmax-1; j >= jmin; j--) {
      for (i = imin; i < imax-1; i++)
        fprintf(ofp, "%.3f ", F[i][j]);
      fprintf(ofp, "%.3f\n", F[imax-1][j]);
    }
  }

  fclose(ofp);
}

/*********************/
/* End of writeout() */
/*********************/


/********************/
/*                  */
/*  allocate2(...)  */
/*                  */
/********************/

/*  Allocation of a two-dimensional double array using pointers to pointers. */

double **allocate2(int m, int n)


{
  int i, tries;
  double **p;

  /* Allocate array of pointers to 1D subarrays first */
  tries = 0;
  while (tries < TRIES_MAX
         && (p = (double**) malloc(m * sizeof(double*))) == NULL) {
    tries++;
    sleep(TRY_WAIT);
  }
  if (tries >= TRIES_MAX) {
    printf("   allocate2:  Memory allocation failed. STOP!\n\n");
    exit(6);
  }

  /* Now allocate space for each of the 1D subarrays  */
  for (i = 0; i < m; i++) {
    tries = 0;
    while (tries < TRIES_MAX
           && (p[i] = (double*) malloc(n * sizeof(double))) == NULL) {
      tries++;
      sleep(TRY_WAIT);
    }
    if (tries >= TRIES_MAX) {
      printf("   allocate2:  Memory allocation failed. STOP!\n\n");
      exit(6);
    }
  }

  return p;
}

/***************************/
/*  End of allocate2(...)  */
/***************************/


/********************/
/*                  */
/*  allocate3(...)  */
/*                  */
/********************/

/* Allocation of a three-dimensional double array using pointers to pointers
   to pointers. */

double*** allocate3(int m, int n, int k)

{
  int i, j, tries;
  double ***p;

  /* Allocate array of pointers to 2D subarrays first */
  tries = 0;
  while (tries < TRIES_MAX
         && (p = (double***) malloc(m * sizeof(double**))) == NULL) {
    tries++;
    sleep(TRY_WAIT);
  }
  if (tries >= TRIES_MAX) {
    printf("   allocate3:  Memory allocation failed. STOP!\n\n");
    exit(7);
  }

  /* Now allocate space for each of the 1D subarrays  */
  for (i = 0; i < m; i++) {
    tries = 0;
    while (tries < TRIES_MAX
           && (p[i] = (double**) malloc(n * sizeof(double*))) == NULL) {
      tries++;
      sleep(TRY_WAIT);
    }
    if (tries >= TRIES_MAX) {
      printf("   allocate3:  Memory allocation failed. STOP!\n\n");
      exit(7);
    }

    /* Finally allocate the space for the m×n 1D subarrays themselves: */
    for (j = 0; j < n; j++) {
      tries = 0;
      while (tries < TRIES_MAX
             && (p[i][j] = (double*) malloc(k * sizeof(double))) == NULL) {
        tries++;
        sleep(TRY_WAIT);
      }
      if (tries >= TRIES_MAX) {
        printf("   allocate3:  Memory allocation failed. STOP!\n\n");
        exit(7);
      }
    }
  }

  return p;
}

/***************************/
/*  End of allocate3(...)  */
/***************************/


/**********************/
/*                    */
/*  deallocate2(...)  */
/*                    */
/**********************/

/* Frees the storage space occupied by a two-dimensional array p[][] and
   by the pointers to it. */

void deallocate2(double **p, int m)

{
  int i;

  for (i = 0; i < m; i++)
    free(p[i]);
  free(p);
  p = NULL;
}

/*****************************/
/*  End of deallocate2(...)  */
/*****************************/


/**********************/
/*                    */
/*  deallocate3(...)  */
/*                    */
/**********************/

/* Frees the storage space occupied by a three-dimensional array p[][][] and
   by the pointers to it. */

void deallocate3(double ***p, int m, int n)

{
  int i, j;

  for (i = m-1; i >= 0; i--) {
    for (j = n-1; j >= 0; j--)
      free(p[i][j]);
    free(p[i]);
  }
  free(p);
  p = NULL;
}

/*****************************/
/*  End of deallocate3(...)  */
/*****************************/


/****************/
/*              */
/*  allocate()  */
/*              */
/****************/

/* Allocates all the two- and three-dimensional arrays used in the program. */

void allocate(void)

{
  int     tries;                /* # failed memory allocation attempts */

  f_old   = allocate3(m, n, 3);
  f_new   = allocate3(m, n, 3);
  src     = allocate3(m, n, 3);

  dx      = allocate2(m, n);
  dy      = allocate2(m, n);
  dA      = allocate2(m, n);
  gx      = allocate2(m, n);
  gy      = allocate2(m, n);
  gz      = allocate2(m, n);
  gz0     = allocate2(m, n);
  G_xy    = allocate2(m, n);
  IIxx    = allocate2(m, n);
  IIyy    = allocate2(m, n);
  IIxy    = allocate2(m, n);
  h       = allocate2(m, n);
  s       = allocate2(m, n);
  u       = allocate2(m, n);
  v       = allocate2(m, n);
  d       = allocate2(m, n);
  p_imp   = allocate2(m, n);
  h_max   = allocate2(m, n);
  s_max   = allocate2(m, n);
  u_max   = allocate2(m, n);
  v_max   = allocate2(m, n);
  p_max   = allocate2(m, n);
  mu      = allocate2(m, n);
  k       = allocate2(m, n);
  z0      = allocate2(m, n);

  if (!strncmp(fmt, "wb", 2)) {
    tries = 0;
    while (tries < TRIES_MAX
           && (data = (float*) malloc(m*n * sizeof(float))) == NULL) {
      tries++;
      sleep(TRY_WAIT);
    }
    if (tries >= TRIES_MAX) {
      printf("   allocate:  Memory allocation failed. STOP!\n\n");
      exit(8);
    }
  }

  if (eromod > 0) {                     /* All erosion models */
    b     = allocate2(m, n);
    b_min = allocate2(m, n);
    if (eromod > 1) {                   /* TJEM, AvaFrame or GOEM */
      tau_c = allocate2(m, n);
      mu_s  = allocate2(m, n);
    }
  }

  if (forest > 0) {                     /* Account for braking by forest */
    nD = allocate2(m, n);
    tD = allocate2(m, n);
    decay_const = allocate2(m, n);
  }

  if (dep > 0) {
    d     = allocate2(m, n);
    d_max = allocate2(m, n);
  }

  if (dyn_surf > 0)
    z      = allocate2(m, n);
}

/***********************/
/*  End of allocate()  */
/***********************/


/******************/
/*                */
/*  deallocate()  */
/*                */
/******************/

/* Frees all the dynamically allocated arrays used in the program. */

void deallocate(void)

{
  if (dyn_surf > 0)
    deallocate2(z, m);

  if (dep > 0) {
    deallocate2(d_max, m);
    deallocate2(d, m);
  }

  if (forest > 0) {                     /* Account for braking by forest */
    deallocate2(decay_const, m);
    deallocate2(tD, m);
    deallocate2(nD, m);
  }

  if (eromod > 0) {
    if (eromod > 1) {                   /* TJEM, AvaFrame, GOEM only */
      deallocate2(mu_s, m);
      deallocate2(tau_c, m);
    }
    deallocate2(b_min, m);              /* All erosion models */
    deallocate2(b, m);
  }

  if (!strncmp(fmt, "wb", 2))
    free(data);

  deallocate3(src, m, n);
  deallocate3(f_new, m, n);
  deallocate3(f_old, m, n);

  deallocate2(z0, m);
  deallocate2(k, m);
  deallocate2(mu, m);
  deallocate2(p_max, m);
  deallocate2(v_max, m);
  deallocate2(u_max, m);
  deallocate2(s_max, m);
  deallocate2(h_max, m);
  deallocate2(p_imp, m);
  deallocate2(d, m);
  deallocate2(v, m);
  deallocate2(u, m);
  deallocate2(s, m);
  deallocate2(h, m);
  deallocate2(IIxy, m);
  deallocate2(IIyy, m);
  deallocate2(IIxx, m);
  deallocate2(G_xy, m);
  deallocate2(gz0, m);
  deallocate2(gz, m);
  deallocate2(gy, m);
  deallocate2(gx, m);
  deallocate2(dx, m);
  deallocate2(dy, m);
  deallocate2(dA, m);
}

/*************************/
/*  End of deallocate()  */
/*************************/


/******************/
/*                */
/*  create_dir()  */
/*                */
/******************/

/** Checks whether the folder formed from the two function arguments exists
    and creates it if not. */

void create_dir(char *main_folder, char *subfolder)

{
  struct stat sb;               /* Structure with file information */
  char   folder[1024];
  int    ec;                    /* Error code */

  sb.st_mode = 0;
  sprintf(folder, "%s%s%s", main_folder, DIRSEP, subfolder);

  /* Check whether target folder exists: */
  stat(folder, &sb);

  /* If the target folder is not found, try to create it: */
  if (S_ISDIR(sb.st_mode) != 1) {
    #ifdef LINUX
      ec = mkdir(folder, 0755);
    #endif
    #ifdef WINDOWS
      ec = mkdir(folder);
    #endif
    if (ec != 0) {              /* Some error has occurred, exit. */
      printf("\n   Failed to create missing target folder %s. STOP!\n", folder);
      printf("   mkdir error code:  %d.\n\n", ec);
      exit(70);
    }
  }
}

/*************************/
/*  End of create_dir()  */
/*************************/

