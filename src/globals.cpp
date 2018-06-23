#include <vic_R.h>

// global variables
int                 flag;
size_t              NR; /* array index for atmos struct that indicates
                           the model step avarage or sum */
size_t              NF; /* array index loop counter limit for atmos
                           struct that indicates the SNOW_STEP values */
FILE *LOG_DEST;
char vic_run_ref_str[MAXSTRING];

double (*funcd)(double z, double es, double Wind, double AirDens, double ZO,
        double EactAir, double F, double hsalt, double phi_r,
        double ushear,
        double Zrh);

global_param_struct global_param;
option_struct       options;
parameters_struct   param;
param_set_struct    param_set;
metadata_struct     out_metadata[N_OUTVAR_TYPES];
