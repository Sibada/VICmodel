/******************************************************************************
 * @section MODIFICATION
 *
 * Modification by Ruida Zhong for the R package VICmodel on Jun 23th, 2018:
 * The stdout/stderr are commented and  `fprintf` are changed to `Rprintf`
 * for the correct output to R terminal.
 * The parts of "Report Bugs and Issues to ..." in print_licence() are commented
 * because this is an R package, to avoid confusion and misunderstanding for the
 * users and trouble and disturbance for the VIC source code developers.
 *
 * @section DESCRIPTION
 *
 * This routine checks the command line for valid program options.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <unistd.h>
#include <getopt.h>

#include <vic_driver_shared_all.h>

char *optstring = "g:vo";

/**********************************************************************
   cmd_proc                  Keith Cherkauer                1997

   This routine checks the command line for valid program options.  If
   no options are found, or an invalid combination of them appear, the
   routine calls usage() to print the model usage to the screen, before
   exiting execution.
**********************************************************************/
void
cmd_proc(int    argc,
         char **argv,
         char  *globalfilename)
{
    /*
    char GLOBAL_SET;
    int  optchar;

    if (argc == 1) {
        print_usage(argv[0]);
        exit(1);
    }

    GLOBAL_SET = false;

    while ((optchar = getopt(argc, argv, optstring)) != EOF) {
        switch ((char)optchar) {
        case 'v':
            display_current_settings(DISP_VERSION);
            exit(EXIT_SUCCESS);
            break;
        case 'o':
            display_current_settings(DISP_COMPILE_TIME);
            exit(EXIT_SUCCESS);
            break;
        case 'g':
            strncpy(globalfilename, optarg, MAXSTRING);
            GLOBAL_SET = true;
            break;
        default:
            print_usage(argv[0]);
            exit(1);
            break;
        }
    }

    if (!GLOBAL_SET) {
        fprintf(stderr,
                "ERROR: Must set global control file using the '-g' flag\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
     */
}

/******************************************************************************
 * @brief    This routine prints out usage details.
 *****************************************************************************/
void
print_usage(char *executable)
{
    Rprintf(
            "Usage: %s [-v | -o | -g <global_parameter_file>]\n", executable);
    Rprintf("  v: display version information\n");
    Rprintf(
            "  o: display compile-time options settings"
            " (set in .h files)\n");
    Rprintf(
            "  g: read model parameters from <global_parameter_file>.\n");
    Rprintf(
            "       <global_parameter_file> is a file that contains all"
            " needed model\n");
    Rprintf(
            "       parameters as well as model option flags, and the names"
            " and\n");
    Rprintf("       locations of all other files.\n");
}

/******************************************************************************
 * @brief    This routine prints out the model Version
 *****************************************************************************/
void
print_version(char *driver)
{
    Rprintf("VIC Driver  : %s\n", driver);
    Rprintf("VIC Version : %s\n", VERSION);
    Rprintf("VIC Git Tag : %s\n", GIT_VERSION);
    Rprintf("Compiled    : by %s on %s (%s) %s %s\n",
            USERNAME, HOSTNAME, PLATFORM, BUILD_DATE, BUILD_TIME);
    Rprintf("Compiler    : %s\n", COMPILER);
    Rprintf(" version    : %s\n", COMPILER_VERSION);

    print_license();
}

/******************************************************************************
 * @brief    This routine prints out license information
 *****************************************************************************/
void
print_license()
{
    Rprintf(
            "\n  Variable Infiltration Capacity (VIC) macroscale hydrologic\n");
    Rprintf(
            "  model version %s, Copyright (C) 2016 Computational\n",
            SHORT_VERSION);
    Rprintf(
            "  Hydrology Group, Dept. of Civil and Environmental Engineering,\n");
    Rprintf(
            "  University of Washington.  VIC comes with ABSOLUTELY NO\n");
    Rprintf(
            "  WARRANTY. This is free software, you may redistribute it\n");
    Rprintf(
            "  under certain conditions; see LICENSE for details.\n\n");

    //Rprintf(
    //        "  Report Bugs and Issues to : https://github.com/UW-Hydro/VIC/issues\n");
    //Rprintf(
    //        "  VIC Users Email Listserve : vic_users@u.washington.edu \n\n");
}
