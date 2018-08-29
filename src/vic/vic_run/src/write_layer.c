/******************************************************************************
*
* @section MODIFICATION
*
* Modification by Ruida Zhong for the R package VICmodel on June 24th, 2018:
* The `Rprintf` are modified to `RRprintf` for It can correctly output to the
* R console.
*
* @section DESCRIPTION
*
* This routine writes soil variables to stdout.
*
* It creates a table of soil moisture values which shows how much liquid water
* and ice are contained in the thawed, frozen and unfrozen sublayers of each
* soil layer.  It also gives the total soil moisture for each layer.
*
* @section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014  The Land Surface Hydrology Group, Department of Civil
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
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine writes soil variables to stdout.
******************************************************************************/
void
write_layer(layer_data_struct *layer,
            int                veg,
            double            *frost_fract)
{
    extern option_struct options;

    size_t               index;
    double               layer_moist;
    double               sum_moist;
    size_t               frost_area;
    double               avg_ice;

    Rprintf("Layer Data for Vegetation Type #%i\n", veg);
    Rprintf("Layer:\t");
    for (index = 0; index < options.Nlayer; index++) {
        Rprintf("\t\t%zu", index + 1);
    }
    Rprintf("\nEvaporation:\t");
    for (index = 0; index < options.Nlayer; index++) {
        Rprintf("\t%f", layer[index].evap);
    }
    Rprintf("\n      Kappa:\t");
    for (index = 0; index < options.Nlayer; index++) {
        Rprintf("\t%f", layer[index].kappa);
    }
    Rprintf("\n         Cs:\t");
    for (index = 0; index < options.Nlayer; index++) {
        Rprintf("\t%f", layer[index].Cs);
    }
    Rprintf(
        "\n\nMoisture Table\n---------------------------------------------------------------------------\n Moist:\t");
    for (index = 0; index < options.Nlayer; index++) {
        Rprintf("\t%f", layer[index].moist);
    }
    Rprintf("\n        Ice:\t");
    for (index = 0; index < options.Nlayer; index++) {
        avg_ice = 0;
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            avg_ice += layer[index].ice[frost_area] * frost_fract[frost_area];
        }
        Rprintf("\t%f", avg_ice);
    }
    Rprintf(
        "\n---------------------------------------------------------------------------\nLayer Moist:\t");
    sum_moist = 0.;
    for (index = 0; index < options.Nlayer; index++) {
        layer_moist = layer[index].moist;
        sum_moist += layer_moist;
        Rprintf("\t%f", layer_moist);
    }
    Rprintf("\n\n-----> Total Moisture = %f\n\n", sum_moist);
}
