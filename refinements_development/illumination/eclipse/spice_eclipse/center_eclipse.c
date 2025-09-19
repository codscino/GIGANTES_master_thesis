#include <stdio.h>
#include <string.h>
#include "SpiceUsr.h"

int main()
{
    /*
    Local constants
    */
    #define MAXWIN       1000
    #define TIMFMT       "YYYY-MON-DD HR:MN:SC.###### TDB::RND"
    #define TIMLEN       41

    /*
    Local variables
    */
    SPICEDOUBLE_CELL ( cnfine, MAXWIN );
    SPICEDOUBLE_CELL ( result, MAXWIN );

    // Body, frame, and shape definitions for the occultation search
    ConstSpiceChar * occtyp      = "ANY";          // Find any eclipse (penumbral or umbral)
    ConstSpiceChar * front_body  = "SATURN";       // The body casting the shadow
    ConstSpiceChar * front_shape = "ELLIPSOID";
    ConstSpiceChar * front_frame = "IAU_SATURN";
    ConstSpiceChar * back_body   = "SUN";          // The light source being blocked
    ConstSpiceChar * back_shape  = "ELLIPSOID";
    ConstSpiceChar * back_frame  = "IAU_SUN";
    ConstSpiceChar * observer    = "ENCELADUS";    // The object being eclipsed
    ConstSpiceChar * abcorr      = "CN+S";         // Converged Newtonian + Stellar Aberration

    // GF search parameters
    SpiceDouble         step        = 60.0; // Step size in seconds. 1 minute is fine for this.

    // Time window variables
    SpiceDouble         et_center, et_start, et_end;
    SpiceDouble         start_time, finish_time;
    SpiceChar           utc_center_str[TIMLEN];
    SpiceChar           utc_result_str[TIMLEN];
    SpiceChar           utc_finish_str[TIMLEN];


    printf("Starting Enceladus eclipse by Saturn search...\n\n");

    /*
    Load the SPICE kernels specified in the meta-kernel.
    Ensure your meta.tm file includes an SPK for Enceladus (e.g., sat*.bsp)
    */
    furnsh_c ( "meta.tm" );
    printf("Kernels loaded from %s\n", "meta.tm");


    /*
    Define the time window for the search.
    */
    strcpy(utc_center_str, "2025-05-10 21:45:23 UTC");
    str2et_c(utc_center_str, &et_center);

    // Search window is +/- 20 minutes from the center time
    et_start = et_center - (20.0 * 60.0);
    et_end   = et_center + (20.0 * 60.0);

    wninsd_c(et_start, et_end, &cnfine);
    
    printf("Search window set from TDB %f to %f\n", et_start, et_end);

    /*
    Execute the Geometry Finder search for the occultation.
    This function finds time intervals when the `front_body` blocks the `back_body`
    as seen from the `observer`.
    */
    printf("Executing gfoclt_c search for eclipse of Enceladus by Saturn...\n");
    gfoclt_c( occtyp,
              front_body,  front_shape, front_frame,
              back_body,   back_shape,  back_frame,
              abcorr,      observer,    step,
              &cnfine,     &result );

    /*
    Check the results.
    */
    if ( wncard_c(&result) == 0 )
    {
        printf("\nNo eclipse event found in the time window.\n");
    }
    else
    {
        printf("\nEclipse event(s) found!\n");
        
        // Loop through all found intervals. In this case, it should be just one.
        for (SpiceInt i = 0; i < wncard_c(&result); i++)
        {
            wnfetd_c(&result, i, &start_time, &finish_time);
            
            // Format the start and end times for printing
            timout_c(start_time, TIMFMT, TIMLEN, utc_result_str);
            timout_c(finish_time, TIMFMT, TIMLEN, utc_finish_str);

            printf("\nEclipse Interval #%d:\n", i + 1);
            printf("   Starts (Enters Penumbra): %s\n", utc_result_str);
            printf("   Ends   (Exits Penumbra) : %s\n", utc_finish_str);
        }
    }

    return 0;
}