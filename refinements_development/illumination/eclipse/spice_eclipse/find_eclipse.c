#include <stdio.h>
#include <string.h>
#include "SpiceUsr.h"

#define MAXWIN  2000
#define TIMFMT  "YYYY-MON-DD HR:MN:SC.### UTC"
#define TIMLEN  40

int main()
{
    // --- Setup ---
    // Declare the SPICE windows for confinement and results
    SPICEDOUBLE_CELL(cnfine, MAXWIN);
    SPICEDOUBLE_CELL(result, MAXWIN);

    // Time window for the search, matching your SPK's coverage
    const char* start_time_str = "2025-05-10 21:00:00 UTC";
    const char* stop_time_str  = "2025-05-10 22:30:00 UTC";

    // --- Load Kernels ---
    // Load the meta-kernel with general planetary data
    furnsh_c("meta.tm");
    // Load the specific SPK for our surface point
    furnsh_c("enceladus_surface_point.bsp");
    printf("Kernels loaded.\n");

    // --- Define the Search Window ---
    SpiceDouble et_start, et_stop;
    str2et_c(start_time_str, &et_start);
    str2et_c(stop_time_str, &et_stop);

    // Insert the time interval into our confinement window
    wninsd_c(et_start, et_stop, &cnfine);

    // --- Define the Eclipse Search Parameters ---
    // We are looking for a PARTIAL occultation (penumbra).
    const char* occtyp = "PARTIAL";

    // The body doing the blocking is Saturn.
    const char* front_body  = "SATURN";
    const char* front_shape = "ELLIPSOID";
    const char* front_frame = "IAU_SATURN";

    // The body being blocked is the Sun.
    const char* back_body   = "SUN";
    const char* back_shape  = "ELLIPSOID";
    const char* back_frame  = "IAU_SUN";

    // The observer is our custom surface point.
    const char* observer = "-602000";

    // Use light-time corrections for apparent geometry.
    const char* abcorr = "LT+S"; // Light time correction plus stellar aberration

    // Set the search step size. 10 minutes (600 seconds) is a safe
    // value that is small enough to not miss the event.
    SpiceDouble step = 1.0;
    
    printf("Searching for first entry into Saturn's penumbra...\n");

    // --- Perform the GF Occultation Search ---
    gfoclt_c(occtyp,
             front_body, front_shape, front_frame,
             back_body,  back_shape,  back_frame,
             abcorr,
             observer,
             step,
             &cnfine,
             &result);

    // --- Process and Display the Results ---
    SpiceInt count = wncard_c(&result);

    if (count == 0)
    {
        printf("\nNo penumbral eclipse was found in the specified timeframe.\n");
    }
    else
    {
        // We found at least one interval. The user wants the first one.
        SpiceDouble start_et, finish_et;
        SpiceChar start_utc[TIMLEN];

        // Fetch the start and end of the FIRST interval (index 0).
        wnfetd_c(&result, 0, &start_et, &finish_et);
        
        // Convert the start time to a readable UTC string.
        timout_c(start_et, TIMFMT, TIMLEN, start_utc);

        printf("\nEvent Found!\n");
        printf("The point enters Saturn's penumbra at: %s\n", start_utc);
    }

    return 0;
}