#include <stdio.h>
#include <string.h>
#include "SpiceUsr.h"

int main()
{
    // --- User-defined Parameters ---
    const char* out_spk       = "enceladus_surface_point.bsp";
    const char* segment_id    = "ENCELADUS SURFACE POINT (0,0)";
    const char* start_time_str= "2025-05-10 21:00:00 UTC";
    const char* stop_time_str = "2025-05-10 22:30:00 UTC";
    
    SpiceInt    target_id     = -602000;
    SpiceInt    center_id     = 602;
    const char* body_fix_frame= "IAU_ENCELADUS"; // Frame of the position

    // Surface point definition
    SpiceDouble lon_deg = 0.0;
    SpiceDouble lat_deg = 0.0;
    SpiceDouble alt     = 0.0;

    // --- Program Logic ---
    SpiceInt    handle;
    SpiceDouble radii[3];
    SpiceDouble re, f;
    SpiceDouble lon_rad, lat_rad;
    SpiceDouble position_vec[3];
    SpiceDouble first_et, last_et;

    // We need TWO state vectors and TWO epochs for degree 1 interpolation
    SpiceDouble states[2][6];
    SpiceDouble epochs[2];
    
    // 1. Load kernels
    furnsh_c("meta.tm");

    // 2. Get Enceladus radii
    bodvrd_c("ENCELADUS", "RADII", 3, (SpiceInt[]){3}, radii);
    re = radii[0];
    f = (radii[0] - radii[2]) / radii[0];
    
    // 3. Calculate the single position vector
    lon_rad = lon_deg * rpd_c();
    lat_rad = lat_deg * rpd_c();
    latrec_c(re, lon_rad, lat_rad, position_vec);

    printf("Calculated position of point in %s frame (km):\n", body_fix_frame);
    printf("  X: %f, Y: %f, Z: %f\n", position_vec[0], position_vec[1], position_vec[2]);

    // 4. Populate the TWO state vectors. They will be identical.
    // First state vector
    states[0][0] = position_vec[0]; // x
    states[0][1] = position_vec[1]; // y
    states[0][2] = position_vec[2]; // z
    states[0][3] = 0.0; // vx
    states[0][4] = 0.0; // vy
    states[0][5] = 0.0; // vz

    // Second state vector (identical to the first)
    states[1][0] = position_vec[0];
    states[1][1] = position_vec[1];
    states[1][2] = position_vec[2];
    states[1][3] = 0.0;
    states[1][4] = 0.0;
    states[1][5] = 0.0;

    // 5. Convert time strings and populate the TWO epochs
    str2et_c(start_time_str, &first_et);
    str2et_c(stop_time_str, &last_et);
    epochs[0] = first_et;
    epochs[1] = last_et;

    // 6. Open a new SPK file
    spkopn_c(out_spk, "SPK for Enceladus Surface Point", 4096, &handle);

    // 7. Write the segment using SPK writer for Type 9 with DEGREE 1
    spkw09_c(handle,
             target_id,            // The ID of our new point
             center_id,            // The body it's relative to (Enceladus)
             body_fix_frame,       // The frame of the input state vectors
             first_et,             // Start time of segment coverage
             last_et,              // Stop time of segment coverage
             segment_id,           // A name for this segment
             1,                    // Degree of polynomial (1 = linear)
             2,                    // Number of states we are writing (TWO)
             (ConstSpiceDouble*)states, // The array of states
             epochs                // The array of epochs
            );

    // Check for errors
    if (failed_c()) {
        printf("\nAn error occurred while writing the SPK file.\n");
        reset_c();
    } else {
        printf("SPK segment written successfully.\n");
    }

    // 8. Close the SPK file
    spkcls_c(handle);
    printf("File '%s' created.\n", out_spk);

    return 0;
}