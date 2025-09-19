#include <stdio.h>
#include "SpiceUsr.h"

int main(void) {
    /* Load the geophysical kernel */
    furnsh_c("/Users/claudioferrara/Downloads/cspice/data/geophysical.ker");

    /* Just print a success message */
    printf("Loaded geophysical.ker OK!\n");

    /* Unload everything and exit */
    kclear_c();
    return 0;
}