#include "SpiceUsr.h"
#include <stdio.h>

int main() {
    /* load the leapseconds kernel */
    furnsh_c(
      "/Users/claudioferrara/Documents/GitHub/GIGANTES-Repo/illumination/eph/kernels/naif0012.tls"
    );

    SpiceDouble et;
    str2et_c("2000 JAN 01 12:00:00", &et);
    printf("ET: %f\n", et);

    /* unload all kernels before exit (optional) */
    kclear_c();

    return 0;
}