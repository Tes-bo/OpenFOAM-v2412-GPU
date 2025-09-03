#ifndef KAHIP_H
#define KAHIP_H

/* *** DUMMY VERSION of kaHIP_interface.h
 *     This file should not be included if you have kahip correctly installed
 *     See: etc/config.sh/kahip
 *     See: src/parallel/decompose/kahipDecomp/Make/options
 */

#warning "Dummy kahip.h - could not find KAHIP installation."

#ifdef __cplusplus
extern "C" {
#endif

void kaffpa
(
    // [inputs]
    int* n,
    int* vwgt,
    int* xadj,
    int* adjcwgt,
    int* adjncy,
    int* nparts,
    double* imbalance,
    bool suppress_output,
    int seed,
    int mode,
    // [outputs]
    int* edgecut,
    int* part
);

#ifdef __cplusplus
}
#endif

#endif
