#ifndef METIS_H
#define METIS_H

/* *** DUMMY VERSION of metis.h
 *     This file should not be included if you have metis correctly installed.
 *     See: etc/config.sh/metis
 *     See: src/parallel/decompose/metisDecomp/Make/options
 */

#warning "Dummy metis.h - could not find METIS installation."

// Integer type: OpenFOAM uses WM_LABEL_SIZE, metis.h uses IDXTYPEWIDTH
#if WM_LABEL_SIZE == 32
  typedef int32_t idx_t;

  #define IDXTYPEWIDTH 32
  #define SCNIDX  SCNd32
  #define PRIIDX  PRId32
#elif WM_LABEL_SIZE == 64
  typedef int64_t idx_t;

  #define IDXTYPEWIDTH 64
  #define SCNIDX  SCNd64
  #define PRIIDX  PRId64
#else
  #error "Incorrect user-supplied value for WM_LABEL_SIZE  (metis IDXTYPEWIDTH)"
#endif


// Float type: OpenFOAM uses WM_SP, WM_DP, metis.h uses REALTYPEWIDTH
#if defined(WM_SP) || defined(WM_SPDP)
  typedef float real_t;
  #define REALTYPEWIDTH 32
#elif defined(WM_DP)
  typedef double real_t;
  #define REALTYPEWIDTH 64
#else
  #error "Incorrect user-supplied value for WM_SP (WM_SPDP) / WM_DP  (metis REALTYPEWIDTH)"
#endif


#ifdef __cplusplus
extern "C" {
#endif
int            METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
                  idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
                  idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options,
                  idx_t *edgecut, idx_t *part);

int            METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
                  idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
                  idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options,
                  idx_t *edgecut, idx_t *part);

#ifdef __cplusplus
}
#endif


#endif
