/*
 * libnsgt
 *
 * Copyright (c) 2015 libnsgt
 *
 * Implementation of the 'Nonstationary Gabor transform' in C
 *
 * Contact: mr.greywater+libnsgt@gmail.com
 *
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 *
 * Implementation using the papers from:
 * - M. Dörfler, N. Holighaus, T. Grill and G. Velasco, “Constructing an Invertible Constant-Q Transform with Nonstationary Gabor Frames,"
 *   Proceedings of the 14th International Conference on Digital Audio Effects (DAFx 11), Paris, France, 2011
 * - P. Balazs, M. Dörfler, N. Holighaus, F. Jaillet and G. Velasco, “Theory, Implementation and Application of Nonstationary Gabor Frames,"
 *   to appear in Journal of Computational and Applied Mathematics
 *
 * For more information see http://univie.ac.at/nonstatgab/
 */

#ifndef NSGT_NSGT_H
#define NSGT_NSGT_H

#include <stdint.h>
#include <stddef.h>

#define nsgt_data_t double

typedef nsgt_data_t nsgt_scalar_t;

//#define NSGT_USE_COMPLEX
#ifdef NSGT_USE_COMPLEX
#include <complex.h>
typedef nsgt_data_t _Complex nsgt_complex_t;
#else
typedef nsgt_data_t nsgt_complex_t[2];
#endif

#define NSGT_SCALE_OCT 0

typedef struct NsgtScaleDesc {
    int32_t scale_type;
    nsgt_scalar_t fmin;
    nsgt_scalar_t fmax;
    size_t num_bins; //bpo
} NsgtScaleDesc;

typedef struct NsgtScale {
    NsgtScaleDesc desc;
    size_t bnds;
    nsgt_scalar_t pow2n;
    nsgt_scalar_t q;
    nsgt_scalar_t *f;
} NsgtScale;
struct fftw_plan;

typedef struct Nsgt {
    NsgtScale *scale;
    nsgt_scalar_t fs;
    size_t Ls;
    int32_t options;
    size_t fbas_size; //size of rfbas and g array
    double *fbas;
    int *M; // len(M) = fbas_size
    nsgt_scalar_t **g; // len(g) = fbas_size, len(g[i]) = M[i]
    nsgt_scalar_t **gd; // len(gd) = fbas_size, len(gd[i]) = M[i]
    int *ranges; // beginning of each range (mod nn) with length of M[i]
    int nn;
    nsgt_complex_t *temp; //size: max(Ls, nn, biggest window size), used as buffer for various calculations
    int M_total;
    struct NsgtResult *forward_output;
    struct NsgtResult *backward_input;
    nsgt_complex_t *forward_input;
    nsgt_complex_t *backward_output;
    void **forwardplans;
    void **backwardplans;
} Nsgt;

typedef struct NsgtStream {
    Nsgt *nsgt;
    int tr_area;
    int odd_forwards;
    int odd_backwards;

    /* tukey window */
    nsgt_scalar_t *tw;
    /* buffer for last slice */
    nsgt_complex_t *past_forward;
    /* buffer for last slice */
    nsgt_complex_t *past_backwards;
    /* temporary slice buffer */
    nsgt_complex_t *f_slice;

    nsgt_complex_t *forward_input;
    nsgt_complex_t *backward_output;
} NsgtStream;

typedef struct NsgtResult {
    nsgt_complex_t **bins;
    size_t *binSize;
    size_t numBins;
} NsgtResult;

#define NSGT_REDUCEDFORM_0 0 /* default */
#define NSGT_REDUCEDFORM_1 1
#define NSGT_REDUCEDFORM_2 2
#define NSGT_SLICED 4
#define NSGT_OPTIMIZE_FFT 8

#define NSGT_GET_REDUCEDFORM(a) ((a) & (NSGT_REDUCEDFORM_0 | NSGT_REDUCEDFORM_1 | NSGT_REDUCEDFORM_2))

/* nsgt result */
size_t nsgt_result_get_size(Nsgt *nsgt);

void nsgt_init_result(Nsgt *nsgt, NsgtResult *result);

NsgtResult *nsgt_alloc_result(Nsgt *nsgt);

void nsgt_free_result(NsgtResult *result);

/* nsgt */

Nsgt *nsgt_create(NsgtScaleDesc *scale, nsgt_scalar_t samplerate, size_t numsamples, int32_t options);

void nsgt_free(Nsgt *nsgt);

void nsgt_forward_set_buffers(Nsgt *nsgt, nsgt_complex_t *input, NsgtResult *output);

void nsgt_backward_set_buffers(Nsgt *nsgt, NsgtResult *input, nsgt_complex_t *output);

/* These functions are for real time usage, so there must be no memory allocation */
void nsgt_forward(Nsgt *nsgt);

void nsgt_backward(Nsgt *nsgt);

/* nsgt stream */

NsgtStream *nsgt_stream_create(NsgtScaleDesc *scale, nsgt_scalar_t fs /* samplerate */, size_t blocksize,
                               int32_t tr_area, int32_t options);

void nsgt_stream_free(NsgtStream *nsgt);

void nsgt_stream_forward_set_buffers(NsgtStream *nsgt, nsgt_complex_t *input, NsgtResult *output);

void nsgt_stream_backward_set_buffers(NsgtStream *nsgt, NsgtResult *input, nsgt_complex_t *output);

/* These functions are for real time usage, so there must be no memory allocation */
void nsgt_stream_forward(NsgtStream *nsgt);

void nsgt_stream_backward(NsgtStream *nsgt);

void nsgt_stream_reset(NsgtStream *nsgt);

/* Utility */
double nsgt_util_calc_error(nsgt_complex_t *a, nsgt_complex_t *b, size_t length);

#endif //NSGT_NSGT_H
