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
//todo: support float and long double (currently only double)
//todo: support real to complex forward and complex to real backwards transform
//todo: support other fft libraries such as muFFT / Intel Mkl / clFFT / CUFFT (FFTW is GPLv2 licensed!)
//      maybe have some way to support any fft library by registering fft functions
//todo: multithreading options
//todo: have a simple small fft library included (MIT licensed), to be able to use this library without dependencies
//todo: support other scaling apart from oct
//todo: export and import fft wisdom (possibly also serialize and store nsgt struct)
//todo: error handling (what if malloc returns 0)

#include "nsgt.h"
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <fftw3.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef M_PI
static const nsgt_scalar_t M_PI = 3.14159265358979323846;
#endif
#if !defined(_DEBUG) && !defined(NDEBUG)
#define NDEBUG
#endif
#define INT_FLOOR_DIVISION(x, y) ( x / y - (x % y < 0)) //y may not be negative
#define INT_MOD(x, y) (((x % y) + y) % y) //modulo towards -inf (C remainder is towards 0)
#define FFT_SHIFT(i, n) ((i + ((n + 1) / 2)) % n)
#define FFT_ISHIFT(i, n) ((i + (n / 2)) % n)

/* restricted declarations */
void nsgt_init(Nsgt *nsgt, nsgt_scalar_t *f, nsgt_scalar_t q, size_t bnds, nsgt_scalar_t sr, size_t Ls, int sliced,
               int min_win, int Qvar);

NsgtScale *nsgt_scale_create(NsgtScaleDesc *desc);

void nsgt_scale_free(NsgtScale *scale);

void hannwin_shift_ifft(nsgt_scalar_t *buf, size_t len);

void arrange(int orientation, int mkk, nsgt_complex_t *t, nsgt_complex_t *ckk);


/* nsgt result */
size_t nsgt_result_get_size(Nsgt *nsgt) {
    size_t len = sizeof(NsgtResult);
    len += nsgt->fbas_size * sizeof(void *);
    len += nsgt->M_total * sizeof(nsgt_complex_t);
    len += nsgt->fbas_size * sizeof(size_t);
    return len;
}

void nsgt_init_result(Nsgt *nsgt, NsgtResult *result) {
    char *c = (char *) result;
    size_t binOffset = sizeof(NsgtResult);
    result->bins = (nsgt_complex_t **) (c + binOffset);
    size_t windowOffset = binOffset + nsgt->fbas_size * sizeof(void *);
    for (size_t i = 0; i < nsgt->fbas_size; i++) {
        result->bins[i] = (nsgt_complex_t *) (c + windowOffset);
        windowOffset += nsgt->M[i] * sizeof(nsgt_complex_t);
    }
    result->binSize = (size_t *) (c + windowOffset);
}

NsgtResult *nsgt_alloc_result(Nsgt *nsgt) {
    NsgtResult *result = fftw_malloc(nsgt_result_get_size(nsgt));
    nsgt_init_result(nsgt, result);
    return result;
}

void nsgt_free_result(NsgtResult *result) {
    fftw_free(result);
}

/* nsgt */

Nsgt *nsgt_create(NsgtScaleDesc *scale, nsgt_scalar_t samplerate, size_t numsamples, int32_t options) {
    Nsgt *nsgt = malloc(sizeof(Nsgt));

    assert(samplerate > 0);
    assert(numsamples > 0);
    int reducedform = NSGT_GET_REDUCEDFORM(options);
    assert(0 <= reducedform && reducedform <= 2);
    int min_win = (options & NSGT_SLICED) ? 16 : 4; //todo: add as option

    nsgt->fs = samplerate;
    nsgt->Ls = numsamples;
    nsgt->options = options;
    nsgt->scale = nsgt_scale_create(scale);

    nsgt_init(nsgt, nsgt->scale->f, nsgt->scale->q, nsgt->scale->bnds, samplerate, numsamples,
              (options & NSGT_SLICED) != 0, min_win, 1);

    nsgt_forward_set_buffers(nsgt, NULL, NULL);
    nsgt_backward_set_buffers(nsgt, NULL, NULL);
    return nsgt;
}

void nsgt_free(Nsgt *nsgt) {
    nsgt_forward_set_buffers(nsgt, NULL, NULL);
    nsgt_backward_set_buffers(nsgt, NULL, NULL);
    nsgt_scale_free(nsgt->scale);
    free(nsgt->M);
    fftw_free(nsgt->g);
    fftw_free(nsgt->gd);
    free(nsgt->ranges);
    fftw_free(nsgt->temp);
    memset(nsgt, 0, sizeof(Nsgt));
    free(nsgt);
}

void nsgt_forward_set_buffers(Nsgt *nsgt, nsgt_complex_t *input, NsgtResult *output) {
    unsigned int optimize = (nsgt->options & NSGT_OPTIMIZE_FFT) ? FFTW_MEASURE : FFTW_ESTIMATE;
    if (nsgt->forward_input != input) {
        nsgt->forward_input = input;
        if (nsgt->forwardplans[0] != NULL) {
            fftw_destroy_plan(nsgt->forwardplans[0]);
        }
        if (input != NULL) {
            nsgt->forwardplans[0] = fftw_plan_dft_1d((int) nsgt->Ls, nsgt->forward_input, nsgt->temp, FFTW_FORWARD,
                                                     optimize);
        } else {
            nsgt->forwardplans[0] = NULL;
        }
    }
    if (nsgt->forward_output != output) {
        nsgt->forward_output = output;
        for (size_t i = 0; i < nsgt->fbas_size; i++) {
            if (nsgt->forwardplans[i + 1] != NULL) {
                fftw_destroy_plan(nsgt->forwardplans[i + 1]);
            }
            if (output != NULL) {
                nsgt->forwardplans[i + 1] = fftw_plan_dft_1d(nsgt->M[i], output->bins[i], output->bins[i],
                                                             FFTW_BACKWARD, optimize);
            } else {
                nsgt->forwardplans[i + 1] = NULL;
            }
        }
    }
}

void nsgt_backward_set_buffers(Nsgt *nsgt, NsgtResult *input, nsgt_complex_t *output) {
    unsigned int optimize = (nsgt->options & NSGT_OPTIMIZE_FFT) ? FFTW_MEASURE : FFTW_ESTIMATE;
    if (nsgt->backward_output != output) {
        nsgt->backward_output = output;
        if (nsgt->backwardplans[0] != NULL) {
            fftw_destroy_plan(nsgt->backwardplans[0]);
        }
        if (output != NULL) {
            nsgt->backwardplans[0] = fftw_plan_dft_1d(nsgt->nn, output, output, FFTW_BACKWARD, optimize);
        } else {
            nsgt->backwardplans[0] = NULL;
        }
    }
    if (nsgt->backward_input != input) {
        nsgt->backward_input = input;
        for (size_t i = 0; i < nsgt->fbas_size; i++) {
            if (nsgt->backwardplans[i + 1] != NULL) {
                fftw_destroy_plan(nsgt->backwardplans[i + 1]);
            }
            if (output != NULL) {
                nsgt->backwardplans[i + 1] = fftw_plan_dft_1d(nsgt->M[i], nsgt->backward_input->bins[i], nsgt->temp,
                                                              FFTW_FORWARD, optimize);
            } else {
                nsgt->backwardplans[i + 1] = NULL;
            }
        }
    }
}

void nsgt_forward(Nsgt *nsgt) {
    nsgt_scalar_t **g = nsgt->g;
    int *wins = nsgt->ranges;
    int nn = nsgt->nn;
    int *M = nsgt->M;
    int Ls = (int) nsgt->Ls; //Count of Frames/Samples
    fftw_execute(nsgt->forwardplans[0]);

    nsgt_complex_t *ft = nsgt->temp;
    if (nn > Ls) {
        assert(0);  //todo: check
        memset(&ft[Ls - 1], 0, (nn - Ls) * sizeof(*ft));
    }

    size_t *binSize = nsgt->forward_output->binSize;
    for (size_t sl = 0; sl < nsgt->fbas_size; sl++) {
        nsgt_complex_t *tc = nsgt->forward_output->bins[sl];
        int mii = M[sl];
        int win = wins[sl];
        nsgt_scalar_t *gi = g[sl];
        for (int i = 0; i < mii; i++) {
            size_t win_i = (win + FFT_ISHIFT(i, mii)) % nn; //todo: consider precomputed index buffer (benchmark)
#ifdef NSGT_USE_COMPLEX
            tc[i] = gi[i] * ft[win_i];
#else
            tc[i][0] = gi[i] * ft[win_i][0];
            tc[i][1] = gi[i] * ft[win_i][1];
#endif
        }
        fftw_execute(nsgt->forwardplans[sl + 1]);
        binSize[sl] = (size_t) mii;
    }

    nsgt->forward_output->numBins = nsgt->fbas_size;
}

void nsgt_backward(Nsgt *nsgt) {
    nsgt_scalar_t **gd = nsgt->gd;
    int nn = nsgt->nn;
    int *M = nsgt->M;

    assert(nsgt->backward_input->numBins == nsgt->fbas_size);

    size_t ln = nsgt->backward_input->numBins;

    //assuming |g[i]| == M[i]

    nsgt_complex_t *fr = nsgt->backward_output;
    memset(fr, 0, sizeof(*fr) * nn);

    for (size_t i = 0; i < ln; i++) {
        nsgt_scalar_t *gdi = gd[i];
        int mii = M[i];
        int win = nsgt->ranges[i];
        //fft
        nsgt_complex_t *fc = nsgt->temp;
        fftw_execute(nsgt->backwardplans[i + 1]);
        //todo: SSE SIMD? This loop might benefit from vectorization
        //todo: index buffer? (might not be faster, considering we're probably already at max memory throughput)
        for (int j = 0; j < mii; j++) {
            size_t win_i = (win + FFT_ISHIFT(j, mii)) % nn;
#ifdef NSGT_USE_COMPLEX
            fr[win_i] += fc[j] * gdi[j];
#else
            fr[win_i][0] += fc[j][0] * gdi[j];
            fr[win_i][1] += fc[j][1] * gdi[j];
#endif
        }
    }

    fftw_execute(nsgt->backwardplans[0]);
}

/* nsgt stream */

NsgtStream *nsgt_stream_create(NsgtScaleDesc *scale, nsgt_scalar_t fs /* samplerate */, size_t blocksize,
                               int32_t tr_area, int32_t options) {
    NsgtStream *nsgt_stream = malloc(sizeof(NsgtStream));
    assert(tr_area >= 0);
    assert(blocksize > tr_area);
    assert(blocksize % 2 == 0); //todo: this shouldn't be a requirement, but it is. fix it!

    int sl_len = (int) (2 * blocksize);
    nsgt_stream->nsgt = nsgt_create(scale, fs, (size_t) sl_len, options | NSGT_SLICED);
    nsgt_stream->tr_area = tr_area;

    int hhop = sl_len / 4;
    int htr = tr_area / 2;

    nsgt_scalar_t *tw = fftw_malloc(sl_len * sizeof(nsgt_scalar_t));
    {
        nsgt_scalar_t *w = fftw_malloc(2 * tr_area * sizeof(nsgt_scalar_t));
        hannwin_shift_ifft(w, (size_t) (2 * tr_area));

        //todo: replace with memset/memcpy
        //tukey window
        for (int i = 0; i < hhop - htr; i++)
            tw[i] = .0;
        for (int i = 0; i < 2 * htr; i++)
            tw[hhop - htr + i] = w[tr_area + i];
        for (int i = hhop + htr; i < 3 * hhop - htr; i++)
            tw[i] = 1.;
        for (int i = 0; i < (3 * hhop + htr) - (3 * hhop - htr); i++)
            tw[3 * hhop - htr + i] = w[i];
        for (int i = 3 * hhop + htr; i < sl_len; i++)
            tw[i] = .0;

        fftw_free(w);
    }

    nsgt_stream->tw = tw;

    nsgt_stream->past_forward = fftw_malloc(2 * hhop * sizeof(*nsgt_stream->past_forward));
    memset(nsgt_stream->past_forward, 0, 2 * hhop * sizeof(*nsgt_stream->past_forward));

    nsgt_stream->past_backwards = fftw_malloc(2 * hhop * sizeof(*nsgt_stream->past_backwards));
    memset(nsgt_stream->past_backwards, 0, 2 * hhop * sizeof(*nsgt_stream->past_backwards));

    nsgt_stream->odd_forwards = 0;
    nsgt_stream->odd_backwards = 0;

    nsgt_complex_t *f_slice = fftw_malloc(sl_len * sizeof(*f_slice));
    nsgt_stream->f_slice = f_slice;

    return nsgt_stream;
}

void nsgt_stream_free(NsgtStream *nsgt) {
    nsgt_free(nsgt->nsgt);
    fftw_free(nsgt->tw);
    fftw_free(nsgt->past_forward);
    fftw_free(nsgt->past_backwards);
    fftw_free(nsgt->f_slice);
    memset(nsgt, 0, sizeof(NsgtStream));
    free(nsgt);
}

void nsgt_stream_forward_set_buffers(NsgtStream *nsgt, nsgt_complex_t *input, NsgtResult *output) {
    nsgt->forward_input = input;
    nsgt_forward_set_buffers(nsgt->nsgt, nsgt->f_slice, output);
}

void nsgt_stream_backward_set_buffers(NsgtStream *nsgt, NsgtResult *input, nsgt_complex_t *output) {
    nsgt->backward_output = output;
    nsgt_backward_set_buffers(nsgt->nsgt, input, nsgt->f_slice);
}

void nsgt_stream_forward(NsgtStream *nsgt_stream) {
    NsgtResult *res = nsgt_stream->nsgt->forward_output;
    nsgt_complex_t *f = nsgt_stream->forward_input;

    int sl_len = (int) nsgt_stream->nsgt->Ls;
    int hhop = sl_len / 4;
    nsgt_complex_t *past = nsgt_stream->past_forward;
    nsgt_complex_t *f_slice = nsgt_stream->f_slice;
    memset(f_slice, 0, sl_len * sizeof(*f_slice));

    for (int j = 0; j < 4; j++) {
        nsgt_complex_t *pi = j < 2 ? past + (j * hhop) : f + ((j - 2) * hhop);
        nsgt_scalar_t *twi = nsgt_stream->tw + (j * hhop);
        int sli_offset = ((j + 3 - nsgt_stream->odd_forwards * 2) % 4) * hhop; //todo: replace with arrange
        for (int sli = 0; sli < hhop; sli++) {
#ifdef NSGT_USE_COMPLEX
            f_slice[sli + sli_offset] = pi[sli] * twi[sli];
#else
            f_slice[sli + sli_offset][0] = pi[sli][0] * twi[sli];
            f_slice[sli + sli_offset][1] = pi[sli][1] * twi[sli];
#endif
        }
    }
    memcpy(past, f, hhop * 2 * sizeof(*past));
    nsgt_forward(nsgt_stream->nsgt); //nsgt->f_slice --> output

    for (size_t j = 0; j < nsgt_stream->nsgt->fbas_size; j++) {
        int mkk = nsgt_stream->nsgt->M[j];
        nsgt_complex_t *ckk = res->bins[j];
        arrange(nsgt_stream->odd_forwards, mkk, nsgt_stream->nsgt->temp, ckk);
    }

    nsgt_stream->odd_forwards ^= 1;
}

void nsgt_stream_backward(NsgtStream *nsgt_stream) {
    NsgtResult *res = nsgt_stream->nsgt->backward_input;

    for (size_t j = 0; j < nsgt_stream->nsgt->fbas_size; j++) {
        int mkk = nsgt_stream->nsgt->M[j];
        nsgt_complex_t *ckk = res->bins[j];
        arrange(!nsgt_stream->odd_backwards, mkk, nsgt_stream->nsgt->temp, ckk);
    }
    int sl_len = (int) nsgt_stream->nsgt->Ls;
    int hhop = sl_len / 4;
    //sl_len out

    nsgt_complex_t *quad = nsgt_stream->f_slice;

    //todo: we cannot use temp as output here, since the backward transform needs it themselves
    nsgt_backward(nsgt_stream->nsgt); //input->nsgt --> nsgt->f_slice

    nsgt_complex_t *output = nsgt_stream->backward_output;
    memset(output, 0, 2 * hhop * sizeof(*output));


    memcpy(output, nsgt_stream->past_backwards, 2 * hhop * sizeof(*output));
    memset(nsgt_stream->past_backwards, 0, 2 * hhop *
                                           sizeof(*nsgt_stream->past_backwards)); // todo: instead of clearing, just assign instead of add in the following loop

    /* basically the same as arrange, but with an add instead of set (also switching buffers in between) */
    for (int j = 0; j < 4; j++) {
        int sli_offset = ((j + 3 - nsgt_stream->odd_backwards * 2) % 4);
        nsgt_complex_t *isl = quad + sli_offset * hhop;
        nsgt_complex_t *osl = j < 2 ? output + (j * hhop) : nsgt_stream->past_backwards + ((j - 2) * hhop);
        for (int sli = 0; sli < hhop; sli++) {
#ifdef NSGT_USE_COMPLEX
            osl[sli] += isl[sli];
#else
            osl[sli][0] += isl[sli][0];
            osl[sli][1] += isl[sli][1];
#endif
        }
    }

    nsgt_stream->odd_backwards ^= 1;
}

void nsgt_stream_reset(NsgtStream *nsgt) {
    nsgt->odd_forwards = 0;
    nsgt->odd_backwards = 0;
    size_t hhop = nsgt->nsgt->Ls / 4;
    memset(nsgt->past_backwards, 0, 2 * hhop * sizeof(*nsgt->past_backwards));
    memset(nsgt->past_forward, 0, 2 * hhop * sizeof(*nsgt->past_forward));
}

double nsgt_util_calc_error(nsgt_complex_t *a, nsgt_complex_t *b, size_t length) {
    double suma = .0, sumb = .0;
    for (size_t h = 0; h < length; h++) {
#ifdef NSGT_USE_COMPLEX
        nsgt_complex_t diff = a[h] - b[h];
        suma += cabs(diff * diff);
        sumb += cabs(b[h] * b[h]);
#else
        nsgt_complex_t diff = {a[h][0] - b[h][0], a[h][1] - b[h][1]};
        suma += diff[0] * diff[0] + diff[1] * diff[1];
        sumb += b[h][0] * b[h][0] + b[h][1] * b[h][1];
#endif
    }
    return sqrt(suma / sumb);
}

NsgtScale *nsgt_scale_create(NsgtScaleDesc *desc) {
    NsgtScale *scale = malloc(sizeof(NsgtScale));
    scale->desc = *desc;
    if (desc->scale_type == NSGT_SCALE_OCT) {
        scale->bnds = ((size_t) ceil(log2(desc->fmax / desc->fmin) * desc->num_bins)) + 1;
        scale->pow2n = pow(2, 1. / (nsgt_scalar_t) desc->num_bins);
        scale->q = sqrt(scale->pow2n) / (scale->pow2n - 1.) / 2.;
        scale->f = malloc(scale->bnds * sizeof(*scale->f));
        for (size_t i = 0; i < scale->bnds; i++) {
            scale->f[i] = desc->fmin * pow(scale->pow2n, (double)i);
        }
    } else {
        memset(scale, 0, sizeof(NsgtScale));
        assert(0);
    }
    return scale;
}

void nsgt_scale_free(NsgtScale *scale) {
    free(scale->f);
    memset(scale, 0, sizeof(NsgtScale));
    free(scale);
}

/* ckk has size mkk, t has size mkk / 4 */
void arrange(int orientation, int mkk, nsgt_complex_t *t, nsgt_complex_t *ckk) {
    int mkk_quarter = mkk / 4;
    int mkk_three_quarter = mkk - mkk_quarter;
    if (orientation) { //ABCD->BCDA
        memcpy(t, ckk, mkk_quarter * sizeof(nsgt_complex_t));
        memmove(ckk, ckk + mkk_quarter, mkk_three_quarter * sizeof(nsgt_complex_t));
        memcpy(ckk + mkk_three_quarter, t, mkk_quarter * sizeof(nsgt_complex_t));
    }
    else { //ABCD->DABC
        memcpy(t, ckk + mkk_three_quarter, mkk_quarter * sizeof(nsgt_complex_t));
        memmove(ckk + mkk_quarter, ckk, mkk_three_quarter * sizeof(nsgt_complex_t));
        memcpy(ckk, t, mkk_quarter * sizeof(nsgt_complex_t));
    }
}

/*
Hann Window (inverse FFT shifted)
*/
void hannwin_shift_ifft(nsgt_scalar_t *buf, size_t len) {
    for (size_t n = 0; n < len; n++) {
        buf[n] = 0.5 * (cos(n * M_PI * 2. / len) + 1);
    }
}

/*
Blackman Harris Window (inverse FFT shifted)
*/
void blackharr_shift_ifft(nsgt_scalar_t *buf, size_t len) {
    double nn = (double)(len / 2 * 2);
    for (size_t n = 0; n < len; n++) {
        size_t k = FFT_ISHIFT(n, len);
        buf[n] = 0.35872 - 0.48832 * cos(k * (2. * M_PI / nn)) + 0.14128 * cos(k * (4. * M_PI / nn)) -
                 0.01168 * cos(k * (6. * M_PI / nn));
    }
}

void nsgt_fixup_frequencies(nsgt_scalar_t **freqs, size_t *bnds, nsgt_scalar_t q,
                            nsgt_scalar_t sr) { //todo: do the fixup inside the actual scale
    size_t fq_len = *bnds;
    nsgt_scalar_t *f = *freqs;

    for (size_t i = 0; i < fq_len; i++) {
        if (f[i] > 0) {
            f = &f[i];
            //q = &q[i];
            fq_len -= i;
            break;
        }
    }

    for (size_t i = 0; i < fq_len; i++) {
        if (2 * f[i] >= sr) { //nyquist
            fq_len = i;
            break;
        }
    }

#ifndef NDEBUG
    nsgt_scalar_t freq = 0;
    for (size_t i = 0; i < fq_len; i++) {
        assert(f[i] > freq);
        freq = f[i];
        assert(q /*[i]*/ > 0);
    }
#endif

    *bnds = fq_len;
    *freqs = f;
}

void nsgt_init_fbas(nsgt_scalar_t *fbas, size_t fbas_size, nsgt_scalar_t *f, size_t Ls, nsgt_scalar_t sr) {
    size_t lbas = fbas_size / 2 - 1;
    // 0 .. f .. nf .. sr - f[(n-1) - i]
    fbas[0] = .0;
    memcpy(&fbas[1], f, lbas * sizeof(*f));
    fbas[lbas + 1] = sr / 2.;
    for (size_t i = 0; i < lbas; i++) {
        fbas[fbas_size - i - 1] = sr - f[i];
    }
    for (size_t i = 0; i < fbas_size; i++) {
        fbas[i] *= (nsgt_scalar_t) Ls / sr;
    }
#ifndef NDEBUG
    nsgt_scalar_t last = .0;
    for (size_t i = 0; i < fbas_size; i++) {
        assert(last <= fbas[i]);
        last = fbas[i];
    }
#endif
}

void nsgt_init_window_sizes(int *M, size_t fbas_size, nsgt_scalar_t *fbas, nsgt_scalar_t q, nsgt_scalar_t q_var,
                            int min_win, size_t Ls, int sliced) {
    size_t bnds = fbas_size / 2 - 1;
    if (sliced) {
#define PREP(a) (((int)(round((a)*(q_var/4.))))*4)
        M[0] = PREP(2 * fbas[1]);
        M[1] = PREP(fbas[1] / q /*[0]*/);
        for (size_t i = 2; i < bnds + 2; i++) {
            M[i] = PREP(fbas[i + 1] - fbas[i - 1]);
        }
        M[bnds] = PREP(fbas[bnds] / q /*[lbas-1]*/);
        for (size_t i = 0; i < fbas_size - (bnds + 2); i++) {
            M[i + bnds + 2] = PREP(M[bnds - i]);
        }
#undef PREP
    }
    else {
        memset(M, 0, fbas_size * sizeof(*M));
        M[0] = (int) round(2 * fbas[1]);
        for (size_t i = 1; i < fbas_size - 1; i++) {
            M[i] = (int) round(fbas[i + 1] - fbas[i - 1]);
        }
        M[fbas_size - 1] = (int) round(Ls - fbas[fbas_size - 2]);
    }
#ifndef NDEBUG
    for (size_t i = 0; i < fbas_size; i++) {
        assert(M[i] > 0);
    }
#endif
    for (size_t i = 0; i < fbas_size; i++) {
        if (M[i] < min_win)
            M[i] = min_win;
    }
}

int nsgt_init_ranges(int *ranges, size_t fbas_size, nsgt_scalar_t *fbas, int *M, size_t Ls, int sliced) {
    int *rfbas = malloc(fbas_size * sizeof(int));
    size_t bnds = fbas_size / 2 - 1;
    if (sliced) {
        for (size_t i = 0; i < fbas_size; i++) {
            rfbas[i] = ((int) round(fbas[i] / 2.)) * 2;
        }
    } else {
        fbas[bnds] = (fbas[bnds - 1] + fbas[bnds + 1]) / 2;
        fbas[bnds + 2] = Ls - fbas[bnds];
        for (size_t i = 0; i < fbas_size; i++) {
            rfbas[i] = (int) round(fbas[i]);
        }
    }

    int timepos = 0;
    for (size_t i = 0; i < fbas_size; i++) {
        int shift = i == 0 ? 0 : rfbas[i] - rfbas[i - 1]; //may be better if integer
        timepos += shift;
        ranges[i] = -(M[i] / 2) + timepos;
    }
    int nn = timepos + INT_MOD(-rfbas[fbas_size - 1], (int) Ls);

    free(rfbas);

    //fix window start so it's positive. (that way, the c-remainder (% nn) is equal to the modulo (mod nn) operation)
    for (size_t i = 0; i < fbas_size; i++) {
        int win = ranges[i] % nn;
        ranges[i] = win < 0 ? win + nn : win;
    }

    assert(nn ==
           Ls); /* todo: find a use case where this is fails and debug the code, WARNING: some buffers may overflow */

    return nn;
}

void nsgt_init_g(nsgt_scalar_t **g, size_t num_bins, int *M, int sliced) {
    nsgt_scalar_t *gi = (nsgt_scalar_t *) (g + num_bins);
    for (size_t i = 0; i < num_bins; i++) {
        g[i] = gi;
        if (sliced) {
            blackharr_shift_ifft(gi, (size_t) M[i]);
        } else {
            hannwin_shift_ifft(gi, (size_t) M[i]);
        }
        gi += M[i];
    }
    if (sliced) {
        for (size_t i = 1; i < num_bins / 2 + 1; i++) {
            if (M[i - 1] > M[i]) {
                for (int j = 0; j < M[i - 1]; j++)
                    g[i - 1][j] = 1.;
                int offset = (M[i - 1] / 2) - (M[i] / 2);
                hannwin_shift_ifft(g[i - 1] + offset, (size_t) M[i]);
            }
        }
    }
}

void nsgt_init_g_inverse(nsgt_scalar_t **gd, size_t fbas_size, nsgt_scalar_t **g, int *M, int *ranges, int nn) {
    nsgt_scalar_t *gdi = (nsgt_scalar_t *) (gd + fbas_size);

    nsgt_scalar_t *x = malloc(nn * sizeof(nsgt_scalar_t));
    memset(x, 0, nn * sizeof(*x));

    for (size_t i = 0; i < fbas_size; i++) {
        nsgt_scalar_t *gi = g[i];
        int mii = M[i];
        nsgt_scalar_t miif = (nsgt_scalar_t) mii;
        int sl = ranges[i];
        for (int j = 0; j < mii; j++) {
            nsgt_scalar_t xa_sqrt = gi[FFT_SHIFT(j, mii)];

            int sl_i = (sl + j) % nn;
            x[sl_i] += xa_sqrt * xa_sqrt * miif;
        }
    }

    for (size_t i = 0; i < fbas_size; i++) {
        nsgt_scalar_t *gi = g[i];
        int mii = M[i];
        int wi = ranges[i];
        gd[i] = gdi;
        for (int j = 0; j < mii; j++) {
            int wi_i = (wi + FFT_ISHIFT(j, mii)) % nn;
            gdi[j] = gi[j] / x[wi_i] / nn;
        }
        gdi += M[i];
    }

    free(x);
}

void nsgt_init(Nsgt *nsgt, nsgt_scalar_t *f, nsgt_scalar_t q /* should be array later */, size_t bnds,
               nsgt_scalar_t sr, size_t Ls, int sliced, int min_win, int Qvar) {

    nsgt_fixup_frequencies(&f, &bnds, q, sr);

#ifndef NDEBUG
    nsgt_scalar_t q_factor = (Ls / (8. * sr));
    for (size_t i = 0; i < bnds; i++) {
        if (q /*[i]*/ > f[bnds - 1 - i] * q_factor) {
            printf("Q-factor too high up to frequency %.2f\r\n", f[bnds - 1 - i]);
            break;
        }
    }
#endif

    size_t fbas_size = 2 * bnds + 2;
    int nn = 0, total_window_len = 0, maximum_window_len = 0;
    nsgt_scalar_t *fbas;
    int *ranges = NULL, *M = NULL;
    {
        fbas = fftw_malloc(fbas_size *
                           sizeof(nsgt_scalar_t)); //todo: better save it as member as it is useful for querying the frequencies of each bin
        nsgt_init_fbas(fbas, fbas_size, f, Ls, sr);

        M = malloc(fbas_size * sizeof(int)); // ALLOCATION
        nsgt_init_window_sizes(M, fbas_size, fbas, q, Qvar, min_win, Ls, sliced);

        for (size_t i = 0; i < fbas_size; i++) {
            total_window_len += M[i];
            if (M[i] > maximum_window_len)
                maximum_window_len = M[i];
        }

        ranges = malloc(fbas_size * sizeof(int)); // ALLOCATION
        nn = nsgt_init_ranges(ranges, fbas_size, fbas, M, Ls, sliced);
    }

    size_t full_window_size = fbas_size * sizeof(void *) + total_window_len * sizeof(nsgt_scalar_t);

    nsgt_scalar_t **g = fftw_malloc(full_window_size); // ALLOCATION
    nsgt_init_g(g, fbas_size, M, sliced);

    nsgt_scalar_t **gd = fftw_malloc(full_window_size); // ALLOCATION
    nsgt_init_g_inverse(gd, fbas_size, g, M, ranges, nn);

    //normalize (if the fft library doesn't scale its result
    for (size_t i = 0; i < fbas_size; i++) {
        for (int j = 0; j < M[i]; j++) {
            g[i][j] /= M[i];
            gd[i][j] *= M[i];
        }
    }

    //allocate temporary buffer
    size_t temp_size = Ls < (size_t) maximum_window_len ? (size_t) maximum_window_len : Ls;
    temp_size = temp_size < (size_t) nn ? (size_t) nn : temp_size;
    nsgt_complex_t *temp = (nsgt_complex_t *) fftw_malloc(sizeof(nsgt_complex_t) * temp_size); // ALLOCATION

    void **forwardplans = malloc((fbas_size + 1) * sizeof(void *));
    memset(forwardplans, 0, (fbas_size + 1) * sizeof(void *));
    void **backwardplans = malloc((fbas_size + 1) * sizeof(void *));
    memset(backwardplans, 0, (fbas_size + 1) * sizeof(void *));

    nsgt->forwardplans = forwardplans;
    nsgt->backwardplans = backwardplans;
    nsgt->fbas_size = fbas_size;
    nsgt->fbas = fbas;
    nsgt->M = M;
    nsgt->M_total = total_window_len;
    nsgt->ranges = ranges;
    nsgt->nn = nn;
    nsgt->g = g;
    nsgt->gd = gd;
    nsgt->temp = temp;
}
