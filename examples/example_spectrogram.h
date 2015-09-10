#ifndef NSGT_EXAMPLE_SPECTROGRAM_H_H
#define NSGT_EXAMPLE_SPECTROGRAM_H_H

#include "../nsgt.h"
#include "example_utils.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#define NSGT_SPECTROGRAM_LINEAR 0
#define NSGT_SPECTROGRAM_DECIBEL 1
#define NSGT_SPECTROGRAM_SONE 2
#define NSGT_SPECTROGRAM_PHON 4

#define NSGT_SPECTROGRAM_GRADIENT_MONO 0
#define NSGT_SPECTROGRAM_GRADIENT_HEAT 1
#define NSGT_SPECTROGRAM_GRADIENT_BLUERED 2
#define NSGT_SPECTROGRAM_GRADIENT_BLACKBLUERED 4
#define NSGT_SPECTROGRAM_GRADIENT_BLACKBLUEREDYELLOW 8
#define NSGT_SPECTROGRAM_GRADIENT_HEATMAP_7 16

void nsgt_util_init_magnitude_map(double *map, size_t dimx, size_t dimy, NsgtResult *result) {
    size_t bins = result->numBins / 2 - 2;
    assert(bins == dimy); //todo: interpolate when necessary
    for (size_t y = 0; y < dimy; y++) {
        for (size_t x = 0; x < dimx; x++) {
            int i = (int) (((int) dimy - 1) - y + 2);

            nsgt_complex_t *bin = result->bins[i];
            int binSize = (int) result->binSize[i];

            double t = (double) x / ((double) dimx) * (binSize - 1);
            int t0 = (int) floor(t);
            int t1 = (int) ceil(t);
            double factor = t - t0;
            nsgt_complex_t r = {
                    (1. - factor) * bin[t0][0] + factor * bin[t1][0],
                    (1. - factor) * bin[t0][1] + factor * bin[t1][1]
            };
            double magnitude = sqrt(r[0] * r[0] + r[1] * r[1]);
            map[y * dimx + x] = magnitude;
        }
    }
}

//todo: cleanup this function (its only a helper but still...)
static void nsgt_gradient(int gradient, double t, double *r, double *g, double *b) {
    assert(0 <= t && t <= 1.);
    if (gradient == NSGT_SPECTROGRAM_GRADIENT_MONO) { //black -> white
        *r = *g = *b = t;
    } else if (gradient == NSGT_SPECTROGRAM_GRADIENT_HEAT) { //black->red->yellow
        if (t < 1. / 2.) {
            *r = t * 2;
        } else {
            *r = 1;
            *g = 2. * t - 1;
        }
    } else if (gradient == NSGT_SPECTROGRAM_GRADIENT_BLUERED) { //blue->red
        *b = 1. - t;
        *r = t;
    } else if (gradient == NSGT_SPECTROGRAM_GRADIENT_BLACKBLUERED) {
        if (t < 1. / 2.) {
            *b = t * 2;
        } else {
            double magn = (t - 1. / 2.) * 2.;
            *b = 1 - magn;
            *r = magn;
        }
    } else if (gradient == NSGT_SPECTROGRAM_GRADIENT_BLACKBLUEREDYELLOW) {
        if (t < 1. / 3.) {
            *b = t * 3.;
        } else if (t < 2. / 3.) {
            double magn = 3 * t - 1;
            *b = 1 - magn;
            *r = magn;
        } else {
            *r = 1.;
            *g = 3 * t - 2;
        }
    } else if (gradient == NSGT_SPECTROGRAM_GRADIENT_HEATMAP_7) {
        if (t < 1. / 6.) { //black->blue
            *b = t * 6;
        } else if (t < 2. / 6.) { //blue->cyan
            *b = 1.;
            *g = t * 6. - 1;
        } else if (t < 3. / 6.) { //cyan->green
            *b = 1. - (t * 6. - 2.);
            *g = 1.;
        } else if (t < 4. / 6.) { //green->yellow
            *g = 1.;
            *r = t * 6. - 3;
        } else if (t < 5. / 6.) { //yellow->red
            *g = 1. - (t * 6. - 4);
            *r = 1.;
        } else { //red->white
            *r = 1.;
            *g = *b = t * 6. - 5.;
        }
    } else {
        assert(0);
    }

    assert(*r <= 1. && *g <= 1. && *b <= 1.);
    assert(*r >= 0. && *g >= 0. && *b >= 0.);
}

typedef struct cubic_interpolator {
    double *xs;
    double *ys;
    double *c1s;
    double *c2s;
    double *c3s;
    int len;
} nsgt_interpolator;

//ported version of https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
//expects xs to be sorted
cubic_interpolator *cubic_interpolator_create(double *xs, double *ys, int len) {
    cubic_interpolator *self = (cubic_interpolator *)malloc(sizeof(cubic_interpolator));
    self->len = len;
    self->xs = (double *) malloc(len * sizeof(double));
    memcpy(self->xs, xs, len * sizeof(double));
    self->ys = (double *) malloc(len * sizeof(double));
    memcpy(self->ys, ys, len * sizeof(double));
    self->c1s = (double *) malloc((len - 1) * sizeof(double));
    self->c2s = (double *) malloc((len - 2) * sizeof(double));
    self->c3s = (double *) malloc((len - 2) * sizeof(double));
    if (len <= 1) return self;

    double *dys = (double *) malloc((len - 1) * sizeof(double));
    double *dxs = (double *) malloc((len - 1) * sizeof(double));
    double *ms = (double *) malloc((len - 1) * sizeof(double));

    for (int i = 0; i < len - 1; i++) {
        double dx = xs[i + 1] - xs[i];
        double dy = ys[i + 1] - ys[i];
        dxs[i] = dx;
        dys[i] = dy;
        ms[i] = dy / dx;
    }
    self->c1s[0] = ms[0];
    for (int i = 0; i < len - 2; i++) {
        double m = ms[i], mNext = ms[i + 1];
        if (m * mNext <= 0) {
            self->c1s[i + 1] = 0;
        } else {
            double dx = dxs[i], dxNext = dxs[i + 1], common = dx + dxNext;
            self->c1s[i + 1] = (3 * common / ((common + dxNext) / m + (common + dx) / mNext));
        }
    }
    self->c1s[len - 2] = ms[len - 2];
    for (int i = 0; i < len - 2; i++) {
        double c1 = self->c1s[i], m = ms[i], invDx = 1 / dxs[i], common = c1 + self->c1s[i + 1] - m - m;
        self->c2s[i] = (m - c1 - common) * invDx;
        self->c3s[i] = common * invDx * invDx;
    }
	free(dys);
	free(dxs);
	free(ms);
    return self;
}

void cubic_interpolator_free(cubic_interpolator *interpolator) {
    free(interpolator->xs);
    free(interpolator->ys);
    free(interpolator->c1s);
    free(interpolator->c2s);
    free(interpolator->c3s);
    free(interpolator);
}

double cubic_interpolator_interpolate(cubic_interpolator *self, double t) {
    if (self->len == 0)
        return .0;
    if (self->len == 1)
        return self->ys[0];
    int i = self->len - 1;
    if (t == self->xs[i])
        return self->ys[i];
    int low = 0, high = self->len - 3;
    while (low <= high) {
        int mid = (low + high) / 2;
        double xc = self->xs[mid];
        if (xc < t) {
            low = mid + 1;
        }
        else if (xc > t) {
            high = mid - 1;
        }
        else {
            return self->ys[mid];
        }
    }
    i = high < 0 ? 0 : high;
    double diff = t - self->xs[i], diffSq = diff * diff;
    return self->ys[i] + self->c1s[i] * diff + self->c2s[i] * diffSq + self->c3s[i] * diff * diffSq;
}

cubic_interpolator *phon_interpolator_create(double freq) {
    static cubic_interpolator *interpolators[7] = {0};
    //values from ISO 226:2003
    //frequencies are linearized so the interpolation algorithm weighs the points correctly
    static double values[8][31] = {
            {.0,     1.,     2.,     3.,     4.,     5.,     6.,     7.,     8.,     9.,     10.,    11.,    12.,   13.,   14.,   15.,   16.,   17.,    18.,    19.,    20.,    21.,   22.,   23.,   24.,   25.,    26.,    27.,    28.,    29.,    30.},
            {76.55,  65.62,  55.12,  45.53,  37.63,  30.86,  25.02,  20.51,  16.65,  13.12,  10.09,  7.54,   5.11,  3.06,  1.48,  0.3,   -0.3,  -0.01,  1.03,   -1.19,  -4.11,  -7.05, -9.03, -8.49, -4.48, 3.28,   9.83,   10.48,  8.38,   14.1,   79.65},
            {83.75,  75.76,  68.21,  61.14,  54.96,  49.01,  43.24,  38.13,  33.48,  28.77,  24.84,  21.33,  18.05, 15.14, 12.98, 11.18, 9.99,  10,     11.26,  10.43,  7.27,   4.45,  3.04,  3.8,   7.46,  14.35,  20.98,  23.43,  22.33,  25.17,  81.47},
            {89.58,  82.65,  75.98,  69.62,  64.02,  58.55,  53.19,  48.38,  43.94,  39.37,  35.51,  31.99,  28.69, 25.67, 23.43, 21.48, 20.1,  20.01,  21.46,  21.4,   18.15,  15.38, 14.26, 15.14, 18.63, 25.02,  31.52,  34.43,  33.04,  34.67,  84.18},
            {99.85,  93.94,  88.17,  82.63,  77.78,  73.08,  68.48,  64.37,  60.59,  56.7,   53.41,  50.4,   47.58, 44.98, 43.05, 41.34, 40.06, 40.01,  41.82,  42.51,  39.23,  36.51, 35.61, 36.65, 40.01, 45.83,  51.8,   54.28,  51.49,  51.96,  92.77},
            {109.51, 104.23, 99.08,  94.18,  89.96,  85.94,  82.05,  78.65,  75.56,  72.47,  69.86,  67.53,  65.39, 63.45, 62.05, 60.81, 59.89, 60.01,  62.15,  63.19,  59.96,  57.26, 56.42, 57.57, 60.89, 66.36,  71.66,  73.16,  68.63,  68.43,  104.92},
            {118.99, 114.23, 109.65, 105.34, 101.72, 98.36,  95.17,  92.48,  90.09,  87.82,  85.92,  84.31,  82.89, 81.68, 80.86, 80.17, 79.67, 80.01,  82.48,  83.74,  80.59,  77.88, 77.07, 78.31, 81.62, 86.81,  91.41,  91.74,  85.41,  84.67,  118.95},
            {128.41, 124.15, 120.11, 116.38, 113.35, 110.65, 108.16, 106.17, 104.48, 103.03, 101.85, 100.97, 100.3, 99.83, 99.62, 99.5,  99.44, 100.01, 102.81, 104.25, 101.18, 98.48, 97.67, 99,    102.3, 107.23, 111.11, 110.23, 102.07, 100.83, 133.73},
    };
    for (int i = 0; i < 7; i++) {
        if (!interpolators[i]) {
            interpolators[i] = cubic_interpolator_create(values[0], values[i + 1], 31);
        }
    }
    static double phon_scale[7] = {.0, 10., 20., 40., 60., 80., 100.};
    double db_scale[7] = {0};
    double freq_index = 10 * log10(freq / 20.);
    for (int i = 0; i < 7; i++) {
        db_scale[i] = cubic_interpolator_interpolate(interpolators[i], freq_index);
    }
    return cubic_interpolator_create(db_scale, phon_scale, 7);
}

double phon_to_sone(double phon) {
    if (phon > 40)
        return pow(2., (phon - 40.) / 10.);
    static cubic_interpolator *intpl = NULL;
	if (!intpl) {
		double p[] = { 1., 9., 11., 14., 19., 25., 32., 40., 50., 60., 70., 80., 90., 100., 110., 120. };
		double s[] = { 0., 1. / 64., 1. / 32, 1. / 16., 1. / 8., 1. / 4, 1. / 2., 1., 2., 4., 8., 16., 32., 64., 128., 256. };
		intpl = cubic_interpolator_create(p, s, 16);
	}
    return cubic_interpolator_interpolate(intpl, phon);
}

void nsgt_scale_map(Nsgt *nsgt, double *map, size_t dimx, size_t dimy, int scalingType) {
    //get min and max
    double min = DBL_MAX, max = -DBL_MAX;
    for (size_t y = 0; y < dimy; y++) {
        for (size_t x = 0; x < dimx; x++) {
            double mag = map[y * dimx + x];
            if (mag < min) min = mag;
            if (mag > max) max = mag;
        }
    }
    //transform to value between 0 and 1
    for (size_t y = 0; y < dimy; y++) {
        nsgt_scalar_t freq = nsgt->fbas[dimy - 1 - y] * nsgt->fs / nsgt->Ls;
        cubic_interpolator *intpl = NULL;
        if (scalingType == NSGT_SPECTROGRAM_SONE || scalingType == NSGT_SPECTROGRAM_PHON) {
            intpl = phon_interpolator_create(freq);
        }
        for (size_t x = 0; x < dimx; x++) {
            double mag = map[y * dimx + x];
            if (scalingType == NSGT_SPECTROGRAM_LINEAR) {
                mag -= min;
                mag /= max - min;
            } else if (scalingType == NSGT_SPECTROGRAM_DECIBEL) {
                mag /= max;
                mag = 10 * log10(mag);
                double range = -10 * log10(min / max);
                mag = (mag + range) / range;
            } else if (scalingType == NSGT_SPECTROGRAM_PHON || scalingType == NSGT_SPECTROGRAM_SONE) {
                mag /= max;
                mag = 10 * log10(mag);
                double range = -10 * log10(min / max);
                mag += range;

                if (mag < 0 || freq < 20 || freq > 20000)
                    mag = 1;
                else
                    mag = cubic_interpolator_interpolate(intpl, mag);

                if (scalingType == NSGT_SPECTROGRAM_SONE) {
                    mag = phon_to_sone(mag);
                }
            } else {
                assert(0);
                return;
            }
            map[y * dimx + x] = mag;
        }
        if (scalingType == NSGT_SPECTROGRAM_SONE || scalingType == NSGT_SPECTROGRAM_PHON) {
            cubic_interpolator_free(intpl);
        }
    }
    //phon and sone needs another pass to find the maximum values and scale all values down to 0->1
    if (scalingType == NSGT_SPECTROGRAM_PHON || NSGT_SPECTROGRAM_SONE) {
        max = -DBL_MAX;
        for (size_t y = 0; y < dimy; y++) {
            for (size_t x = 0; x < dimx; x++) {
                double mag = map[y * dimx + x];
                if (mag > max) max = mag;
            }
        }
        for (size_t y = 0; y < dimy; y++) {
            for (size_t x = 0; x < dimx; x++) {
                map[y * dimx + x] /= max;
            }
        }
    }
}

void nsgt_util_write_spectrogram_ppm(Nsgt *nsgt, char *filename, double *map, size_t dimx, size_t dimy, int scalingType,
                                     int colorType) {
    FILE *fp = NULL;
#ifdef _MSC_VER
    fopen_s(&fp, filename, "wb");
#else
    fp = fopen(filename, "wb");
#endif
    if (!fp) {
        fprintf(stderr, "Unable to open file");
        return;
    }

	double *magnitudes = (double *)malloc(dimx * dimy * sizeof(double));
	memcpy(magnitudes, map, dimx * dimy * sizeof(double));
	nsgt_scale_map(nsgt, magnitudes, dimx, dimy, scalingType);

	fprintf(fp, "P6\n%d %d\n255\n", (int) dimx, (int) dimy);
    for (size_t y = 0; y < dimy; y++) {
        for (size_t x = 0; x < dimx; x++) {
            double mag = magnitudes[y * dimx + x];
            if (mag > 1.)
                mag = 1.;
            if (mag < .0)
                mag = .0;
			double r = .0, g = .0, b = .0;
            nsgt_gradient(colorType, mag, &r, &g, &b);
            static unsigned char color[3] = {0};
            color[0] = (unsigned char) (r * 255);
            color[1] = (unsigned char) (g * 255);
            color[2] = (unsigned char) (b * 255);
            fwrite(color, 1, 3, fp);
        }
    }
    fclose(fp);
	free(magnitudes);
}


void example_spectrogram(double *input, size_t len, int samplerate) {
    printf("Calculating nsgt for spectrograms\r\n");

    NsgtScaleDesc scale = {0};
    scale.scale_type = NSGT_SCALE_OCT;
    scale.fmin = 50.;
    scale.fmax = 22050.;
    scale.num_bins = 48;
    Nsgt *nsgt = nsgt_create(&scale, samplerate, len, 0);

    NsgtResult *result = nsgt_alloc_result(nsgt);

    nsgt_complex_t *input_complex = (nsgt_complex_t *) new nsgt_scalar_t[len * 2];
    nsgt_forward_set_buffers(nsgt, input_complex, result); //lots of time is spent here, initializing the fftw library
    copy_real_to_complex(input, input_complex, len);

    nsgt_forward(nsgt);

    size_t dimx = (size_t) (25 * len / samplerate);
    if (dimx < 2000) dimx = 2000;
    if (dimx > 25000) dimx = 25000;
    size_t dimy = result->numBins / 2 - 2;
    double *magnitude = (double *) malloc(dimx * dimy * sizeof(double));
    nsgt_util_init_magnitude_map(magnitude, dimx, dimy, result);

    printf("Writing spectrograms\r\n");

    int gradient_type = NSGT_SPECTROGRAM_GRADIENT_HEATMAP_7; //NSGT_SPECTROGRAM_GRADIENT_HEAT is also quite nice

    nsgt_util_write_spectrogram_ppm(nsgt, (char *) "spectrogram_sone.ppm", magnitude, dimx, dimy,
                                    NSGT_SPECTROGRAM_SONE, gradient_type);
    nsgt_util_write_spectrogram_ppm(nsgt, (char *) "spectrogram_phon.ppm", magnitude, dimx, dimy,
                                    NSGT_SPECTROGRAM_PHON, gradient_type);
    nsgt_util_write_spectrogram_ppm(nsgt, (char *) "spectrogram_linear.ppm", magnitude, dimx, dimy,
                                    NSGT_SPECTROGRAM_LINEAR, gradient_type);
    nsgt_util_write_spectrogram_ppm(nsgt, (char *) "spectrogram_decibel.ppm", magnitude, dimx, dimy,
                                    NSGT_SPECTROGRAM_DECIBEL, gradient_type);
    free(magnitude);
    nsgt_free_result(result);
    nsgt_free(nsgt);
}

#endif //NSGT_EXAMPLE_SPECTROGRAM_H_H
