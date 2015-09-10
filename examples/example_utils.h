#ifndef NSGT_EXAMPLE_UTILS_H
#define NSGT_EXAMPLE_UTILS_H

#include "../nsgt.h"
#ifdef __unix__

#include <time.h>

typedef struct timespec time_point_t;

time_point_t timer_start() {
    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    return start_time;
}

double timer_end(time_point_t start_time) {
	time_point_t end_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    double diffms = (end_time.tv_nsec - start_time.tv_nsec) / 1000000.;
    diffms += (end_time.tv_sec - start_time.tv_sec) * 1000;
    return diffms;
}

#elif _MSC_VER

#include <Windows.h>

typedef LARGE_INTEGER time_point_t;

time_point_t timer_start() {
	LARGE_INTEGER start_time;
	QueryPerformanceCounter(&start_time);
	return start_time;
}

double timer_end(time_point_t start_time) {
	time_point_t end_time;
	QueryPerformanceCounter(&end_time);
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	return static_cast<double>(end_time.QuadPart - start_time.QuadPart) / frequency.QuadPart * 1000;
}

#else
#error "unsupported platform"
#endif

void copy_real_to_complex(nsgt_scalar_t *real, nsgt_complex_t *cmplx, size_t len) {
    typedef nsgt_scalar_t c[2];
    c *cmplx_arr = reinterpret_cast<c *>(cmplx);
    for (size_t i = 0; i < len; i++) {
        cmplx_arr[i][0] = real[i];
        cmplx_arr[i][1] = .0;
    }
}

void copy_complex_to_real(nsgt_complex_t *cmplx, nsgt_scalar_t *real, size_t len) {
    typedef nsgt_scalar_t c[2];
    c *cmplx_arr = reinterpret_cast<c *>(cmplx);
    for (size_t i = 0; i < len; i++) {
        real[0] = cmplx_arr[i][0];
    }
}

#endif //NSGT_EXAMPLE_UTILS_H
