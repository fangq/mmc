#ifndef _MMC_FAST_MATH_FUNCTIONS_H
#define _MMC_FAST_MATH_FUNCTIONS_H

static inline float fast_expf9(float x) {
        return (362880.f+x*(362880.f+x*(181440.f+x*(60480.f+x*(15120.f+x*(3024.f+x*(504.f+x*(72.f+x*(9.f+x)))))))))*2.75573192e-6f;
}

#endif
