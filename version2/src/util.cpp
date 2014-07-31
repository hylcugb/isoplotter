#include "util.h"

#include <math.h>

bool equals(float a, float b, float epsilon) {
    return fabs(a - b) < epsilon;
}
