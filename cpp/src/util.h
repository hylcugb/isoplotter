#pragma once

#include <stdio.h>

#define db(x) dbf("%s", x)
#define dbf(x...) printf(x); printf("\n")
#define err(msg...) {fprintf(stderr, msg); fprintf(stderr, "\n"); exit(1);}
#define errif(expr, msg...) if(expr) { err(msg); }
#define errif_(expr) errif(expr, #expr)

bool equals(float a, float b, float epsilon);
