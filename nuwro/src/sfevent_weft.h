#ifndef _sfevent_weft_h_
#define _sfevent_weft_h_
#include "event1.h"
#include "nucleus.h"
#include "params.h"
#include <map>
#include <random>

class CSpectralFunc;

double sfevent_weft(params &par, event &e, nucleus &t);
static inline vect on_shell_4vec(particle& N, double M);

#endif