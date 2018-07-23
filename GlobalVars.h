#ifndef __GLOBALVARS__
#define __GLOBALVARS__

// Simulation time constants
#define dt .1
#define duration 10.0E3
#define Tst 0.0E3
#define Tw .10E3
#define Tl 4.0E3

// LIF parameters
#define Vr  0. // Membrane Resting Potential
#define Vth  1. // Voltage Threshold
#define Vpeak 20. // Spikes Peak

#define m0 1.0E-2 // External rate
const double Tm[4] = {20.,10.,10.,10.} ; // Membrane time constants

#define IF_Nk 0

#define argIext 0

#define IF_Prtr 0
#define PrtrPop 1
#define IF_GAUSS 0

#define IF_RING 0
#define L 2.0*M_PI
#define IF_SPEC 0

#define PULSE 40
#define IF_PULSE 0

#define IF_ADAPTATION 0
#define Tad 100.
#define Gad 1.

#define PROFILE 2

#define IF_OPSIN 0
#define OpsPb 1.

#endif