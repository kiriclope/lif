#ifndef __GLOBALVARS__
#define __GLOBALVARS__

// Simulation time constants
#define dt .1
#define duration 10.E3 

#define Tst 1.0E3 
#define Tw 1.E3 
#define Tl 100.0E3 

#define nbPref 10000.0 
/* #define nbPref 10971.5 */ 
#define IF_Nk 0 

// LIF parameters 
#define Vr -70. // Membrane Resting Potential 
#define Vth -50. // Voltage Threshold 
#define Vpeak 20. // Spikes Peak 

#define m0 20.0E-2 // External rate c 
// Membrane time constants 
const double Tm[4] {20.,10.,20.,20.} ; 

#define argIext 0 

#define IF_JabLoop 0 
#define Ax 0 
#define By 0 

#define IF_Prtr 0
#define IF_GAUSS 0 
#define PrtrPop 1 

#define IF_RING 0 
#define L 1.5 
#define DIM 1
#define IF_SPEC 0 

#define PULSE 40 
#define IF_PULSE 0 

#define IF_ADAPTATION 0 
#define Tad 100. 
#define Gad 1. 

#define PROFILE 2 

#define IF_OPSIN 0 
#define OpsPb 1.

#define IF_TIMECOURSE 0 
#define Tc 0.E3 
#define dTc 0.5E3 
#define dPrtr 0.2
#endif
