#ifndef __GLOBALVARS__ 
#define __GLOBALVARS__ 

// Simulation time constants
#define DT .01
#define duration 10.0E3 // 10.E3  //

#define Tst 10.0E3 //  2.0E3 // 10.E3 // 
#define Tw .250E3 // 1.00E3 // 10.E3 //  
#define Tl 2.0E3 // 2.0E3 // 30.0E3 // 

#define nbPref 10000
#define IF_Nk 1

#define IF_BENCHMARK 1 
#define IF_INTERPOLATION 1
#define IF_EULER 1 
#define IF_RK2 0

// LIF parameters 
#define Vr -70. // Membrane Resting Potential 
#define Vth -50. // Voltage Threshold 
#define Vpeak 20. // Spikes Peak 

#define m0 1.0E-2 // External rate 
// Membrane time constants 
const double Tm[4] = {20.,10.,20.,20.} ; 
#define argIext 0 

#define IF_TRANSIENT_IEXT 0 
#define T_IEXT 1.0E3 
#define nbIext 1.0E3 

#define IF_JabLoop 0
#define Ax 0 
#define By 0 

#define IF_Prtr 0
#define IF_GAUSS 0 
#define PrtrPop 1

#define IF_RING 0 
#define L 1.0 
#define DIM 2 
#define IF_SPEC 0 
#define PROFILE 2 // 2 
#define PHI0 0. 

#define PULSE 40 
#define IF_PULSE 0 

#define IF_ADAPTATION 0 
#define Tad 100. 
#define Gad 1. 

#define IF_OPSIN 0 
#define OpsPb 1.

#define IF_TIMECOURSE 0
#define Tc 1.E3 
#define dTc 1.E3 
#define dPrtr .25 

#define IF_SHARED 0
#define PROP_SHARED .25
#define CLUSTER_SIZE .5

#endif
