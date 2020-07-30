//Created by M.Bazhenov 1997-2009
//All rights reserved
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define time2(a) (2)

#ifndef _CURRENTS_H_
#define _CURRENTS_H_

#define PS 1 //Plasticity switch. 1 is ON and 0 is OFF
#define HHH      0.02   

// testing periods
extern int test1_start, test1_end, test2_start, test2_end, test3_start, test3_end, test4_start, test4_end;
extern int test5_start, test5_end, test6_start, test6_end, test7_start, test7_end, test8_start, test8_end;
extern int test9_start, test9_end, test10_start, test10_end, test11_start, test11_end, test12_start, test12_end;
extern int test13_start, test13_end, test14_start, test14_end, test15_start, test15_end;

// training periods
extern int seq1_trainingStart, seq1_trainingEnd, seq2_trainingStart, seq2_trainingEnd, seq3_trainingStart, seq3_trainingEnd;

//awake to stage 3
extern int awake_end;
extern int stage3_end;

/// CX neuron -gKl
extern double fac_gkl; 
extern double fac_gkl_TC; 
extern double fac_gkl_RE; 

extern double fac_gk_oI; 
extern double fac_gkv_cx; 
extern double fac_gkca_cx;
extern double fac_gkm_cx;

extern double fac_gh_TC; 

extern double fac_AMPA_D2;
extern double fac_AMPA_TC;

extern double fac_GABA_D2;
extern double fac_GABA_TC;

////////

class BaseCell{

 public:

  
  virtual double get_v_dend() =0;
  virtual double get_v_soma() =0;
  virtual void calc(double,double, double*, double*) = 0; 

  BaseCell(){
  }
};

class OtherBase:public BaseCell{
 public:
  virtual double get_s_cx_dend()=0;
  OtherBase():BaseCell(){
  }
};

//Identifiers for cell types
enum Cell_Type{
  E_RE = 0,
  E_REa = 1,
  E_TC = 2,
  E_TCa =3,
  E_CX = 4,
  E_CXa = 5,
  E_CX6 = 6,
  E_IN = 7,
  E_INa = 8,
  E_IN6 = 9
};

//Synaptic types
enum Syn_Type{
  E_AMPA = 11,
  E_GABA_A = 12,
  E_GABA_B = 13,
  E_AMPA_D2 = 14,
  E_NMDA_D1 = 15,
  E_GABA_A_D2 = 16,
  E_AMPA_D1 = 17,
  E_AMPA_D3 = 18,
  E_GAP = 19,
  E_AMPAMap_D1 = 20,
  E_NMDAMap_D1 = 21,
  E_GABAMap_A_D2 = 22,
  E_GiantAMPA=23,
  E_GiantAMPAMap=24,
  ENUM_Syn_Types=25 
};

typedef struct{
  double strength;
  double mini_s;
  double mini_fre;
  enum Cell_Type to_type;
  enum Cell_Type from_type;
  int to_x;
  int to_y;
  int from_x;
  int from_y;  
  enum Syn_Type type;
}Syn_Info;

class GiantSyn;
int get_cell_index(int type, int m, int n);

class Syn{

 public:

  double strength;
  double mini_s;
  double mini_fre;
  enum Cell_Type to_type;
  enum Cell_Type from_type;
  int from_cell; //from which cell we receive input
  enum Syn_Type type;
  double lastspike;
  GiantSyn *giant; //accumulator for spiking data, one per cell. NULL if not present for the particular synapse

  Syn(Syn_Info *syn_info, GiantSyn *g=NULL){ 
    this->strength = syn_info->strength;
    this->to_type = syn_info->to_type;
    this->from_type = syn_info->from_type;
    this->from_cell = get_cell_index(syn_info->from_type,syn_info->from_x,syn_info->from_y);
    this->type = syn_info->type;
    this->mini_s = syn_info->mini_s;
    this->mini_fre = syn_info->mini_fre;
    this->giant = g;

    delete syn_info;
  }

  virtual double calc(double,double,int){return 0;};

};

// --- Specific classes ----


class RS:public OtherBase{
 double xp, xpp, mu, sigma_n, beta_n;
 double sigma_dc, beta_dc, Tcr;
 double yLimit;
 
 public:
 double x, y, alpha, sigma, sigma_e, beta_e, Idc, S_CX_DEND,v_DEND,v_SOMA;
 int spike;

 RS():OtherBase(){ 
   init();
 }

 double get_v_soma(){
   return this->v_SOMA;
 }
 double get_v_dend(){
   return this->v_DEND;
 }

 double get_s_cx_dend(){
   return this->S_CX_DEND;
 }

 void init(){ 
   mu=0.0005;
   spike=0;
   alpha=3.65;
   sigma=0.02;
   sigma_e=1.0;
   sigma_dc=1.0;
   beta_e=0.133;
   beta_dc=0.133;
   Idc = 0;
   S_CX_DEND = 165.0e-6;
   xpp=-1+sigma;
   xp=-1+sigma;
   x=-1+sigma;
   y= x - alpha/(1-x);
   Tcr=-100;
   v_DEND=-60;
   v_SOMA=-60;
 }

 void calc(double, double, double*, double*); 
};  

class FS1:public OtherBase{
 double xp, xpp, mu, sigma_n, beta_n, ii, gg; 
 double sigma_dc, beta_dc, beta_ii, sigma_ii, Tcr; 
public: 
 double x, y, alpha, sigma, sigma_e, beta_e, Idc, S_CX_DEND,v_DEND, v_SOMA; 
 int spike; 

 FS1():OtherBase(){  
       init();
       } 

 double get_v_soma(){
   return this->v_SOMA;
 }
 double get_v_dend(){
   return this->v_DEND;
 }
 
 double get_s_cx_dend(){
   return this->S_CX_DEND;
 }

 void init(){  
   ii=0;
   mu=0.002;
   spike=0; 
   alpha=3.8;     
   sigma=-1.75E-2;
   sigma_e=0.5; ; 
   sigma_dc=1.0; 
   beta_e=0.05; 
   beta_dc=0.1;
   Idc = -0.; 
   S_CX_DEND = 165.0e-6; 
   xpp=-1+sigma; 
   xp=-1+sigma;
   x=-1+sigma; 
   y= x - alpha/(1-x); 
   gg = 0.5;
   beta_ii = 0.5;
   sigma_ii = 0.5;
   Tcr=-100;
   v_DEND=-60;
   v_SOMA=-60;
       } 
 void calc(double, double, double*, double*);  
};

//---------------Low-threshold Ca2+ current (RE cell)---------------------
class IT_RE {
 static double Shift, Ca_0, Cels;
 double m_inf, tau_m, h_inf, tau_h, ratio, eca, Phi_m, Phi_h, eca0;
public:
 double iT, m0, h0, Qm, Qh;
 double G_Ca;
 IT_RE(double v) {
   G_Ca = 1.75;
   Qm = 5;  
   Qh = 3; 
   Phi_m = pow(Qm,((Cels-24)/10));
   Phi_h = pow(Qh,((Cels-24)/10));
   m0 = 1/(1 + exp(-(v + Shift + 50)/7.4));
   h0 = 1/(1 + exp((v + Shift + 78)/5));
   eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);
   } 
 void calc(double m, double h, double &fm, double &fh, 
           double v, double cai, double x);
};

//--------------fast Na and K current (RE and TC cells)------------------
class INaK {
 static double Cels;
 double Alpha1, Beta1, Alpha2, Beta2, Alpha3, Beta3, v2, v2K, Phi;
 double tau_m, m_inf, tau_h, h_inf, tau_n, n_inf;
public:
 static double E_Na, E_K;
 double iK, iNa, m0, h0, n0;
 double G_Na, G_K, Vtr, VtrK;
 INaK(double v) {
   G_K = 10;///////////////////////
   G_Na = 100;/////////////////////
   Vtr = -50;
   VtrK = -50;
   v2 = v - Vtr;
   v2K = v - VtrK;
   Phi = pow(3,((Cels-36)/10));
   Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
   Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
   m0 = Alpha1/(Alpha1 + Beta1);

   Alpha2 = 0.128*exp((17 - v2)/18);
   Beta2 = 4/(exp((40 - v2)/5) + 1);
   h0 = Alpha2/(Alpha2 + Beta2);

   Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
   Beta3 = 0.5*exp((10 - v2K)/40);
   n0 = Alpha3/(Alpha3 + Beta3);     } 
 void calc(double m, double h, double n, double &fm, double &fh, double &fn, 
           double v, double x);
};

//------------------Ca-dynamics------------------------------------
class ICa {
 static double Ca_inf, K_T, K_d;
 double drive, drive0;                                 
public:
 double Taur, D;
 ICa() {Taur = 5; D = 1.; 
        drive0 = 10.0/(2.*96489.); }
 void calc(double cai, double &fcai, double iT, double x);
};


//------------------Low-theshold Ca2+ current (TC cell)-----------------
class IT_TC {
 static double Ca_0, Cels, Qm, Qh, Shift;
 double m_inf, tau_m, h_inf, tau_h, Phi_h, Phi_m, ratio, eca, eca0; 
public:
 double iT, m0, h0;
 double G_Ca;
 IT_TC(double v) {
    G_Ca = 2; 
    Phi_m = pow(Qm,((Cels-24)/10));
    Phi_h = pow(Qh,((Cels-24)/10));
    m0 = 1 / (1+exp(-(v+59)/6.2));
    h0 = 1 / (1+exp((v+83)/4.0));
    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489);  } 
 void calc(double m, double h, double &fm, double &fh,  
           double v, double cai, double x);
};

//----------------- h-current (TC cell) -----------------------------------
class Ih_TC {
 static double E_h, Shift, Cels, k2, nca, nexp, taum;
 double h_inf, tau_s, alpha, beta, k3p, cc, Phi;
public:
 double ih, p10, o10, o20;
 double G_h, k1ca, ginc, cac, pc, k4;

 Ih_TC(double v, double cai) {
    G_h = 0.02; 
    ginc = 1.5;  
    cac = 0.0015;
    pc = 0.01;
    k4 = 0.001;
    Phi = pow(3,((Cels-36)/10));
    h_inf = 1/(1 + exp((v + 75 - Shift)/5.5));
    tau_s = (taum + 1000 / (exp((v + 71.5 - Shift)/14.2) + 
                         exp(-(v + 89 - Shift)/11.6))) / Phi;
    alpha = h_inf/tau_s;
    beta = (1 - h_inf)/tau_s;
    p10 = 1/(1 + pow((cac/cai),nca));
    o10 = 1/(1 + beta/alpha + pow((p10/pc),nexp));
    o20 = pow((p10/pc),nexp) * o10;
 }
 void calc(double o1, double p1, double o2,  
                double &fo1, double &fp1, double &fo2,  
                double v, double cai, double x);
};

//----------------------Potassium A-current (TC cell)------------------------
class IA_TC {
 static double E_K, Cels;     
 double m_inf, tau_m, h_inf, tau_h, Tad;                                 
public:
 double iA, m0, h0, G_A;
 IA_TC(double v) {
   G_A = 1; 
   Tad = pow(3,((Cels-23.5)/10));
   m0 = 1.0 / (1+exp(-(v+60)/8.5));
   h0 = 1.0/(1+exp((v+78)/6)); } 
 void calc(double m, double h, double &fm, double &fh, double v, double x);
};

//---------------------Hight-threshold Ca2+ current (CX cell)----------------
class IHVA_CX {
 static double Shift, Ca_0, Cels, Qm, Qh, E_Ca;
 double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
 double ratio, eca0, eca;
public:
 double iHVA, m0, h0;
 double G_HVA;
 IHVA_CX(double v) {
   G_HVA = 0.03;
   Phi_m = pow(Qm,((Cels-23)/10));
   Phi_h = pow(Qh,((Cels-23)/10));
   vm = v + Shift;
   a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
   b = 0.94*exp((-75-vm)/17);
   m0 = a/(a+b);
   a = 0.000457*exp((-13-vm)/50);
   b = 0.0065/(exp((-vm-15)/28) + 1);
   h0 = a/(a+b);
   } 
 void calc(double m, double h, double &fm, double &fh, double v, double cai, double x);
};


//--------------Ca-dependent potassium current (CX cell)-----------------------
class IKCa_CX {
 static double E_KCa, Ra, Rb, Cels, Q, caix;     
 double m_inf, tau_m, Tad, a, b;                                 
public:
 double iKCa, m0;
 double G_KCa;
 IKCa_CX(double cai) {
   G_KCa = 0.3; 
   Tad = pow(Q,((Cels-23)/10));
   a = Ra * cai;  
   b = Rb;
   m0 = a/(a+b);
 }
 void calc(double m, double &fm, double v, double cai, double x);
};


//--------------------Potassium M-current (CX cell)-------------------------
class IKm_CX {
 static double E_Km, Ra, Rb, Cels, Q, tha, qa;     
 double m_inf, tau_m, Tad, a, b;                                 
public:
 double iKm, m0;
 double G_Km;
 IKm_CX(double v) {
   G_Km = 0.01;
   Tad = pow(Q,((Cels-23)/10));
   a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
   b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
   m0 = a/(a+b);
 }
 void calc(double m, double &fm, double v, double x);
};

//--------------------Fast potassium current (CX cell)-------------------
class IKv_CX {
 static double Ra, Rb, Cels, Q, tha, qa;     
 double m_inf, tau_m, Tad, a, b;                                 
public:
 static double E_Kv;
 double iKv, g_Kv, m0;
 double G_Kv;
 IKv_CX(double v) {
   G_Kv = 150; 
   Tad = pow(Q,((Cels-23)/10));

   a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
   b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
   m0 = a/(a+b);
 }
 void calc(double m, double &fm, double v, double x);
};


//----------------Fast sodium current (CX cell)--------------------------
class INa_CX {
 static double Shift, Ca_0, Cels, Qm, Qh;
 static double tha, qa, Ra, Rb, thi1, thi2, qi, thinf, qinf, Rg, Rd; 
 double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
 double trap0(double v, double th, double a, double q) {
       if (fabs(v/th) > 1.0e-6) {
               return ( a * (v - th) / (1 - exp(-(v - th)/q)) );
       } else {
	  return (a * q ); }
       }
public:
 static double E_Na;
 double iNa, g_Na, m0, h0;
 double G_Na;
 INa_CX(double v) {
   G_Na = 3000;
   Phi_m = pow(Qm,((Cels-23)/10));
   Phi_h = pow(Qh,((Cels-23)/10));
   vm = v + Shift;
   a = trap0(vm,tha,Ra,qa);
   b = trap0(-vm,-tha,Rb,qa);
   m0 = a/(a+b);

   a = trap0(vm,thi1,Rd,qi);
   b = trap0(-vm,-thi2,Rg,qi);
   h0 = 1/(1+exp((vm-thinf)/qinf));
   } 
 void calc(double m, double h, double &fm, double &fh, 
           double v, double x);
};


//-------------------Persist Sodium current (CX cell)-------------------
class INap_CX {
 double Tet, Sig, f, Cels, Q10, Phi_m, tau_m, m_inf;    
public:
 double m0, iNap, G_Nap, g_Nap;
 INap_CX(double v) {
   G_Nap = 2; 
   Tet = -42;
   Sig = 5;  
   f = 0.02;
   Cels = 36;
   Q10 = 2.7;
   Phi_m = pow(Q10,((Cels-22)/10));      
   tau_m = 0.8/Phi_m;
   m0 = f/(1 + exp(-(v - Tet)/Sig));
 }
 void calc(double, double&, double, double);
 double cond();
};

/// Synapse classes

class GAP:public Syn{

 public:
  GAP(Syn_Info *syn_info):Syn(syn_info)  {

  }
  double calc_gap(double time,double v_dend,double input_dend){
    return strength*(v_dend - input_dend);
  }
  double calc(double,double,int) { printf("Error: GAP::calc called\n"); exit(1);};

};

//Giant class represents all synapses of a particular type for a given neuron.
class GiantSyn{
 protected:
  int input_spikes; //input spikes during the last time step from all connected neurons
  double I;   //current activation based on the previous time step

 public:
  enum Syn_Type syn_type;
  enum Cell_Type from_type;
  virtual void add_spike() {input_spikes++;}; //successively adding spikes from all synapses
  virtual void reset_spikes() {input_spikes=0;};
  virtual double calc() = 0;
  GiantSyn():input_spikes(0),I(0) {};
};

class GiantSynAMPA: public GiantSyn{
 public:
  double alpha;     //decay constant 
  double gamma;     //conductance
  double mini_s;
  double strength;
  GiantSynAMPA(double a=1, double g=1): alpha(a), gamma(g){};
  virtual double calc(){
    return I = alpha*I + input_spikes*gamma;
  }
};

class GiantSynAMPAMap: public GiantSyn{
 public:
  double alpha;     //decay constant 
  double gamma;     //conductance
  double mini_s;
  double strength;
  GiantSynAMPAMap(double a=1, double g=1): alpha(a), gamma(g){};
  virtual double calc(){
    return I = alpha*I + input_spikes*gamma;
  }
};



//---------SYNAPCES DESCRIPTION-------------------------------------------
//---------first order kinet model for GABA-A synapse---------------------
class GABA_A :public Syn{
 static double Cdur, Cmax, Deadtime, Prethresh;  
 static double de;
 double lastrelease;
 double q, Rinf, Rtau;
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double I, E_GABA, Alpha, Beta,old_v;
 double R, C, R0, R1;
 GABA_A(Syn_Info *syn_info):Syn(syn_info)  {
   old_v = 0;
   E_GABA = -70;
   R = 0, C = 0, R0 = 0, R1 = 0;
   Alpha = 10.5;
   Beta = 0.166;
   lastrelease = -100;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
 }
 void init() {
   E_GABA = -70;
   R = 0, C = 0, R0 = 0, R1 = 0;
   Alpha = 10.5;
   Beta = 0.166;
   lastrelease = -100;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
 } 

 double calc(double x, double y_pre, int y_post);
}; 


//------second order kinet model (including G-proteins) for GABA-B synapse----
class GABA_B :public Syn{
 static double E_GABA, Cmax, Deadtime, Prethresh;   
 static double Kd, n; 
 double Gn, q;
public:
 double C, lastrelease;
 double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4,old_v;
 GABA_B(Syn_Info *syn_info):Syn(syn_info)  {
   old_v = 0;
   Cdur = 0.3; 
   K1 = 0.52; 
   K2 = 0.0013; 
   K3 = 0.098;
   K4 = 0.033;
   lastrelease = -10000000;
   C = 0, r0 = 0, g0 = 0;
 }
 void calc_gaba_b(double r, double g, double &fr, double &fg, 
		  double g_GABA_B, double x, double y_pre, int y_post);
 double calc(double,double,int) { printf("Error: GABA_B::calc called\n"); exit(1);};
};
//------------
//GB is integrator of GABA_B
class GB: public GABA_B {
public:
 GB(Syn_Info *syn_info) :GABA_B(syn_info) { }
 void init(double *y, int location){

   y[0+location] = 0;
   y[1+location] = 0;
   lastrelease = -10000000;
   C = 0;

 }   
 double calc_gb( double, double*, double*, double, int);
};


//------------first order kiner model for AMPA synapse---------------------
class AMPA :public Syn{
 int s;
 double q, Rinf, Rtau;
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
 static double Alpha, Beta;
 static double E_AMPA;
 static double de;
 double I,old_v;
 double R, C, R0, R1;
 double lastrelease, lastspike;
 static double Cconst;
 AMPA(Syn_Info *syn_info):Syn(syn_info)  {
   old_v = 0;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -100;
   lastspike = -100;
   s = 1;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
 }
 double calc(double x, double y_pre, int y_post);
};


//------------first order kiner model for AMPA synapse WITH depression------
class AMPA_D1:public Syn {
 static double E_AMPA;
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
 static double Alpha, Beta; 
 int s;
 double R, C, R0, R1;
 double lastrelease, lastspike, Use, Tr;
 double q, Rinf, Rtau;
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double I, E, old_v;
 AMPA_D1(Syn_Info *syn_info):Syn(syn_info)  {
   old_v = 0;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -1000;
   lastspike = -1000;
   s = 1;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
   E = 1;
   Use = 0.073; 
   Tr = 700;
 }
 void iii(unsigned int seek) {
 }
 double calc(double x, double y_pre, int y_post);
};



//--first order kiner model for AMPA synapse WITH depression & spont releases------
class AMPA_D2:public Syn {
 static double E_AMPA;
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
 static double Alpha, Beta; 
 int s;
 double C, R0, R1;
 double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr;
 double q, q1, Rinf, Rtau;
 double Tau, Period, S, SS;
 drand48_data *rand_buffer;  
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double R, I, E, g1, old_v;
 AMPA_D2(Syn_Info *syn_info):Syn(syn_info)  {
   rand_buffer = new drand48_data;
   static int seed = time2(NULL);
   seed = seed + 100;
   srand48_r(seed,rand_buffer);
   old_v = 0;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -10000;
   lastrelease1 = -10000;
   newrelease = 0;
   s = 1;
   g1=0.00006;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
   E = 1;
   Use = 0.07; 
   Tr = 700;
   Period = 8000;
   Tau = 50;
 }
 void iii(unsigned int seek) {
 }
 double calc(double x, double y_pre, int y_post);
};


//--------------------D3-------------
//--first order kiner model for AMPA synapce WITH depression & spont releases & Learning------                                                                                 
class AMPA_D3:public Syn {

  static double E_AMPA;
  static double Cdur, Cmax, Deadtime, Prethresh, Prethresh1, Cdel;
  static double Alpha, Beta;
  int s;
  double R, C, R0, R1;
  double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr;
  double spike_pre, spike_post, keyP, keyD, MinFactor, TauL, cP, cD;
  double q, q1, Rinf, Rtau;
  double Tau, Period, S, SS, factor;
  drand48_data *rand_buffer;  
  double exptable(double z)
    {
      if((z > -10) && (z < 10)) return( exp(z) );
      else return( 0 );
    }
 public:
  double I, E, g1, g_AMPA0, g_AMPA00;
  double weight_max,weight_min,minis0;
  AMPA_D3(Syn_Info *syn_info):Syn(syn_info)  {
    rand_buffer = new drand48_data;
    static int seed = time2(NULL);
    seed = seed+100;
    srand48_r(seed,rand_buffer);
    R = 0, C = 0, R0 = 0, R1 = 0;
    lastrelease = -10000;
    lastrelease1 = -10000;
    spike_post=-100;
    spike_pre=-100;
    newrelease = 0;
    TauL=20;
    s = 1;
    g1=0.00006;
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
    Rtau = 1 / ((Alpha * Cmax) + Beta);
    E = 1;
    Use = 0.07;
    Tr = 700; 
    Period = 8000;
    Tau = 50;
    factor = 1;
    cP=0.002;
    cD=0.002;
    keyP=0, keyD=0;
}
  void iii(unsigned int seek) {
  }
  double calc(double time, double y_post, int y_pre);
  double calc_ampad3(double time, double v_post, int y_pre,int y_post);
};


//------------first order kiner model for NMDA synapse WITH depression------
class NMDA_D1 :public Syn{
 static double E_NMDA;
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel; 
 static double Alpha, Beta; 
 int s;
 double R, C, R0, R1;
 double lastrelease, lastspike, Use, Tr;
 double q, Rinf, Rtau, fn;
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double I, E, old_v;
 NMDA_D1(Syn_Info *syn_info):Syn(syn_info) {
   old_v = 0;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -100;
   lastspike = -100;
   s = 1;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
   E = 1;
   Use = 0.0; 
   Tr = 750;
 }
 void iii(unsigned int seek) {
 }
 double calc(double x, double y_pre, int y_post);
};



//---first order kinet model for GABA-A synapse with DEPRESSION & spont IPSPs--
class GABA_A_D2 :public Syn{
 static double Cdur, Cmax, Deadtime, Prethresh; 
 static double Alpha, Beta;  
 double R, C, R0, R1;
 double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr;
 double q, q1, Rinf, Rtau;
 double Tau, Period, S, SS;
 drand48_data *rand_buffer;  
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double I, E_GABA, E, g1,old_v;
 GABA_A_D2(Syn_Info *syn_info):Syn(syn_info)  {
   rand_buffer = new drand48_data;
   static int seed = time2(NULL);
   seed = seed + 103;
   srand48_r(seed,rand_buffer);
   old_v = 0;
   E_GABA = -70;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -10000;
   lastrelease1 = -10000;
   newrelease = 0;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
   E = 1;
   Use = 0.07;
   Tr = 700;
   Period = 8000;
   Tau = 50;
 }
 void iii(unsigned int seek) {
 }
 double calc(double x, double y_pre, int y_post);
}; 



//=------------------------------
class AMPAmapD1:public Syn {
  static double E_AMPA;
  static double d_dep, d_rec;
  double d, gamma;   // depression variable
  double lastrelease, lastrelease1, newrelease;
  double Tau, S, SS, factor, Tcr, h;
  drand48_data *rand_buffer;  
 public:
  double I, g, old_v;
  AMPAmapD1(Syn_Info *syn_info):Syn(syn_info)  {
    rand_buffer = new drand48_data;
    static int seed = time2(NULL);
    seed = seed + 102;
    srand48_r(seed,rand_buffer);
    old_v = 0;
    h=HHH;
    I=0;
    g=0;
    d=1;
    lastrelease = -10000;
    lastrelease1 = -10000;
    newrelease = -10000;
    Tau = 50;
    Tcr=-100;
    gamma=exp(-HHH/2);
  }
  void iii(unsigned int seek) {
  }
 double calc(double x, double y_pre, int y_post);
};


//=------------------------------
class AMPAmapD:public Syn {
  static double E_AMPA;
  static double d_dep, d_rec;
  double d, gamma;   // depression variable
  double Tcr, h;
 public:
   double I, g, old_v;
   AMPAmapD(Syn_Info *syn_info):Syn(syn_info)  {
     old_v = 0;
     h=HHH;
     I=0;
     g=0;
     d=1;
     Tcr=-100;
     gamma=exp(-HHH/2);
   }
   void iii(unsigned int seek) {
   }
 double calc(double x, double y_pre, int y_post);
};

//=------------------------------
class AMPAmap :public Syn{
   static double E_AMPA;
   double Tcr, gamma, h;
public:
   double I, g, old_v;
   AMPAmap(Syn_Info *syn_info):Syn(syn_info) {
     old_v = 0;
       h=HHH;
       I=0;
       g=0;
       Tcr=-100;
       gamma=exp(-HHH/2);
  }
   void iii(unsigned int seek) {
   }
 double calc(double x, double y_pre, int y_post);
};


//=------------------------------
class NMDAmapD1 :public Syn{
   static double E_NMDA;
   static double d_dep, d_rec;
   double fn, d, gamma;   // depression variable
   double Tcr, h;
public:
   double I, g,old_v;
   NMDAmapD1(Syn_Info *syn_info):Syn(syn_info) {
     old_v = 0;
       h=HHH;
       I=0;
       g=0;
       d=1;
       Tcr=-100;
       gamma=exp(-HHH/50);
  }
   void iii(unsigned int seek) {
   }
 double calc(double x, double y_pre, int y_post);
};

//=------------------------------
class GABAAmapD1 :public Syn{
   double Tcr, gamma, h;
public:
   static double E_GABAA;
   double I, g,old_v;
   GABAAmapD1(Syn_Info *syn_info):Syn(syn_info) {
     old_v = 0;
     I=0;
     g=0;
     h=HHH;
     gamma=exp(-HHH/0.379);
   }
   void iii(unsigned int seek) {
   }
 double calc(double x, double y_pre, int y_post);
};

//-----first order kiner model for AMPA synapse used for external stimulation----
class Extern_ampa:public Syn{
 static double Cdur, Cmax, Deadtime, Prethresh; 
 double R, C, R0, R1;
 double lastrelease;
 double q, Rinf, Rtau;
 double TR, w, wom, RRR;
 double exptable(double z)
   {
   if((z > -10) && (z < 10)) return( exp(z) );
   else return( 0 );
   }
public:
 double g, Alpha, Beta,old_v;
 Extern_ampa(Syn_Info *syn_info):Syn(syn_info) {
   old_v = 0;
   Alpha = 0.94;
   Beta = 0.18;
   R = 0, C = 0, R0 = 0, R1 = 0;
   lastrelease = -100;
   Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
   Rtau = 1 / ((Alpha * Cmax) + Beta);
   TR = 1000, w=0.01, wom=0;
 }
 void iii(unsigned int seek) {
 }
 void calc(double g_Extern_ampa, double x);
};

//-------------------basic RE CELL-----------------------------------------------
class RE: public BaseCell, public IT_RE, public INaK, public ICa {// public IKCa_RE, public ICAN_RE,
 static double Cai0, V0;
public:
 double G_kl, G_l, E_l, S_RE,v_SOMA,v_DEND;

 RE() :BaseCell(),IT_RE(V0), INaK(V0), ICa() {
       E_l = -70;
       G_l = 0.05;
       G_kl = 0.018;
       S_RE = 1.43e-4;
       v_DEND = v_SOMA = V0;
       }

 double get_v_soma(){
   return this->v_SOMA;
 }
 double get_v_dend(){
   return this->v_DEND;
 }


/* Reduced model modif */
 void init(double *y) {

   y[0] = V0;
   y[1] = Cai0;
   y[2] = IT_RE::m0;
   y[3] = IT_RE::h0;
   y[4] = INaK::m0;
   y[5] = INaK::h0;
   y[6] = INaK::n0;

 }
void calc(double , double, double* , double*); 
};  

//-------------------basic TC CELL-----------------------------------------------
class TC: public BaseCell, public IT_TC, public Ih_TC, public INaK, public ICa, public IA_TC {
 static double G_l, Cai0, V0;
public:
 double G_kl, S_TC, E_l, DC,v_SOMA,v_DEND;

 TC() :BaseCell(),IT_TC(V0),Ih_TC(V0,Cai0), INaK(V0), ICa(), IA_TC(V0) {
       G_kl = 0.01;
       E_l = -70;
       DC = 0;
       INaK::G_Na = 90;
       INaK::G_K = 10;
       S_TC = 2.9e-4;
       v_DEND = v_SOMA = V0;
       }

 double get_v_soma(){
   return this->v_SOMA;
 }
 double get_v_dend(){
   return this->v_DEND;
 }

 void init(double *y) {

   y[0] = V0;
   y[1] = Cai0;
   y[2] = IT_TC::m0;
   y[3] = IT_TC::h0;
   y[4] = Ih_TC::o10;
   y[5] = Ih_TC::p10;
   y[6] = Ih_TC::o20;
   y[7] = INaK::m0;
   y[8] = INaK::h0;
   y[9] = INaK::n0;
   y[10] = IA_TC::m0;
   y[11] = IA_TC::h0;

 }
void calc(double, double, double*, double*); 
};  


//-------------------CX CELL (DENDRITE)-------------------------------------
class CX_DEND: public IHVA_CX, public IKCa_CX, public IKm_CX, 
              public INa_CX, public INap_CX, public ICa {
public:
  double iDEND, I_Stim1, E_l, G_kl, G_l;
  CX_DEND(double V0, double Cai0) :IHVA_CX(V0), IKCa_CX(Cai0), IKm_CX(V0),
                                  INa_CX(V0), INap_CX(V0), ICa() { 
    E_l = -70;
    G_kl = 0.005; 
    INa_CX::G_Na = 0.8;
    I_Stim1 = 0; 
    ICa::Taur = 165;
    INap_CX::G_Nap = 2.5;
        
    IKCa_CX::G_KCa = 0.05;
    IHVA_CX::G_HVA = 0.01;;

    IKm_CX::G_Km = 0.015; 
    G_l = 0.02;//  mS/cm^2 

  }
  void init(double *y, double V0, double Cai0) {

    y[0] = V0;
    y[1] = Cai0;
    y[2] = IHVA_CX::m0;
    y[3] = IHVA_CX::h0;
    y[4] = IKCa_CX::m0;
    y[5] = IKm_CX::m0;
    y[6] = INa_CX::m0;
    y[7] = INa_CX::h0;
    y[8] = INap_CX::m0;
  }
  void calc(double, double*, double*); 
};  

//-------------------CX CELL (SOMA)-------------------------------------
class CX_SOMA: public IKv_CX, public INa_CX, public INap_CX {
public:
  double v_SOMA, iSOMA, g1_SOMA, g2_SOMA, I_Stim2;
  CX_SOMA(double V0, double Cai0) :IKv_CX(V0), INa_CX(V0), INap_CX(V0){ 
    I_Stim2 = 0; 
    IKv_CX::G_Kv = 200; 
    INa_CX::G_Na = 3000; 
    INap_CX::G_Nap = 15;
  }
  void init(double *y, double V0, double Cai0, int start_location) {
    v_SOMA = V0;
    y[0+start_location] = IKv_CX::m0;
    y[1+start_location] = INa_CX::m0;
    y[2+start_location] = INa_CX::h0;
    y[3+start_location] = INap_CX::m0;
  }
 void calc(double, double*, double*); 
}; 

//------------CX CELL (connect DENDRITE and SOMA)---------------------------
class CX: public OtherBase, public CX_DEND, public CX_SOMA {
 static double Cai0, V0, C;
 int N_DEND;
public:
 double field_effect_dend, axil_current,field_effect_soma;
 double kappa, rho, S_CX_SOMA, S_CX_DEND, v_DEND; 
 CX() :OtherBase(),CX_DEND(V0,Cai0), CX_SOMA(V0,Cai0) {
       kappa = 10.0e3;      // kOm: to get mS=1/kOm 
       rho = 165; 
       field_effect_dend = 0.0;
       field_effect_soma = 0.0;
       axil_current = 0.0;
       v_DEND = v_SOMA = V0;
       }

 double get_s_cx_dend(){
   return this->S_CX_DEND;
 }

 double get_v_soma(){
   return this->v_SOMA;
 }
 double get_v_dend(){
   return this->v_DEND;
 }

 void init(double *y, int n_dend) {
       N_DEND=n_dend;
       CX_DEND::init(y,V0, Cai0);
       CX_SOMA::init(y,V0, Cai0,n_dend);
       S_CX_SOMA = 1.0e-6;  // cm^2
       S_CX_DEND = S_CX_SOMA * rho;   
       }
void calc(double, double, double*, double*); 
};  
#endif
