#ifndef CellSyn_h
#define CellSyn_h

#include "currents.h"    //header for all classes describing currents and basic cell types
#include <list>

using namespace std;

//it is important this comes before the other includes
#define I_HH    1     // 1 -- HH model; 0 -- for Maps;

//------------Number of ODE for each cell -------------------------------

#define N_RE 7
#define N_TC 12 
#define N_GB 2
#define N_TCa 12 

#if (I_HH==1)
#define N_DEND   9
#define N_SOMA   4 
#else 
#define N_DEND   0
#define N_SOMA   0
#endif

#define N_CX     (N_DEND + N_SOMA)
#define N_CXa     (N_DEND + N_SOMA)
#define N_CX6     (N_DEND + N_SOMA)
#define N_IN     N_CX   
#define N_INa     N_CX 	
#define N_IN6     N_CX 	

//other defines
#define TAU      0.02 //integration time step
#define TAU_Map  0.5  // time step for maps
#define TAUr     25   //(TAU_Map/TAU)

void print_spike_time(double time,int type,int m,int n);
extern double input_strength, input_start_neuron, input_end_neuron;
extern int seq1_trainingStart, seq1_trainingEnd, seq2_trainingStart, seq2_trainingEnd, seq3_trainingStart, seq3_trainingEnd;
extern int test1_start, test1_end, test2_start, test2_end, test3_start, test3_end, test4_start, test4_end;
extern int test5_start, test5_end, test6_start, test6_end, test7_start, test7_end, test8_start, test8_end;
extern int test9_start, test9_end, test10_start, test10_end, test11_start, test11_end, test12_start, test12_end;
extern int test13_start, test13_end, test14_start, test14_end, test15_start, test15_end;
extern int awake_end,stage3_end;

//information about each cell as read from network file. CellSyn contain members relater rather to computation
typedef struct{
  int *num_connects_norm; //number of incoming dendrites of each type
  int *num_connects_gap; //number of gap connections
  int *num_connects_total;

  //same as 3 above but does not include repeated connections involving the same two cells
  int *num_cell_cons; //includes every thing including short range and gap 
  int *num_cell_cons_gap; //does not include short range
  int *num_cell_cons_norm; //does not include short range

  int total_connects; //total numper of incoming 
  int cell_index; //the index identifing where in "cells" array this cell is 
  Syn_Info **syn_info;
}Cell_Info;

typedef struct{
  int x;
  int y;
}Pair;



//abtract class for the combination of a cell and its synapsis
//7 classes derive from it. 
class CellSyn{

 protected:
  int lhs_size; // how big y,yk,f and s are
  int num_b_eq; //used to decide how far into y an f you have to go to find the indexes for I_GB
  double  *yk; 
  double *s; 
  double cur_time;
  double xk;
  double tau; //time step size
  int current_step; //which of the 4 integration steps we are on

  Cell_Info *cell_info; //information about the cell sent from the neuro file

  double AMPA_Rsum; //Partial sums for R variables in synapse
  //list of giant synapses for this cell.
  //individual synapses are responsible for having proper pointer to this list
  list<GiantSyn*> giant_syns;
  //find or create new giant synapse in case it was not created yet
  GiantSyn* get_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //clear spikes from last time step
  void reset_GiantSyns();
  //call calcs for all giant syns and return sum of currents
  double calc_GiantSyns();

 public:

  Syn **syns; //synapses
  int num_syns; // number of synapses
  drand48_data *rand_buffer;
  double *y, *f;
  int m; // cell's location
  int n; // cell's location
  int flip; // if cell is spiking
  int flip_for_maps; // check if the cell is spiking for maps
  double old_v; // voltage before last calc of cell
  int type; // type of this cell
  int num_spikes; // how many spikes this window
  int ismap; // map or not

  int cell_type;
  int x_index;

  BaseCell *base_cell; //Either an RE,CX, or TC cell defined in currents

  CellSyn(int i,int j,int b_size,int t_size,Cell_Info *cell_info,int type);

  //this function creates all the cells(CellSyns) and does some basic initialization
  //all arguments are inputs.
  static CellSyn* initialize(int type,int i, int j, Cell_Info ****cells_info, string output_location);

  double zero_one_rand(){
    double my_rand;
    my_rand = rand()/(RAND_MAX + 1.0);
    return my_rand;
  }
  int get_flip(){
    return this->flip;
  }
  int get_flip_for_maps(){
    return this->flip_for_maps;
  }

  //manages spiking event (signal & store needed info)
  void signal_spike();
  void reset_flip_maps();

  //RK computation
  void zero_temps();
  void reset_y();
  void step();

  //step - 1 for first step of runge-kutta 0 for all other steps
  void calc(double x, double *y_ini, double *f_ini, int step);
  //create synapses and init base cell
  virtual void memory(Syn_Info **syn_info,int num_syns) = 0;
  //each type of neuron needs to define its own matrix of used Giant synapses
  //we ditched virtuality and use generic version for all from->to types right now
  //otherwise those need to be defined in particular neurons types and used in ::memory()
  GiantSyn* define_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //convenience function for searching
  GiantSyn* find_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //note that after calling this, syn_info parameter is already freed from memory
  Syn* initialize_synapses(Syn_Info *syn_info);
};





///////////////////////////////////////////////////////////////=================================
//we now define several classes that are derived from CellSyn.
class REsyn:public CellSyn{

 public:

 REsyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new RE();
    //type = E_RE;
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    ((RE*)base_cell)->E_l = -77; 
    ((RE*)base_cell)->G_Ca = 2.2;
    ((RE*)base_cell)->G_kl = 0.012;
    this->AMPA_Rsum=0;
    this->ismap=0;

    this->cell_type=E_RE;
    this->x_index=i;
  }

  void memory(Syn_Info **syn_info,int num_syns){

    ((RE*)base_cell)->init(y);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength/((RE*)base_cell)->S_RE;

      if(syns[i]->type == E_GiantAMPAMap){
        ((GiantSynAMPA*) syns[i]->giant) -> alpha = 1;
      }
      if(syns[i]->to_type == E_RE){
        ((GABA_A*)syns[i])->E_GABA = -70;
      }
    } 
  }

};

//-------------------TC CELL core type and all Synapses from other cells-----------------------------------------------
class TCcore: public CellSyn {

 public:

 TCcore(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new TC();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;

    this->cell_type=E_TC;
    this->x_index=i;

    ((TC*)base_cell)->G_A = 0;
    ((TC*)base_cell)->ginc = 2.0; 
    ((TC*)base_cell)->G_Ca = 2.5; 
   
    ((TC*)base_cell)->D = 2;
    ((TC*)base_cell)->pc = 0.007;
    ((TC*)base_cell)->k4 = 0.001;
    ((TC*)base_cell)->Vtr = -40;
    ((TC*)base_cell)->VtrK = -28;
    ((TC*)base_cell)->G_K = 12;
    ((TC*)base_cell)->G_h = 0.016;  
  
    ((TC*)base_cell)->G_kl = 0.024;
    ((TC*)base_cell)->DC = 0;
    
  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((TC*)base_cell)->init(y);
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_AMPA){
      }else if(syns[i]->type == E_GABA_A){
        ((GABA_A*)syns[i])->E_GABA = -83;
      }else if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->K1 = 0.5; 
        ((GB*)syns[i])->K2 = 0.0012;
        ((GB*)syns[i])->K3 = 0.1; 
        ((GB*)syns[i])->K4 = 0.034;
        ((GB*)syns[i])->init(y,12+2*k);
        k = k + 1;
      }
      syns[i]->strength = syns[i]->strength/((TC*)base_cell)->S_TC;
    } 
  }
  
};  

//-------------------TC CELL matrix type (TCa) and all Synapses from other cells-----------------------------------------------
class TCmatrix: public CellSyn {

 public:

 TCmatrix(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new TC();
    //type = E_TCa;
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;

    this->cell_type=E_TCa;
    this->x_index=i;

    ((TC*)base_cell)->G_A = 0;
    ((TC*)base_cell)->ginc = 2; 
    ((TC*)base_cell)->G_Ca = 2.2;
   
    ((TC*)base_cell)->D = 2;
    ((TC*)base_cell)->pc = 0.007;
    ((TC*)base_cell)->k4 = 0.001;
    ((TC*)base_cell)->Vtr = -40;
    ((TC*)base_cell)->VtrK = -28;
    ((TC*)base_cell)->G_K = 12;
    ((TC*)base_cell)->G_h = 0.015;  
  
    ((TC*)base_cell)->G_kl = 0.025; 
    ((TC*)base_cell)->DC = 0; 

  }

  void memory(Syn_Info **syn_info,int num_syns){

    ((TC*)base_cell)->init(y);
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_A){
        ((GABA_A*)syns[i])->E_GABA = -83; 
      }else if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->K1 = 0.5; 
        ((GB*)syns[i])->K2 = 0.0012;
        ((GB*)syns[i])->K3 = 0.1; 
        ((GB*)syns[i])->K4 = 0.034;
        ((GB*)syns[i])->init(y,12+2*k);
        k = k + 1;
      }

      syns[i]->strength = syns[i]->strength/((TC*)base_cell)->S_TC;
    } 
  }
};  


//Below we describe first headers for specific cortical cell classes based on the HH type basic class and then we describe very similar headers for the same cortical classes but based on map type basic class. The switch I_HH select which headers are used. The functions associated with these specific cell classes are the same and are describe dbelow
//-------------------CX CELL and all Synapses from other cells-----------------------------------------------
class CXsyn: public CellSyn {

 public:
  
 CXsyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
 
    this->cell_type=E_CX;
    this->x_index=i;

    ((CX*)base_cell)->G_kl = 0.011; 
    ((CX*)base_cell)->E_l = -67.0; 
    ((CX*)base_cell)->G_Km = 0.02; 
    ((CX*)base_cell)->G_l = 0.011 + (((double) rand() / (RAND_MAX)) + 1) * 0.003;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
          
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

  void print_c_stren(FILE *connection_record) { 
        
    double total = 0.0;
    int total_num = 0;
    int i = 0;

    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        total = total + ((AMPA_D2*)syns[i])->strength; 
        total_num = total_num + 1;
      }
    }    
    double average = total / total_num;
    fprintf(connection_record,"%lf %lf %d %d",cur_time,average,m,n);
    
    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){ 
        fprintf(connection_record," %lf",((AMPA_D2*)syns[i])->strength); 
      }
    }    
    fprintf(connection_record,"\n");
  }
  

};  

//-------------------CXa CELL and all Synapses from other cells------------------
class CXasyn: public CellSyn {

 public:

 CXasyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    ((CX*)base_cell)->G_kl = 0.01; 
    ((CX*)base_cell)->E_l = -68; 
    ((CX*)base_cell)->G_Km = 0.02; 
    ((CX*)base_cell)->G_l =  0.022; 
 
    this->cell_type=E_CXa;
    this->x_index=i;

    this->ismap=0;
    this->AMPA_Rsum=0;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);
    
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

 
    }
  }


};  

//CX6 layer
//-------------------Cx6 CELL and all Synapses from other cells-----------------------------------------------
class CXsyn6: public CellSyn {

 public:

 CXsyn6(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);

    this->cell_type=E_CX6;
    this->x_index=i;

    ((CX*)base_cell)->G_kl = 0.01; 
    ((CX*)base_cell)->E_l = -68;
    ((CX*)base_cell)->G_Km = 0.02; 
    ((CX*)base_cell)->G_l =  0.022;
    this->ismap=0;
    this->AMPA_Rsum=0;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

// TB
//-------------------IN CELL core type and all Synapses from other cells-----------------------------------------------
class INsynCore: public CellSyn {
  
 public:  

 INsynCore(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;
    
    ((CX*)base_cell)->G_kl = 0.009; 
    ((CX*)base_cell)->G_l = 0.009 + (((double) rand() / (RAND_MAX)) + 1) * 0.003;

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    this->cell_type=E_IN;
    this->x_index=i;
    
   
  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();



    }
  }
};  

//-------------------IN CELL matrix type (INa) and all Synapses from other cells---------
class INsynMatrix: public CellSyn {

 public:

 INsynMatrix(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new CX();
    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    this->cell_type=E_INa;
    this->x_index=i;

    
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 

//-------------------IN CELL matrix type (INa) and all Synapses from other cells--------
class INsyn6: public CellSyn {
  
 public:

 INsyn6(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();

    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;
    
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;

    this->cell_type=E_IN6;
    this->x_index=i;

    
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);
    
    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 


//------------------ Map cells---------------------------------------
class CXsyn_Map: public CellSyn {

 public:
   
 CXsyn_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;
		
    
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
          
    ((RS*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

  void print_c_stren(FILE *connection_record){
        
    double total = 0.0;
    int total_num = 0;
    int i = 0;

    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        total = total + ((AMPA_D2*)syns[i])->strength; 
        total_num = total_num + 1;
      }
    }    
    double average = total / total_num;
    fprintf(connection_record,"%lf %lf %d %d",cur_time,average,m,n);
    
    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        fprintf(connection_record," %lf",((AMPA_D2*)syns[i])->strength);
      }
    }    
    fprintf(connection_record,"\n");
  }
  
};  

class CXasyn_Map: public CellSyn {

 public:

 CXasyn_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((RS*)base_cell)->init();
    
    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

};  

class CXsyn6_Map: public CellSyn {

 public:

 CXsyn6_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((RS*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){


      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

class INsynCore_Map: public CellSyn {
  
 public:  

 INsynCore_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new FS1();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;
    
  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

class INsynMatrix_Map: public CellSyn {

 public:

 INsynMatrix_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new FS1();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

    
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 

class INsyn6_Map: public CellSyn {
  
 public:

 INsyn6_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new FS1();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();
    
    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);  
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
    
    }
  }
}; 



int get_cell_index(int type, int m, int n);
void apply_field(double ***field_effect, double **cx5_local_field, double **cx6_local_field,double **cx5_soma_local_field, double **cx6_soma_local_field,int ii,double time,FILE **field_file);
void allocate_state_save(double ***cx_base_v_SOMA, double ***cxa_base_v_SOMA, double ***cx5_local_field,double ***cx6_local_field,double ***cx5_soma_local_field,double ***cx6_soma_local_field, double ****field_effect);
void spike_fre_calc(int *total_region,double *frequencies);
void boost_activity();
void start_critical();
void end_critical();
void root_critical();
int receive_spike(int m, int n, enum Cell_Type type);
int receive_spike_for_maps(int m, int n, enum Cell_Type type);
double receive_dend(int my_type,int m, int n, enum Cell_Type type);
extern int **max_connect;
extern int **max_connect_gap;

#endif //CellSyn_h
