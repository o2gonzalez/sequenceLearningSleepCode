//Created by M.Bazhenov, 1997-2009
//All rights reserved
//Most recently edited by Peter Lonjers(just search the net if you need to find me)
//I added all the MPI stuff but don't know much about the bio
#include <omp.h>
#include "CellSyn.h"    //header for all classes describing currents and basic cell types
#include "io.h"
#include "network.h"

int num_mp_threads; //number of openMP threads

int print_c_sten; // if 1 prints all the connection strengths for cx every once in awile note this only works with openmp
int fre_print_cs; // how often to print connection strengths(every how many milliseconds)
int total_rec; // used for debugging number of non normal messages received in a time step
int total_rec_norm; // used for debugging number of normal(basic spiking messages) messages received in a time step

string output_location; // keeps track of what folder to put outputs in

Homeostasis homeo;

//block size for openmp
int block_size = 13;

//boolean for if we are doing the very first step and send functions
int first_round = 1;
//main track of milliseconds of simulation so far. ii tracks steps of simulation locally
double t = 0;

//how long we run the simulation
double tmax;
//time at which we print information from all cells
double t3D;
//smaller time interval to print information from some cells
double ttime;
//how long main process runs for set in input
double run_time = 0;

double input_strength,input_start_neuron,input_end_neuron;

// number of types of cells
int Num_Types;
//gives the dimensions of the arrays of cells of each type
Pair *cell_sizes;
//total number of cells
int num_cells;
//total number of cells of each type
int *cell_numbers;
//cells info is an array of all cells containing the color of each cell and the numbers of different types of connections it has.
Cell_Info ****cells_info = NULL;
//array of all the pointers to cells belonging to this process
CellSyn **cells =NULL;

int *mapns, *hhns;

int mapi,hhi;

int one_state=0;
int test_duration=50000;
int seq1_duration=250000;
int seq2_duration=250000;
int seq3_duration=300000;
int n3_duration=1000000;

//--- test baselines ---//
// test seq1 baseline
int test1_start=12001;
int test1_end=test1_start+test_duration-1;

// test seq2 baseline;
int test2_start = test1_end+1;
int test2_end = test2_start+test_duration-1;

// test seq3 baseline
int test3_start = test2_end+1;
int test3_end = test3_start+test_duration-1;

//--- traing sequence 1 ---//
// sequence 1 training
int seq1_trainingStart=test3_end+1;
int seq1_trainingEnd=seq1_trainingStart+seq1_duration-1;

//--- tests after sequence 1 training ---//
// test seq1
int test4_start = seq1_trainingEnd+1;
int test4_end = test4_start+test_duration-1;

// test seq2
int test5_start = test4_end+1;
int test5_end = test5_start+test_duration-1;

// test seq3
int test6_start = test5_end+1;
int test6_end = test6_start+test_duration-1;

//--- training sequence 2 ---//
// sequence 2 training
int seq2_trainingStart=test6_end+1;
int seq2_trainingEnd=seq2_trainingStart+seq2_duration-1;

//--- tests after sequence 2 training ---//
// test seq1
int test7_start = seq2_trainingEnd+1;
int test7_end = test7_start+test_duration-1;

// test seq2
int test8_start = test7_end+1;
int test8_end = test8_start+test_duration-1;

// test seq3
int test9_start = test8_end+1;
int test9_end = test9_start+test_duration-1;

//--- training sequence 3 ---//
// sequence 3 training
int seq3_trainingStart = test9_end+1;
int seq3_trainingEnd = seq3_trainingStart+seq3_duration-1;

//--- tests after sequence 3 training ---//
// test seq1
int test10_start = seq3_trainingEnd+1;
int test10_end = test10_start+test_duration-1;

// test seq2
int test11_start = test10_end+1;
int test11_end = test11_start+test_duration-1;

// test seq3
int test12_start = test11_end+1;
int test12_end = test12_start+test_duration-1;

//--- sleep ---//
int awake_end = test12_end+1;
int stage3_end=awake_end+n3_duration;

//--- test after sleep ---//
// test seq1
int test13_start = stage3_end+1;
int test13_end = test13_start+test_duration-1;

// test seq2
int test14_start = test13_end+1;
int test14_end = test14_start+test_duration-1;

// test seq3
int test15_start = test14_end+1;
int test15_end = test15_start+test_duration-1;

/// CX neuron -gKl
double fac_gkl; 
double fac_gkl_TC; 
double fac_gkl_RE; 

double fac_gkv_cx; 
double fac_gkca_cx;
double fac_gkm_cx;

double fac_gh_TC;

double fac_AMPA_D2;
double fac_AMPA_TC;
double fac_GABA_D2;
double fac_GABA_TC;


//run for openMP
//all == whether both RK and MAP should be called at the same time
//the order is important in order to ensure signaling and receiveing spikes between different steps

void runMP(int all){
  int iter,i;
#pragma omp parallel  private(i,iter)  num_threads(num_mp_threads)
  {
    // Run step for maps
    if (all){
#pragma omp for schedule(dynamic, block_size)
      for(i=0; i<mapi; i++){
        cells[mapns[i]]->step();
      }
    }

    for(iter=0; iter<4; iter++){
#pragma omp for schedule(dynamic, block_size)
      for(i=0; i<hhi; i++){
        cells[hhns[i]]->step();
      }
    }

#pragma omp for schedule(dynamic, block_size)
    for(i=0; i<hhi; i++){
      cells[hhns[i]]->reset_y();
    }

    // Signal and/or reset spikes here
    if (all) {
#pragma omp for schedule(dynamic, block_size)
      for(i=0; i<mapi; i++) {
        //Reset flip for maps which should be consumed in this time step (used for HH->Map connection)
        cells[mapns[i]]->reset_flip_maps();
      }

#pragma omp for schedule(dynamic, block_size)
      for(i=0; i<hhi; i++) {
        cells[hhns[i]]->reset_flip_maps();
      }
    }

#pragma omp for schedule(dynamic, block_size)
    for(i=0; i<mapi; i++) {
      cells[mapns[i]]->signal_spike();
    }

#pragma omp for schedule(dynamic, block_size)
    for(i=0; i<hhi; i++) {
      cells[hhns[i]]->signal_spike();
    }


  }
}

//gets index for a particular cell in a group(cells on the same process)
int get_cell_index(int type, int m, int n){
  return cells_info[type][m][n]->cell_index;
}

//all arguments are inputs
//this functions decides what cells a process should simulate and starts them
void initialize_cells(CellSyn** &cells){
  int total=0;
  for(int x=0; x<Num_Types; x++)
    for(int y=0; y<cell_sizes[x].x; y++)
      for(int z=0; z<cell_sizes[x].y; z++)
        cells_info[x][y][z]->cell_index = total++;

  int type = 0;
  int m = 0;
  int n = 0;
  int total_set = 0;

  cells = new CellSyn*[num_cells];
  for(m = 0; m<num_cells; m++){
    cells[m]= NULL;
  }

  for(type = 0; type<Num_Types; type++){
    printf("Initializing cell %d for %d, %d\n",type, cell_sizes[type].x, cell_sizes[type].y);
    for(m=0; m<cell_sizes[type].x; m++){
      for(n=0; n<cell_sizes[type].y; n++){
        cells[get_cell_index(type,m,n)]=CellSyn::initialize(type,m,n,cells_info,output_location);
        total_set++;
      }
    }
  }

  //sanity checks
  if(total_set != num_cells){
    printf("process not given enough cells. This should not happen\n");
    printf("number given: %d\n",total_set);
    printf("number needed: %d\n",num_cells);
    exit(1);
  }

  for(m=0;m<num_cells;m++){
    if(cells[m]==NULL){
      printf("graph coloring probably failed\n");
      exit(1);
    }
  }
  //all checks passed
  return;
}


// Scale connection strength by number of inputs
void scale_synapses(CellSyn** &cells, int num_cells){
  printf("Scaling Synapses ... \n");
  double factor=1;
  double rr=50;
  for(int i=0; i<num_cells; i++) {
    if(i>=200&&i<200+rr)
    {
      factor=1;
    }
    else if(i>=700-rr&&i<700)
      factor=1;
    else factor=1;
    cout<<i<<":"<<cells[i]->type<<" factor:"<<factor<<endl;
    // -------------- Scaling Method ----------------------
    for(int k=0; k < cells[i]->num_syns; k++) {
      if (cells[i]->syns[k]->type==E_GAP) {
        if ((max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type])==0) {
          printf("Error during scaling GAP synapse for %d from %d =0 \n",cells[i]->type, cells[i]->syns[k]->from_type);
          exit(1);
        }
        cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type];
        cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/max_connect_gap[cells[i]->type][cells[i]->syns[k]->from_type]*factor;
      }
      else{
        if ((max_connect[cells[i]->type][cells[i]->syns[k]->from_type])==0) {
          printf("Error during scaling synapse for %d from %d =0 \n",cells[i]->type, cells[i]->syns[k]->from_type);
          exit(1);
        }
        cells[i]->syns[k]->strength= cells[i]->syns[k]->strength/max_connect[cells[i]->type][cells[i]->syns[k]->from_type];
        cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s/max_connect[cells[i]->type][cells[i]->syns[k]->from_type]*factor;
      }
    }
  }
  return;
}

void send_IStim1(int i, int j,double cx_stim, double cxa_stim){
  int group_index_CXa = get_cell_index(E_CXa,i,j);
  int group_index_CX = get_cell_index(E_CX,i,j);
  ((CX*)cells[group_index_CX]->base_cell)->I_Stim1 = cx_stim;
  ((CX*)cells[group_index_CXa]->base_cell)->I_Stim1 = cxa_stim;
}

//receive for getting things to the root printing
double print_receive(int m , int n, enum Cell_Type type){
  double message = 0;
  int group_index = get_cell_index(type,m,n);
  message = cells[group_index]->base_cell->get_v_soma();
  return message;
}

//receive for getting things to the root printing
void print_receive_gsynapse_index(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp){
  int group_index = get_cell_index(type,m,n);
  // Run through all synapses
  for(int i =0; i < cells[group_index]->num_syns; i++){
    if ((cells[group_index]->syns[i]->type == stype) && ((cells[group_index]->syns[i]->from_type == from_type))){   
      int from=cells[cells[group_index]->syns[i]->from_cell]->m;
      fprintf(fp,"%d %d ", from ,m);
    }
  }
}

//receive for getting things to the root printing
void print_receive_gsynapse(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp){
  int group_index = get_cell_index(type,m,n);
  // Run through all synapses
  for(int i =0; i < cells[group_index]->num_syns; i++){
    if (cells[group_index]->syns[i]->type == stype){   
      fprintf(fp,"%d %lf ", m, ((AMPA_D3*)(cells[group_index]->syns[i]))->g_AMPA0); //D3
    }
  }
}

//receive minis for printing
void print_receive_minis(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp){
  int group_index = get_cell_index(type,m,n);
  // Run through all synapses
  for(int i =0; i < cells[group_index]->num_syns; i++){
    if (cells[group_index]->syns[i]->type == stype){   
      fprintf(fp,"%d %lf ", m, ((AMPA_D3*)(cells[group_index]->syns[i]))->mini_s);
    }
  }
}

//receive for getting how many times a cell has spiked during the fre_window
int fre_receive(int m , int n, enum Cell_Type type){
  int message = 0;
  int group_index = get_cell_index(type,m,n);
  message = cells[group_index]->num_spikes;
  return message;
}

//receives dedritic voltage from connecting cell and returns it
//all args inputs
//m n and type are the index of the sending cell, my_type is the index of the receiving cell
double receive_dend(int my_type,int m, int n, enum Cell_Type type){
  int group_index = get_cell_index(type,m,n);
  return cells[group_index]->base_cell->get_v_dend();
}

//receives voltage from connecting cell and returns it
//all args inputs
//m n and type are the index of the sending cell, my_type is the index of the receiving cell
int receive_spike(int m, int n, enum Cell_Type type){
  int group_index = get_cell_index(type,m,n);
  return cells[group_index]->get_flip();
}

int receive_spike_for_maps(int m, int n, enum Cell_Type type){
  int group_index = get_cell_index(type,m,n);
  return cells[group_index]->get_flip_for_maps();
}

extern FILE *f28; 

void load_input_data(int argc, char *argv[]){ 
  //checks inputs
  if (argc < 4) {
    puts("Bad command: Should be");
    puts("input_file_name output_directory connections_file_name");
    exit(1);
  }
  output_location = argv[2];
  output_location = output_location + "/";

  //open up the file of connections
  FILE *connections_file;
  if (!(connections_file=fopen(argv[3],"r"))) {
    printf("%s file for connections does not exist or something\n",argv[3]);
    exit(1);
  }

  set_cell_sizes(connections_file,Num_Types,cell_sizes,cell_numbers); //scans connections file for how many cells of each type we have
  create_connections_from_file(&cells_info,cell_sizes,connections_file,Num_Types);

  //finds total number of cells
  num_cells = 0;
  for(int iter = 0; iter<Num_Types; iter++){
    num_cells = num_cells+cell_numbers[iter];
  }
  printf("number cells: %d\n",num_cells);

}

double correct1_func(double value, double facp){
  if (value<1){
    return value*(1-facp);
  } else{
    return value*(1+facp);    
  }
}

double linear_change(double start_value, double tdiff, double end_value){
  return start_value  + (tdiff * (end_value - start_value));
}

FILE *fmanipulations;

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++
int main(int argc,char **argv){

  LocalFieldPotential LFP;

  load_input_data(argc,argv);
  load_input_params(argc,argv,tmax,t3D,ttime,num_mp_threads,print_c_sten,fre_print_cs,LFP.local_field_effect,LFP.lfp_scale,LFP.num_field_layers,homeo.boost,homeo.amp_boost,homeo.con_boost,homeo.fre_boost,homeo.target_f,homeo.fre_window,homeo.num_regions,input_strength,input_start_neuron,input_end_neuron, block_size);

  // seed random number generator with current second of execution
  srand(2000);

  homeo.allocate();

  initialize_cells(cells);

  scale_synapses(cells, num_cells);

  LFP.init();
  open_files(output_location,LFP.field_file,LFP.num_field_layers);
  LFP.allocate_state_save(cell_sizes);

  //------Main Loop--------
  printf("timer started");
  printf("\n starting main loop: time= %lf: end_time= %lf\n", t,tmax);
  int print_count = 0;
  int ii = 0;
  int i=0;
  hhns=new int[num_cells];
  for(i=0; i<num_cells; i++){
    if (cells[i]->ismap){
      mapns[mapi]=i;
      mapi=mapi+1;
    }
    else{
      hhns[hhi]=i;
      hhi=hhi+1;
    }
  }

  print_used_mem();
  time_t start_time = time(NULL);
  printf("root set up tau:%lf time:%lf tmax:%lf\n",TAU,t,tmax);

  double s3_scale=2.0;

  /// gKl 
  double gkl_awake_fix     = 0.19; 
  double gkl_s3            = gkl_awake_fix*s3_scale;

  double gkl_TC_awake_fix  = 0.8;
  double gkl_TC_s3         = gkl_TC_awake_fix*s3_scale;

  double gkl_RE_awake_fix  = 0.9;
  double gkl_RE_s3         = gkl_RE_awake_fix*((2-s3_scale/2)-0.5);

  double awake_AMPAd2_fix  =0.19; //wake_ach_fac;
  double s3_AMPAd2         =awake_AMPAd2_fix*s3_scale*1.28;

  double gh_TC_awake       =-8.0;
  double gh_TC_s3          =-1.0;

  double awake_GABAd2_fix  =0.22;
  double s3_GABAd2         =awake_GABAd2_fix*s3_scale;

  double awake_GABA_TC_fix =0.6; //awake_gaba_fac;
  double s3_GABA_TC        =awake_GABA_TC_fix*s3_scale; 

  fmanipulations = fopen("manipulation.txt","w");

  double awake_AMPA_TC     =0.5; //awake_ach_fac;
  double s3_AMPA_TC        =0.5;  

  double gk_cx_slow_awake  = 1.0; 
  double gk_cx_slow_s3     = 1.0; 

  double gk_cx_spike_awake = 1.0;
  double gk_cx_spike_s3    = 1.0;


  print_connectivity(cell_sizes);

  while( t < tmax){
    total_rec = 0;

    if ((t<=awake_end)){

      // Fix all values for awake state
      fac_AMPA_D2 = awake_AMPAd2_fix*0.7;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix*0.5;
      fac_gkl     = gkl_awake_fix*0.7;
      fac_gh_TC   = gh_TC_awake*3;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;

    }else if ((t>awake_end&&t<=stage3_end)){ 
      // Fix all values for S3
      fac_AMPA_D2 = s3_AMPAd2;
      fac_AMPA_TC = s3_AMPA_TC;
      fac_GABA_D2 = s3_GABAd2;
      fac_GABA_TC = s3_GABA_TC;
      fac_gkl_RE  = gkl_RE_s3;
      fac_gkl_TC  = gkl_TC_s3;
      fac_gkl     = gkl_s3;
      fac_gh_TC   = gh_TC_s3;
      fac_gkca_cx = gk_cx_slow_s3;
      fac_gkm_cx  = gk_cx_slow_s3;
      fac_gkv_cx  = gk_cx_spike_s3;
    }
    else { //awake after sleep
      fac_AMPA_D2 = awake_AMPAd2_fix*0.7;
      fac_AMPA_TC = awake_AMPA_TC;
      fac_GABA_D2 = awake_GABAd2_fix;
      fac_GABA_TC = awake_GABA_TC_fix;
      fac_gkl_RE  = gkl_RE_awake_fix;
      fac_gkl_TC  = gkl_TC_awake_fix*0.5;
      fac_gkl     = gkl_awake_fix*0.7;
      fac_gh_TC   = gh_TC_awake*3;
      fac_gkca_cx = gk_cx_slow_awake;
      fac_gkm_cx  = gk_cx_slow_awake;
      fac_gkv_cx  = gk_cx_spike_awake;
    }

    if(t>print_count){
      printf("milliseconds of sim = %lf \n", t);
      print_count = print_count + 200;
    }

    ii = ii + 1; //number iterations

    if (((ii)/(TAUr))*(TAUr) == (ii))
      runMP(1);
    else
      runMP(0);

    //applies field effects
    if(LFP.local_field_effect == 1){
      LFP.apply_field(ii,t,ttime,t,cell_sizes);
    }

    //prints all 1D cell Voltages or all x cells for one y in 2d voltages occasionally
    if((t > ttime) && ((ii/(50))*(50) == ii)){ //multiples of 50 print
      print_freq(LFP.cx_base_v_SOMA,LFP.cxa_base_v_SOMA,cell_sizes,t);
    }

    //prints all cell voltages occationally
    if((t > t3D) && (((ii)/(50000))*(50000) == (ii))){ //multiples of 500 print
      print_occ(cell_sizes);
      fprintf(fmanipulations, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", t, fac_AMPA_D2, fac_AMPA_TC, fac_GABA_D2, fac_GABA_TC, fac_gkl_RE, fac_gkl_TC, fac_gkl, fac_gh_TC, fac_gkca_cx, fac_gkm_cx, fac_gkv_cx);
    }

    //prints strength of connections
    if((ii >= fre_print_cs) && (ii % fre_print_cs) == 0 && print_c_sten == 1){
      int i =0;
      //TODO we are able to enumerate E_CX directly
      for(i=0; i<num_cells; i++){
        if(cells[i]->type == E_CX){
          ((CXsyn*)cells[i])->print_c_stren(f28);
        }
      }
    }

    //implements homeostatic mechanisms
    if((ii >= homeo.fre_window) && (ii % homeo.fre_window) == 0){
      homeo.spike_fre_calc(cell_sizes,cell_numbers); //figures out frequency of spiking
      int i = 0;
      for(i=0; i<num_cells; i++){
        if(cells[i]->type == E_CX || cells[i]->type == E_CXa || cells[i]->type == E_CX6){
          cells[i]->num_spikes=0;
        }
      }
      homeo.boost_activity(cell_sizes,cells,num_cells); //causes homeostatic changes
    }

    t = t + TAU; //increase time
  } //end while loop

  printf("Computation time: %d s\n",(int)difftime(time(NULL),start_time));
  close_files(LFP.field_file,LFP.num_field_layers);

  fclose(fmanipulations);
  return 0;
}
