#include "io.h"
#include <sstream>
#include "params.h"
#include <string.h>

double print_receive(int m , int n, enum Cell_Type type); 
void print_receive_gsynapse(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp);
void print_receive_minis(int m , int n, enum Cell_Type type, enum Syn_Type stype, FILE *fp);
void print_receive_gsynapse_index(int m , int n, enum Cell_Type from_type, enum Cell_Type type, enum Syn_Type stype, FILE *fp);

void load_input_params(
  int argc,char *argv[],

  double& tmax,
  double&  t3D,
  double&  ttime,
  int&  num_mp_threads,
  int&  print_c_sten,
  int&  fre_print_cs,
  int&  LFP_local_field_effect,
  double&  LFP_lfp_scale,
  int&  LFP_num_field_layers,
  int&  homeo_boost,
  double&  homeo_amp_boost,
  double&  homeo_con_boost,
  double&  homeo_fre_boost,
  double&  homeo_target_f,
  int&  homeo_fre_window,
  int&  homeo_num_regions,
  double& input_strength,
  double& input_start_neuron,
  double& input_end_neuron,
	int& block_size

 ){
  
  add_double_param(  tmax  );
  add_double_param(  t3D   );
  add_double_param(  ttime );
  add_int_param(  num_mp_threads );
  add_int_param(  print_c_sten   );
  add_int_param(  fre_print_cs   );
  add_int_param(  LFP_local_field_effect );
  add_double_param(  LFP_lfp_scale );
  add_int_param(  LFP_num_field_layers );
  add_int_param(  homeo_boost  );
  add_double_param(  homeo_amp_boost  );
  add_double_param(  homeo_con_boost  );
  add_double_param(  homeo_fre_boost  );
  add_double_param(  homeo_target_f   );
  add_int_param(  homeo_fre_window    );
  add_int_param(  homeo_num_regions   );

  add_double_param(input_strength);
  add_double_param(input_start_neuron);
  add_double_param(input_end_neuron);
	add_int_param(block_size);

  assert(load_parameters(argv[1]));
  assert(cmdline_parameters(argc,argv));
  print_parameters();

} //load_input_params

//giant list of all the output files 
FILE *flocal, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16, *f17, *f18, *f19,*f20,*f21,*f22,*f23,*f24,*f25,*f26,*f27,*f28,*f30,*f32; 
FILE *fconn;

void print_connectivity(Pair *cell_sizes){

  int i =0;
  int j =0;
  
  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      print_receive_gsynapse_index(i, j, E_CX, E_CX, E_AMPA_D3, fconn);
    }    
  }
  fprintf(fconn,"\n");
}

void print_occ(Pair *cell_sizes){

  int i =0;
  int j =0;

  for(i = 0; i < cell_sizes[E_RE].x; ++i){
    for(j = 0; j < cell_sizes[E_RE].y; ++j){
      fprintf(f6,"%lf ", print_receive(i,j,E_RE));
    }
    fprintf(f6,"\n");
  }

  for(i = 0; i < cell_sizes[E_TC].x; ++i){
    for(j = 0; j < cell_sizes[E_TC].y; ++j){
      fprintf(f8,"%lf ", print_receive(i,j,E_TC));
    }
    fprintf(f8,"\n");
  }

  for(i = 0; i < cell_sizes[E_TCa].x; ++i){
    for(j = 0; j < cell_sizes[E_TC].y; ++j){
      fprintf(f16,"%lf ", print_receive(i,j,E_TCa));
    }
    fprintf(f16,"\n");
  }

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      fprintf(f10,"%lf ", print_receive(i,j,E_CX));
    }
    fprintf(f10,"\n");
  }
  
  for(i = 0; i < cell_sizes[E_CXa].x; ++i){
    for(j = 0; j < cell_sizes[E_CXa].y; ++j){
      fprintf(f14,"%lf ", print_receive(i,j,E_CXa));
    }
    fprintf(f14,"\n");
  }

  for(i = 0; i < cell_sizes[E_IN].x; ++i){
    for(j = 0; j < cell_sizes[E_IN].y; ++j){
      fprintf(f12,"%lf ", print_receive(i,j,E_IN));
    }
    fprintf(f12,"\n");
  }

  // TB IN layers
  for(i = 0; i < cell_sizes[E_INa].x; ++i){
    for(j = 0; j < cell_sizes[E_INa].y; ++j){
      fprintf(f18,"%lf ", print_receive(i,j,E_INa));
    }
    fprintf(f18,"\n");
  }
  

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
     print_receive_gsynapse(i, j, E_CX, E_AMPA_D3, f20);
    }    
  }
  fprintf(f20,"\n");
  

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
     print_receive_minis(i, j, E_CX, E_AMPA_D3, f30);
    }    
  }
  fprintf(f30,"\n");  

  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    for(j = 0; j < cell_sizes[E_CX6].y; ++j){
      fprintf(f24,"%lf ", print_receive(i,j,E_CX6));
    }
    fprintf(f24,"\n");
  }

  for(i = 0; i < cell_sizes[E_IN6].x; ++i){
    for(j = 0; j < cell_sizes[E_IN6].y; ++j){
      fprintf(f26,"%lf ", print_receive(i,j,E_IN6));
    }
    fprintf(f26,"\n");
  }
}

// This works best for sparse firing
void print_spike_time(double time, int type, int m, int n){
  fprintf(f32, "%lf %d %d %d \n",time,type,m,n);
 }


void print_freq(double **cx_base_v_SOMA, double **cxa_base_v_SOMA, Pair *cell_sizes, double const t){

  int i = 0;
  
  fprintf(f7,"%lf ", t);
  for(i = 0; i < cell_sizes[E_RE].x; ++i){
    fprintf(f7,"%lf ",print_receive(i,cell_sizes[E_RE].y/2,E_RE));
  }
  fprintf(f7,"\n");
  
  fprintf(f9,"%lf ", t);
  for(i = 0; i < cell_sizes[E_TC].x; ++i){
    fprintf(f9,"%lf ", print_receive(i,cell_sizes[E_TC].y/2,E_TC));
  }
  fprintf(f9,"\n");

  fprintf(f17,"%lf ", t);
  for(i = 0; i < cell_sizes[E_TCa].x; ++i){
    fprintf(f17,"%lf ", print_receive(i,cell_sizes[E_TCa].y/2,E_TCa));
  }
  fprintf(f17,"\n");
  
  //print time
  fprintf(f11,"%lf ", t);
  
  double total = 0.0;
  double temp_cx=0.0;

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    temp_cx = print_receive(i,cell_sizes[E_CX].y/2,E_CX);
    
    //print individual cell data
    fprintf(f11,"%lf ", temp_cx);
    total = total + temp_cx;
  }

  //print average 
  fprintf(f11,"%lf ", total/cell_sizes[E_CX].x);

  fprintf(f11,"\n");
  
  //print time
  fprintf(f15,"%lf ", t);
  
  double temp_cxa  = 0.0;
  double total_cxa = 0.0;

  
  for(i = 0; i < cell_sizes[E_CXa].x; ++i){
    temp_cxa = print_receive(i,cell_sizes[E_CXa].y/2,E_CXa);
    total_cxa = total_cxa + temp_cxa;

    //print individual cell data
    fprintf(f15,"%lf ", temp_cxa);
  }
  //print average 
  fprintf(f15,"%lf ", total_cxa/cell_sizes[E_CXa].x);
  fprintf(f15,"\n");

  fprintf(f13,"%lf ", t);
  for(i = 0; i < cell_sizes[E_IN].x; ++i){
    fprintf(f13,"%lf ", print_receive(i,cell_sizes[E_IN].y/2,E_IN));
  }
  fprintf(f13,"\n");

  fprintf(f19,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_INa].x; ++i){
    fprintf(f19,"%lf ", print_receive(i,cell_sizes[E_INa].y/2,E_INa));  // TB IN layers
  }
  fprintf(f19,"\n"); // TB IN layers

  fprintf(f25,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    fprintf(f25,"%lf ", print_receive(i,cell_sizes[E_CX6].y/2,E_CX6));  // TB IN layers
  }
  fprintf(f25,"\n"); // TB IN layers

  
  fprintf(f27,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_IN6].x; ++i){
    fprintf(f27,"%lf ", print_receive(i,cell_sizes[E_IN6].y/2,E_IN6));  // TB IN layers
  }
  fprintf(f27,"\n"); // TB IN layers
}

void open_files(string output_location,FILE **field_file, int num_field_layers){
  printf("open files\n");
  int i = 0;
  for(i=0; i<num_field_layers; i++){
    stringstream ss;
    ss << i;
    field_file[i] = fopen((output_location+"field_file_" + ss.str()  ).c_str(), "w");
  }
  
  f2 = fopen((output_location+"dat").c_str(), "w");
  f6 = fopen((output_location+"graf_re").c_str(), "w");
  if (!(f7=fopen((output_location+"time_re").c_str(), "w"))){
      printf("probably out put folder doesn't exist\n");
      exit(1); 
  }
  //f7 = fopen((output_location+"time_re").c_str(), "w");
  f8 = fopen((output_location+"graf_tc").c_str(), "w");
  f9 = fopen((output_location+"time_tc").c_str(), "w");
  f10 = fopen((output_location+"graf_cx").c_str(), "w");
  f11 = fopen((output_location+"time_cx").c_str(), "w");
  f12 = fopen((output_location+"graf_in").c_str(), "w");
  f13 = fopen((output_location+"time_in").c_str(), "w");
  f14 = fopen((output_location+"graf_cxa").c_str(), "w");
  f15 = fopen((output_location+"time_cxa").c_str(), "w");
  f16 = fopen((output_location+"graf_tca").c_str(), "w");
  f17 = fopen((output_location+"time_tca").c_str(), "w");
  f18 = fopen((output_location+"graf_ina").c_str(), "w"); // TB IN layers
  f19 = fopen((output_location+"time_ina").c_str(), "w"); // TB IN layers 
  f20 = fopen((output_location+"time_G_AMPA0_CX_CX").c_str(), "w");
  f21 = fopen((output_location+"time_G_AMPA0_CXa_CX").c_str(), "w");
  f22 = fopen((output_location+"time_G_AMPA0_CX_CXa").c_str(), "w"); // TB IN layers
  f23 = fopen((output_location+"time_G_AMPA0_CXa_CXa").c_str(), "w"); // TB IN layers

  f24 = fopen((output_location+"graf_cx6").c_str(), "w"); // TB IN layers
  f25 = fopen((output_location+"time_cx6").c_str(), "w"); // TB IN layers
  f26 = fopen((output_location+"graf_in6").c_str(), "w"); // TB IN layers
  f27 = fopen((output_location+"time_in6").c_str(), "w"); // TB IN layers

  f28 = fopen((output_location+"cx_cx_g_ampa0").c_str(), "w");

  f30 =fopen((output_location+"time_minis_AMPA_CX_CX").c_str(), "w"); 
  f32 = fopen((output_location+"spike_time").c_str(),"w");
  fconn = fopen((output_location+"conn_index_cx_cx_G").c_str(), "w"); // 

  printf("files open for write\n");
}


void close_files(FILE **field_file, int num_field_layers){

  int i = 0;
  for(i=0; i<num_field_layers; i++){
    fclose(field_file[i]);
  }

  if(f2!=NULL){
    fclose(f2);
  }
  if(f6!=NULL){
    fclose(f6);
  }
  if(f7!=NULL){
    fclose(f7);
  }
  if(f8!=NULL){
    fclose(f8);
  }
  if(f9!=NULL){
    fclose(f9);
  }
  if(f10!=NULL){
    fclose(f10);
  }
  if(f11!=NULL){
    fclose(f11);
  }
  if(f12!=NULL){
    fclose(f12);
  }
  if(f13!=NULL){
    fclose(f13);
  }
  if(f14!=NULL){
    fclose(f14);
  }
  if(f15!=NULL){
    fclose(f15);
  }
  if(f16!=NULL){
    fclose(f16);
  }
  if(f17!=NULL){
    fclose(f17);
  }
  if(f18!=NULL){
    fclose(f18);
  }  // TB IN layers
  if(f19!=NULL){
    fclose(f19);
  }
  if(f20!=NULL){
    fclose(f20);
  }
  if(f21!=NULL){
    fclose(f21);
  }
  if(f22!=NULL){
    fclose(f22);
  }
  if(f23!=NULL){
    fclose(f23);
  }
  if(f24!=NULL){
    fclose(f24);
  }
  if(f25!=NULL){
    fclose(f25);
  }
  if(f26!=NULL){
    fclose(f26);
  }
  if(f27!=NULL){
    fclose(f27);
  }
  if(f28!=NULL){
    fclose(f28);
  }
  if(f30!=NULL){
    fclose(f30);
  }
  if(fconn!=NULL){
    fclose(fconn);
  }
  if(f32!=NULL){
    fclose(f32);
  }
}

//returns MB
void print_used_mem(){
  FILE* proc = fopen("/proc/self/status", "r");
  string res;
  char line[128];

  while (fgets(line, 128, proc)){
      if (!strncmp(line, "VmRSS:", 6)){
	  printf("Resident RAM: %s",line);
          break;
      }
  }
  fclose(proc);
} 
