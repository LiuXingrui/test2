#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <chrono>

#include"BP.h"
#include"modified_BP.h"
#include "Sparse.h"

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


//GlobalRNG_reset (1);
int main(int argc, char **argv){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//For timing.
  GlobalRNG_randomize ();
  double pmax;
  double pmin;
  long long int num_of_cws=0;
  long long int max_num_of_cws=1e18;//when use number_of_decoding_failure, set max_num_of_cws=a very large number
  string H_file;
  string L_file;
  string P_file;
  string data_file;
  int lmax;
  int wt=0;
  int channel=0;
 
  int MMM=0;

  int d=-1;
  int debug=0;
  double pavg;
  double range;
  double alpha;
  int max_row_freeze;
  int num_of_failed_cws=0;
  int max_failed_cws;
  int schedule=0;
  double decode_p,decode_prange,decode_pmin,decode_pmax;
  int other_decoder=0;
  int max_num_del_cols=0;
  //can input pmin and pmax, or the weight of error_vectors and dec_method, which is same_p or diff_p: 
  //diff p decode: pmin=wt/n*0.5, pmax=wt/n*1.5,   same p decode: p=wt/n
  
  cout<<"input files:"<<endl;
  vector<string> input;
  input.push_back("./CBP");
   input.push_back("H_file/Hx_file");
    input.push_back("L_file/Hz_file");
	 input.push_back("P_file/P");
	  input.push_back("max_failures/max_trails(for negative number)");
	   input.push_back("max number of iterations");
	    input.push_back("data_file");
		input.push_back("debug");
		input.push_back("alpha");
		input.push_back("schedule");
		input.push_back("max_num_del_cols");
		input.push_back("max_row_freeze");
		input.push_back("OSD_order");
  
   for (int i=1;i<4;i++)
	  {
		  cout<<i<<": "<<input[i]<<"="<<argv[i]<<endl;
	  }
  if (argc!=13)
  {
	  cout<<" need 12 parameters: ./CBP  H_file  L_file P_file max_failed_cws lmax data_file debug alpha  schedule max_num_del_cols max_row_freeze OSD_order"<<endl;
	  cout<<"you input "<<argc-1<<"parameters"<<"they are"<<endl;	 
	  return 1;
	}
  //get the parameters: 
   
  H_file=argv[1];
  L_file=argv[2];
  P_file=argv[3];
  data_file=argv[6];
 
  istringstream argv4( argv[4] );
  if ( argv4>> max_failed_cws){}
  else
    {
      cout<<"max_failed_cws should be an int"<<endl;
      return 1;
    }

  //cout<<max_failed_cws<<endl;
  if (max_failed_cws<0)
    {
      max_num_of_cws=-max_failed_cws;
      max_failed_cws=100000000;// when use max_num_of_cws, set max_failed_cws to a large number.
    }
 

  istringstream argv5( argv[5] );
  if ( argv5 >> lmax){}
  else
    {
      cout<<"lmax should be an int"<<endl;
      return 1;
    }
 
  istringstream argv7( argv[7] );
  if ( argv7 >> debug){}
  else
    {
      cout<<"debug should be an int"<<endl;
      return 1;
    }

   istringstream argv8( argv[8] );
  if ( argv8 >> alpha){}
  else
    {
      cout<<"alpha should be a double"<<endl;
      return 1;
    }

    istringstream argv9( argv[9] );
    if ( argv9 >> schedule){}
    else
      {
	cout<<"schedule should be an int, 0 for parallel, 1 for serial V, 2 for serial C"<<endl;
	return 1;
      }

    istringstream argv10( argv[10] );
    if ( argv10 >> max_num_del_cols){}
    else
      {
		cout<<"max_num_del_cols should be an int"<<endl;
		return 1;
      }

	  
	   istringstream argv11( argv[11] );
    if ( argv11 >> max_row_freeze){}
    else
      {
		cout<<"max_row_freeze should be an int"<<endl;
		return 1;
      }
	 // if(max_row_freeze>0){cout<<"use freeze rows methods: try at most "<<freeze_rows_suc<<" times"<<endl;}
	  //cout<<"input is ok \n\n"<<endl;
	  
	     if (schedule==0)
	{
	  cout<<"parallel decoding"<<endl;
	}
      else if (schedule==1)
	{
	  cout<<"serial V decodeing"<<endl;
	}
	      else if (schedule==2)
	{
	  cout<<"serial C decodeing"<<endl;
	}
	else if (schedule==3)
			{
	  cout<<"serial V Random decodeing"<<endl;
	}
	else {
		cout<<"schedule should be :0 for parallel, 1 for serial V, 2 for serial C, 3 for random serial V"<<endl;
		return 1;
	}
    int n1,n,r1,n2,r2;
	GF2mat H=read_matrix( n1,r1, H_file);
	n=n1;
  cout<<"the size of H is "<<r1<<"*"<<n1<<endl;
  
int max_rep=1;

int OSD_order=0;
istringstream argv12( argv[12] );
    if ( argv12 >> OSD_order){}
    else
      {
		cout<<"OSD_order should be an int"<<endl;
		return 1;
      }
	  
	  if(OSD_order>0){cout<<"use OSD_2 up to "<<OSD_order<<" bits"<<endl;}
	  else if (OSD_order==0) {cout<<"OSD_1"<<endl;}
	  else if(debug&8){cout<<" OSD_0"<<endl;}
	  else {cout<<" only BP"<<endl;}

  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  bool H_suc=false;

  //read the parity check matrix:

	double error_rate=-1;
    double gate_p=-1;
	double decode_Prange=0;
    vec p(n);
	GF2mat L;
	GF2mat Hz=read_matrix( n2,r2, L_file);
   	 
	istringstream argvP( P_file);
	if (  argvP>> error_rate)
	{
		cout<<"BSC: p="<<error_rate<<endl;
		gate_p=error_rate;
		L=get_gen(Hz);
		r2=L.rows();
		n2=L.cols();
		for (int i=0;i<n;i++) {p(i)=gate_p;}  
		//cout<<p<<endl;
		//istringstream argvL( L_file);
		//if (  argvL>> decode_Prange){}
		//else{cout<<"decode_Prange should be a double"<<endl;return 1;}	
	}
	else
	{
		std::ifstream file(P_file);
		if (!file) {
		std::cerr << "Could not open the file " << P_file << std::endl;
		return 1;
		}
		
		std::vector<double> numbers;
		double num;			
		while (file >> num) {numbers.push_back(num);}
			
		//the last entry of numbers is the depolarzing error rate of each gate
		if (numbers.size()!=n+1) {cout<<"error: numbers.size()!=n+1 numbers.size()="<<numbers.size()<<" n="<<n<<endl;return 0;}
		gate_p=numbers[n];
		for (int i=0;i<n;i++) {p(i)=numbers[i];}  
		cout<<"circuit level decoder, after gates depolarzing rate="<<gate_p<<endl; 
		
		L=Hz;		
		if (n1!=n2)
		{
			cout<<" the two matrices H and L donot match"<<endl;	  
			cout<<"H_file is "<<H_file<<endl;
			cout<<"L_file is "<<L_file<<endl;
			if (n1<100 and n2<100) {cout<<"H"<<H<<endl;cout<<"L is "<<L<<endl;  }
			return 1;
		}  
	}
		 	
// this is an useless matrix, 
  GF2mat G(1,1);
  
  
  Sparse_GF2 sparseH(H);
  Sparse_GF2 sparseL(L);
  Sparse_GF2 Sparse_zero_rvec(r2,1);
  GF2mat zero_mat1(H.rows(),1);
  
  //sparseH.print_cols();
  Sparse_GF2 modified_SH=sparseH;
  Sparse_GF2 modified_SL=sparseL;
  vec modified_p=p;
  
  vector<vector<int>> Positions;
  modified_SH.del_identical_cols(modified_SL,max_num_del_cols, Positions, modified_p) ;
  
  int n_MH=modified_SH.cols();
  
  

  vec LR(n_MH);
  GF2mat modified_H=modified_SH.to_GF2mat();
  GF2mat modified_L=modified_SL.to_GF2mat();
      nodes  checks[r1];//checks for Hx and Z errors
  nodes  errors[n_MH];
  int E1=0; //number of edges in Hx factor graph, but this parameter is not used in this prog
  initialize_checks (modified_H, checks,  E1);
  initialize_errors(modified_H, errors);
  
  	  cout <<"\nMBP_method: delete "<<max_num_del_cols<<" identical columns"<<endl;
	  
	  	  cout <<"new H size=("<<modified_SH.rows()<<","<<modified_SH.cols()<<")"<<endl;
		  
		  
  std::chrono::high_resolution_clock::time_point prep = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> prep_time = (prep - start)/1000;
	 //std::cout << "initialization took " << prep_time.count() << " seconds.\n";
	 
	 if (debug&1)
	 {
				/*
				GF2mat test_e(n,1); 
				test_e.set(20,0,1);
				Sparse_GF2 Stest_e(test_e);
				Stest_e.error_del_cols(Positions);
				cout<<"transformed test_e is"<<endl;
				Stest_e.print_cols();
				*/
	 		cout<<"\noriginal Hx is \n "<<endl;
			//sparseH.print_rows();
			//cout<<"\n "<<endl;
			sparseH.print_cols();
			
						cout<<"identical cols are (after deleting, index of cols are changing) \n"<<endl;
			    for (int i = 0; i < Positions.size(); i++) {
					
					for (int j = 0; j < Positions[i].size(); j++) {
						std::cout << Positions[i][j] << ' ';
					}
				std::cout << std::endl;
				}
			cout<<"\noriginal Gx (zero space of Hz) is \n "<<endl;
			//sparseL.print_rows();
			//cout<<"\n "<<endl;
			sparseL.print_cols();
			
			cout<<"\ntransformed Hx is \n "<<endl;
				//modified_SH.print_rows();
				//cout<<"\n "<<endl;
				modified_SH.print_cols();
				
			cout<<"\ntransformed Lx is \n "<<endl;
				//modified_SL.print_rows();		
				//cout<<"\n "<<endl;
				modified_SL.print_cols();		
	 }
	 
  while (  num_of_failed_cws<max_failed_cws&&num_of_cws<max_num_of_cws)
    {
		num_of_cws++;
		H_suc=false;
	  
		int wt_real_e=0;
		GF2mat real_e(n,1); 
		wt_real_e=error_channel(real_e, p);

		if (wt_real_e==0)
		{
			num_no_error++;
			H_suc=true;			
		}
		else
		{	GF2mat real_Me=real_e; 
			Sparse_GF2  Sparse_real_Me(real_Me);
			Sparse_real_Me.error_del_cols(Positions);
			
			H_suc=MBP_del_cols(modified_SH,modified_SL,sparseH,sparseL,Sparse_real_Me,real_e,Positions,checks,errors,modified_p,num_iter, 
 lmax,wt,  debug, schedule,OSD_order, Sparse_zero_rvec,  zero_mat1);
		}
   
		if (H_suc==true)
		{       	 
			num_of_suc_dec++;
		}
		else
		{
			num_of_failed_cws++;
		}  
    }

  //if (MMM==214748364){MMM=0;}
  //  cout<<"lambda="<<lambda<<endl;
   
    cout<<"\n"<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws "<<endl;
	cout<<"logical error rate="<<1.0-1.0*num_of_suc_dec/num_of_cws<<endl;
   cout<<"average iterations="<<num_iter/num_of_suc_dec<<"\n"<<endl;
   cout<<"num of trivial error is "<<num_no_error<<endl;
    cout<<" BP_suc="<< BP_suc<<endl;
	cout<<"OSD_suc ="<<OSD_suc<<endl;
   cout<<"OSD_0_suc="<<OSD_0_suc<<endl;
   cout<<"OSD_1_suc="<<OSD_1_suc<<endl;
   cout<<"OSD_2_suc="<<OSD_2_suc<<endl;
   cout<<"OSD_fail="<<OSD_fail<<endl;
   cout<<"syn_fail="<<syn_fail<<endl;
   cout<<"max_fail="<<max_fail<<endl;

   // cout<<"num of zero errors is about "<<pow(p,n)*num_of_cws<<endl;
      
   ofstream myfile;
   myfile.open (data_file,ios::app);
  int failed_trials=num_of_cws-num_of_suc_dec;
   int before_OSD_fail=num_of_cws-BP_suc-num_no_error;
   int OSD_0_fail=num_of_cws-BP_suc-num_no_error-OSD_0_suc;
    int OSD_1_fail=num_of_cws-BP_suc-num_no_error-OSD_0_suc-OSD_1_suc;
     int OSD_2_fail=num_of_cws-BP_suc-num_no_error-OSD_0_suc-OSD_1_suc-OSD_2_suc;
     //  int OSD_higher_order_fail=num_of_cws-BP_suc-num_no_error-OSD_0_suc-OSD_1_suc-OSD_2_suc-OSD_higher_order_suc;
   int OSD_trails=(OSD_suc+OSD_fail);
  
   myfile<<n<<" "<<gate_p<<" "<<num_of_cws<<" "<<failed_trials<< " "<<before_OSD_fail<<" "<<OSD_0_fail<<" "<<OSD_1_fail<<" "<<OSD_2_fail<< " "<<OSD_trails<<endl;
  
     //myfile << n<<" "<<d<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<gate_p<<" "<<1.0*num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<" "<<num_of_cws<<"  "<<syn_fail<<" "<<max_fail<<" "<<1.0*syn_fail/num_of_cws<<" "<<1.0*max_fail/num_of_cws<<" "<<alpha<<" "<<OSD_suc<<" "<<MMM<<"  "<<LLR_suc<<"  "<<LLR_fail<<endl;
    
   myfile.close();
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "\n Run-time " << time_span.count() << " seconds.\n";
      cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  return 0;

}
  
