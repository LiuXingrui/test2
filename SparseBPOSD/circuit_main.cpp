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
  double decode_p,decode_prange,decode_pmin,decode_pmax;
  int other_decoder=0;
  //can input pmin and pmax, or the weight of error_vectors and dec_method, which is same_p or diff_p: 
  //diff p decode: pmin=wt/n*0.5, pmax=wt/n*1.5,   same p decode: p=wt/n
  
  cout<<"your input:"<<endl;
  vector<string> input;
  input.push_back("./CBP");
   input.push_back("H_file");
    input.push_back("L_file");
	 input.push_back("P_file");
	  input.push_back("max_failures/max_trails");
	   input.push_back("max number of iterations");
	    input.push_back("data_file");
		input.push_back("debug");
		input.push_back("alpha");
		input.push_back("other_decoder (start row reduction at the begining or not, it is useless now)");
		input.push_back("M");
		input.push_back("max_row_freeze");
		input.push_back("OSD_order");
  
   for (int i=1;i<argc;i++)
	  {
		  cout<<i<<": "<<input[i]<<"="<<argv[i]<<endl;
	  }
  if (argc!=13)
  {
	  cout<<" need 12 parameters: ./CBP  H_file  L_file P_file max_failed_cws lmax data_file debug alpha  other_decoder M max_row_freeze OSD_order"<<endl;
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
    if ( argv9 >> other_decoder){}
    else
      {
	cout<<"other_decoder should be a double"<<endl;
	return 1;
      }

    istringstream argv10( argv[10] );
    if ( argv10 >> MMM){}
    else
      {
		cout<<"M should be an int"<<endl;
		return 1;
      }
	  
	   istringstream argv11( argv[11] );
    if ( argv11 >> max_row_freeze){}
    else
      {
		cout<<"max_row_freeze should be an int"<<endl;
		return 1;
      }
	  if(max_row_freeze>0){cout<<"use freeze rows methods: try at most "<<freeze_rows_suc<<" times"<<endl;}
	  cout<<"input is ok \n\n"<<endl;
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
	  else {cout<<" OSD_0"<<endl;}



// if (MMM==0){MMM=214748364;}
  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  // int num_of_x_suc_dec=0;//number of Hx successfully decoded results
  
  bool H_suc=false;


  //read the parity check matrix:
  int n1,n,r1,n2,r2;
  
  
  GF2mat H=read_matrix( n1,r1, H_file);
   GF2mat L=read_matrix( n2,r2, L_file);
   
    if (n1!=n2)
  {
	  cout<<" the two matrices donot match"<<endl;
	  cout<<"H_file is "<<H_file<<endl;
	   cout<<"L_file is "<<L_file<<endl;
	    cout<<"H"<<H<<endl;
		  cout<<"L is "<<L<<endl;  
	  return 1;}  

  n=n1;

  int rankH=GF2mat_rank(H);
  bool start_row_reduction=false; 
  //GF2mat G=get_gen(H);
  GF2mat G(1,1);
 // int rankG=GF2mat_rank(G);
   
  nodes  checks[r1];//checks for Hx and Z errors
  nodes  errors[n];

  int E1=0; //number of edges in Hx factor graph, but this parameter is not used in this prog

  initialize_checks (H, checks,  E1);
  initialize_errors(H, errors);
  
  vec p(n);
    // Create an ifstream object for the file
    std::ifstream file(P_file);
    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Could not open the file " << P_file << std::endl;
        return 1;
    }
    // Create a vector to store the doubles
    std::vector<double> numbers;

    // Variable to store each number
    double num;

    // While there is still data to read
    while (file >> num) {
        // Add the number to the vector
        numbers.push_back(num);
    }
	
	//the last entry of numbers is the depolarzing error rate of each gate
	if (numbers.size()!=n+1) {cout<<"error: numbers.size()!=n+1 numbers.size()="<<numbers.size()<<endl;return 0;}
	double gate_p=numbers[n];
for (int i=0;i<n;i++) {p(i)=numbers[i];}
   
cout<<"circuit level decoder, after gates depolarzing rate="<<gate_p<<endl;
  vec LR(n);


  //  cout<<"at the begining, Gz=\n"<<Gz<<"\n Hx=\n"<<Hx<<endl;
 int LLR_fail=0;
 int LLR_suc=0;
 int num_fixed_LLR=0;
 
 
 cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"M="<<MMM<<", lmax="<<lmax<<", alpha="<<alpha<<endl;

     if (debug&2)
	{
	  cout<<"parallel decoding";
	}
      else
	{
	  cout<<"serial V decodeing";
	}
    

  if (debug&8)
    {
      cout<<", apply OSD "<<endl;
    }
	else{ cout<<", didnot apply OSD"<<endl;}
  
	if (debug&16) {cout<<", apply algebraic decoder at the end"<<endl;}
	if (debug&32) {cout<<", apply avgLLR for BP and OSD "<<endl;}
	//else {cout<<"protocol A"<<endl;}
	
   if (other_decoder==1 and MMM>=0){cout<<"start row reduction at the begining, didnot cout num_convergence"<<endl;}
  else if (other_decoder==0 and MMM>=0){cout<<"start row reduction when BP converged to a trapping set"<<endl;}
  else if(MMM<0){cout<<"M<0 donot perform any row reduction just standard BP"<<endl;} 
  if (MMM==-10 and debug&16){cout<<" but apply algebraic decoder at the end"<<endl;}
  if (MMM==-10 and debug&8){cout<<" but apply OSD at the end"<<endl;}
  cout<<"the size of H is"<<r1<<"*"<<n1<<endl;
  
  Sparse_GF2 sparseH(H);
  Sparse_GF2 sparseL(L);
  Sparse_GF2 Sparse_zero_rvec(L.rows(),1);
  GF2mat zero_mat1(H.rows(),1);
  
  std::chrono::high_resolution_clock::time_point prep = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> prep_time = (prep - start)/1000;
	 std::cout << "initialization took " << prep_time.count() << " seconds.\n";
  while (  num_of_failed_cws<max_failed_cws&&num_of_cws<max_num_of_cws)
    {
      num_of_cws++;
	  H_suc=false;
	  
		int wt_real_e=0;
		GF2mat real_e(H.cols(),1); 
		wt_real_e=error_channel(real_e, p);

		if (wt_real_e==0)
		{
			num_no_error++;
			H_suc=true;			
		}
		else
		{	
			H_suc=Fcircuit_decoder(sparseH,sparseL,H, L, G,real_e,checks, errors,p, num_iter, lmax,wt,  debug, OSD_order,Sparse_zero_rvec, zero_mat1);
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
   
    cout<<" p=   "<<gate_p<<", "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws "<<endl;
	cout<<"logical error rate="<<1.0-1.0*num_of_suc_dec/num_of_cws<<endl;
   cout<<"average iterations="<<num_iter/num_of_suc_dec<<"\n\n"<<endl;
   cout<<"syn_fail="<<syn_fail<<endl;
   cout<<"max_fail="<<max_fail<<endl;
   cout<<" BP_suc="<< BP_suc<<endl;
    cout<<" avg_BP_suc="<< avg_BP_suc<<endl;
	cout<<" freeze_rows_suc="<< freeze_rows_suc<<endl;
	
  
   cout<<"OSD_suc(include avg_OSD_suc) ="<<OSD_suc<<endl;
   cout<<"avg_OSD_suc ="<<avg_OSD_suc<<endl;
   cout<<"OSD_0_suc="<<OSD_0_suc<<endl;
   cout<<"OSD_1_suc="<<OSD_1_suc<<endl;
   cout<<"OSD_2_suc="<<OSD_2_suc<<endl;
   cout<<"OSD_higher_order_suc="<<OSD_higher_order_suc<<endl;
   cout<<"OSD_fail="<<OSD_fail<<endl;
   cout<<"num of trivial error is "<<num_no_error<<endl;
   cout<<"LLR_fail="<<LLR_fail<<endl;
   cout<<"LLR_suc="<<LLR_suc<<endl;
   cout<<"num_of_fixed_qubits/num_of_trails="<<1.0*num_fixed_LLR/num_of_cws<<endl;
   cout<<"fix_suc="<<fix_suc<<endl;
    cout<<"fix_fail="<<fix_fail<<endl;

	cout<<"protocol_B_suc="<<protocol_B_suc<<endl;
	cout<<"protocol_B_syn_fail (decoder ouput the error with right syndrome but fails)="<<protocol_B_syn_fail<<endl;
	cout<<"protocol_B_max_fail (decoder did not ouput error with the right syndrome)="<<protocol_B_max_fail<<endl;
	cout<<"num of converging to a trapping set="<<converge_to_trapping_set<<endl;;
	cout<<"num of iterations exceeds half of max iterations (most likely meets an oscillation or maybe it's a matter of calculation accuracy)="<<did_not_converge_for_large_iterations<<endl;;

  
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
  int total_red_qubits=avg_red_qubits;
   

 
   myfile<<n<<" "<<gate_p<<" "<<num_of_cws<<" "<<failed_trials<< " "<<before_OSD_fail<<" "<<OSD_0_fail<<" "<<OSD_1_fail<<" "<<OSD_2_fail<< " "<<OSD_trails<<endl;
  
     //myfile << n<<" "<<d<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<gate_p<<" "<<1.0*num_iter/num_of_suc_dec<<"  "<<num_of_suc_dec<<" "<<num_of_cws<<"  "<<syn_fail<<" "<<max_fail<<" "<<1.0*syn_fail/num_of_cws<<" "<<1.0*max_fail/num_of_cws<<" "<<alpha<<" "<<OSD_suc<<" "<<MMM<<"  "<<LLR_suc<<"  "<<LLR_fail<<endl;
    
 
   myfile.close();
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "\n Run-time " << time_span.count() << " seconds.\n";
     
  return 0;

}
  
