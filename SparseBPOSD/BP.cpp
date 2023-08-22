#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <chrono>
#include <thread>


#include"BP.h"
#include"modified_BP.h"

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

void initialize_checks (const GF2mat &H, nodes checks[], int & E){
  int r=H.rows();
  int c=H.cols();
  
  for (int i=0; i<r;i++)
    {
    checks[i].degree=0;
    for (int j=0; j<c;j++)
      { 
      if (H(i,j)==1)
	{
	checks[i].degree++;
	(checks[i].neighbors).ins(0,j);
	E++;
      }     
    }    
  }
}



//find the neighbors of variable(error) nodes
void initialize_errors(const GF2mat &H, nodes errors[]){

    int r=H.rows();
    int c=H.cols();
    //int index=0;
 
    for (int i=0; i<c;i++)
      {
      errors[i].degree=0;
      for (int j=0; j<r ;j++)
	{
	if (H(j,i)==1)
	  {
	  errors[i].degree++;
	  (errors[i].neighbors).ins(0,j);
	  //index++;	
	}	
      }
    }
}

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H){
    int r=H.rows();
    int c=H.cols();
    for (int i=0;i<r;i++)
      {
      for (int j=0;j<c;j++)
	{
	    if (H(i,j)==1){
	      mcv.set(i,j,1);
	      mvc.set(i,j,1);	 
	     }
	}	     
    }  
}


//parallel schedule for classical BP
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e,int debug){
  double ipr=(1-p)/p;//initial probability ratios.

      //update all variables-to-checks massages:
  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//fix c_i,update all v-to-c_i massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);	
      }
    }
  
    //update all checks-to-variables massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=ipr;
   for (int i=0;i<vj_degree;i++)
     {
       int cnode=(errors[j].neighbors)(i);
       if( (debug/16)%2==1){ update_ci_to_vj_MS( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));}
       else {update_ci_to_vj( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));}
       
       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
      output_e.set(j,0,final_pr<1? 1:0);        
    }  
}


//the serial V schedule 
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat &mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e,int debug)
{
    double ipr=(1-p)/p; 
    double  final_pr;  
    
    //fix j,  for every v_j, do the following:
    for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 final_pr=ipr;

	 //ci is the ith neibor of vj:
	 for (int i=0;i<vj_degree;i++)
	   {
	     //  update the  v_j to c_i massage:      
	     int ci=(errors[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);          
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     if( (debug/16)%2==1){ update_ci_to_vj_MS( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));}
	     else {update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));}
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
}


//the serial C schedule
void s_C_update(const nodes checks[],const nodes errors[],mat &mcv,mat &mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e)
{
  /*
    double ipr=(1-p)/p; 
    double  final_pr;  
    
    //fix j,  for every c_j, do the following:
    for (int j=0;j<c;j++)
      {
	 int cj_degree=checks[j].degree;
	 final_pr=ipr;

	 //ci is the ith neibor of vj:
	 for (int i=0;i<cj_degree;i++)
	   {
	     //  update the  v_i to c_j massage:      
	     int vi=(checks[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,ci,j, ipr);          
	   }
	 
	 for (int i=0;i<cj_degree;i++)
	   {     
	     int vi=(checks[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,j,vi,syndrome(j,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(j,vi);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
  */
}

//update vi to cj massage:
 
      
//clasical decoding
int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,double p,int debug)
{
  int n=v; 
 
  int wt=0;
  GF2mat real_e(v,1);    
  wt=cla_error_channel(real_e, p);//the weight of the error
  //GF2mat zero_vec(v,1);
  
   //if no error, return suc
  if (wt==0)
    {
      return 1;
    }
  
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{	    
	  er=er+ distance(zero_mat2, real_e, n); // for calculating bit error rate after decoding
	  return 0;    //failed		   
	}
   
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      GF2mat output_e(v,1);
   
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  if (debug==0)
	    {
	      p_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e,debug);//parallel schedule
	    }
	  else
	    {
	      s_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e,debug);//parallel schedule
	    }
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	     
	      if(output_e==real_e)
		{
		  //cout<<"success! iteration number="<<l<<endl;
		  return 1;
		  
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  //cout<<"error!,get the wrong e"<<endl;
		  // er=er+ distance(output_e, real_e, n);	  
		  return 0;		       
		}	    	  
	    }
	  if(l==lmax)
	    {
	      //er=er+ distance(output_e, real_e, n);
		//cout<<"real_e"<<endl;   
		return 0;	   	 
	    }
	}
 }

  //the distance between 2 cws
int distance(const GF2mat &output_e, const GF2mat &real_e,int n)
{ 
  int er=0;
  for (int i=0;i<n;i++)
    {
    if( (output_e(i,0)!=real_e(i,0)))
      {
      er++;
      }
    }
  return er;
}


//read a parity check matrix
GF2mat read_matrix (int& n,int &r, string & file_name)
{
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int temp;
  int temp2=0;
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;   
    }
 
  GF2mat H(r,n);
  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp)
	{
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H.set(row_ind,temp-1,1);
	    // cout<< H(row_ind,temp-1)<<endl;
	  }
	else if(temp2==0)
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;
	    temp2++;
	  }	
	}
      row_ind++;
    }
  parity_check.close();
  // cout<<H<<endl;
  return H;  
}

GF2mat read_matrix_from_0 (int& n,int &r, string & file_name)
{
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int temp;
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;   
    }
 
  GF2mat H(r,n);
  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp)
	{
	if (temp>=0&&temp<n&&row_ind<r)
	  {
	    H.set(row_ind,temp,1);
	    // cout<< H(row_ind,temp-1)<<endl;
	  }
	else
	  {
	    cout<<"temp="<<temp<<" n="<<n<<"  row_ind="<<row_ind<<"  r="<<r<<endl;
	    cout<<"the format of the parity check is wrong, the first element is 0 rathar than 1"<<endl;	    
	  }	
	}
      row_ind++;
    }
  parity_check.close();
  // cout<<H<<endl;
  return H;  
}


void write_matrix(string file_name, GF2mat &H)
{
  ofstream Hx;
  Hx.open (file_name,ios::trunc);
  int n=H.cols();
  int r=H.rows();

  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<r;j++)
    {
      for (int i=0; i<n;i++)
	{
	  if (H(j,i)!=0)
	    {
	      Hx<<i+1<<" ";
	    }	  
	}
      Hx<<endl;
    }
  Hx.close();
}


GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right)
{
  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      GF2mat error(1,1); 
      return error;
    }

  else
    {
      int r=left.rows();
      int c=left.cols()+right.cols();
      int c1=left.cols();
      int c2=right.cols();
      GF2mat m(r,c);

      for (int i=0;i<r;i++)
	{
	  for (int j1=0;j1<c1;j1++){
	    m.set(i,j1,left(i,j1));
	  }
	  for (int j2=0;j2<c2;j2++){
	    m.set(i,c1+j2,right(i,j2));
	  }	  
	}
       return m;
    }
}

//the bsc error channel


// Function to generate a batch of random numbers
std::vector<double> generateRandomBatch(int batchSize) 
{
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    std::vector<double> randomNumbers(batchSize);
    for(int i = 0; i < batchSize; ++i) 
    {
        randomNumbers[i] = distribution(generator);
    }
    return randomNumbers;
}

int error_channel(GF2mat &cw, const vec &p)
{
    bin one = 1;
    int temp = 0;
    int r = cw.rows();
    int p_size = p.size();

    if (r != p_size)
    {
        cout << "the size of p and cw do not match" << endl;
        return 0; 
    }

    // Generate random numbers in batch
    std::vector<double> randomNumbers = generateRandomBatch(r);

    for (int i = 0; i < r; i++)
    {
        if (randomNumbers[i] < p[i])
        {
            cw.set(i, 0, cw(i, 0) + one);
            temp++;
        }
    }

    return temp;
}


int depolarizing(GF2mat &cw, double p)
{
  double temp;
  bin one=1;
  int temp2=0;
  int r=cw.rows();
  int temp3=r%3;
  int n=r/3;
  if (temp3!=0){cout<<"depolarizing error, length of e is not 3n"<<endl;return 1;}
  
    for (int i=0;i<n;i++)
      {
		temp=randu();
		//X-error
		if(temp<p)
		{	 
			cw.set(i,0,cw(i,0)+one);
			temp2++;
		}
			//Z-error
		else if (temp>p and temp<2*p)
		{	 
			cw.set(i+n,0,cw(i+n,0)+one);
			temp2++;
		}
			//Y-error
		else if (temp>2*p and temp<3*p)
		{	 
			cw.set(i+2*n,0,cw(i+2*n,0)+one);
			temp2++;
		}
     }
 
  return temp2;
}

int depolarizing2(GF2mat &ex, GF2mat &ez, double p)
{
  double temp;
  bin one=1;
  int temp2=0;
  int n=ex.rows();
 
  if (ex.rows()!=ez.rows()){cout<<"depolarizing channel error, length of ex is not equal to length of ez"<<endl;return 1;}
  
    for (int i=0;i<n;i++)
      {
		temp=randu();
		//X-error
		if(temp<p)
		{	 
			ex.set(i,0,ex(i,0)+one);
			temp2++;
		}
			//Z-error
		else if (temp>p and temp<2*p)
		{	 
			ez.set(i,0,ez(i,0)+one);
			temp2++;
		}
			//Y-error
		else if (temp>2*p and temp<3*p)
		{	 
			ex.set(i,0,ex(i,0)+one);
			ez.set(i,0,ez(i,0)+one);
			temp2++;
		}
     }
 
  return temp2;
}


int cla_error_channel(GF2mat &cw, double p)
  {
  double temp2;
  bin one=1;
  int temp=0;
  int r=cw.rows();
 
    for (int i=0;i<r;i++)
      {
	temp2=randu();
	if(temp2<p)
	  {	 
	    cw.set(i,0,cw(i,0)+one);
	    temp++;
	  }	
      }
    
  return temp;
}
// error channel for fixed weight errors
void error_channel2(GF2mat &error, int wt)
{ 
  double temp2;
  bin one=1;
  int n=error.rows();

    for (int i=0;i<wt;i++)
      {    
	temp2=randi(0,n-1);
	error.set(temp2,0,one);	      
      }
    if (s_weight(error)!=wt)
      
      { GF2mat error2(n,1);
	error=error2;
	//cout<<"!=wt, try again"<<endl;
	error_channel2(error,wt);
      }
}

//weight of a row vector
int weight(const GF2mat &cw)
{
  int n=cw.cols();
  int wt=0;
  for (int i=0;i<n;i++)
    {
      if(cw(0,i)==1){wt++;}
    }
  return wt;
}

//weight of a column vector
int s_weight(const GF2mat &s)
{
  int k=s.rows();
  int wt=0;
  for (int i=0;i<k;i++)
    {
      if(s(i,0)==1){wt++;}
    }
  return wt;
}


//get the error rate distribution
void pro_dist(double pmin,double pmax, vec& pv)
{ 
  double pdiff=pmax-pmin;
  int pvsize=pv.size();
  double temp;
 
  for (int i=0;i<pvsize;i++)
    {
      temp=pdiff*randu();
      pv(i)=pmin+temp;
    }
}

void pro_dist2(double pmin,double pmax, vec& pv)
{ 
  double p_avg=(pmax+pmin)/2;
  int range=(p_avg-pmin)/p_avg;
  int pvsize=pv.size();
 
 
  for (int i=0;i<pvsize;i++)
    {
      double temp=randu();
	  if (temp>0.5){ pv(i)=(1+range)*p_avg;}
	  else{pv(i)=(1-range)*p_avg;}
     
    }
}

bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,int& num_no_error,double alpha,int other_decoder)
{
  int v=H.cols();
  int c=H.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1); 

  if (wt==0)  
    {
      wt_real_e=error_channel(real_e, pv);
      if (wt_real_e==0)
	{
	  num_no_error++;
	  return true;
	}
    }
  else
    {
      error_channel2(real_e,wt);
    }

  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat tempmat333=G*real_e;
  GF2mat zero_rvec2(tempmat333.rows(),1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix

      GF2mat output_e(v,1);
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  if (other_decoder==0)
	    {
	      if ((debug/2)%2==1)
		{
		  quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);
		}
	      else
		{
		  quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);
		}
	    }
	  else if (other_decoder==1){quan_s_C_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);}
	    
	   if (H*output_e==syndrome)
		{		  
		  if(G*(output_e+real_e)==zero_rvec2)
		    {
		      num_iter= num_iter+l;		  
		      return true;
		    }
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
		      syn_fail++;		
		      return false;      	        
		    }	    	  
		}  
	}
      	            
      if (debug%2==1)
	{
	  mcv.zeros();
	  mvc.zeros();
	  initialize_massages( mcv,mvc, H);
	  vec pv2(v);
	 
	  // GF2mat output_e2(v,1);
	  // GF2mat syndrome2=H*(real_e+output_e);
	  //int wt_syn=s_weight(syndrome);
	  //int wt_new_syn=s_weight(syndrome2);
	  vec LR2=LR;
	  int l2=0;
	  while (l2<2)
	    {
	      pro_dist(pmin,pmax,pv2);
	      for (int l=1;l<=lmax;l++)
		{
		  if (other_decoder==0)
		    {
		      if ((debug/2)%2==1)
			{
			  quan_p_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,LR2,debug,alpha);
			}
		      else
			{
			  quan_s_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,LR2,debug,alpha);
			}
		    }
		  else if (other_decoder==1){quan_s_C_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,LR,debug,alpha);}
		  
		  if (H*output_e==syndrome)
		    {		  
		      if(G*(output_e+real_e)==zero_rvec2)
			{
			  num_iter= num_iter+l;		  
			  return true;
			}
		      else
			{
			  // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
			  syn_fail++;		
		      return false;      	        
			}	    	  
		    }		 
		}
	      l2++;
	      // output_e=output_e+output_e2;
	      //syndrome2=H*(real_e+output_e);
	      //wt_syn=wt_new_syn;
	      //wt_new_syn=s_weight(syndrome2);
	    }
	  // output_e=output_e+output_e2;
	}
   
      if ((debug/8)%2==1)
	{
		int suc_order=0;
	  OSD(LR,H,syndrome,output_e,suc_order);
	  if(G*(output_e+real_e)==zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;
	      return true;
	    }	 
	}
      
       if ((debug/4)%2==1)
	 {
	   cout<<"output e is \n"<<endl;
	   err_pos2(output_e);
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<LR<<endl;
	 }
      max_fail++;
      return false;
 }


bool  OSD_decode(GF2mat &H, GF2mat &G,int &wt, vec& pv,int debug, vec &LR,int rankH,int OSD_order)
{
  int v=H.cols();
  int c=H.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1); 
  

  if (wt==0)  
    {
      wt_real_e=error_channel(real_e, pv);
      if (wt_real_e==0)
	{
	  num_no_error++;
	  return true;
	}
    }
  else
    {
      error_channel2(real_e,wt);
    }

  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat tempmat333=G*real_e;
  GF2mat zero_rvec2(tempmat333.rows(),1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	
   GF2mat output_e(v,1);
	int suc_order=0;
	
	   full_OSD(LR, H, syndrome,output_e,OSD_order,suc_order);
	
	  if(G*(output_e+real_e)==zero_rvec2)
	    {	     
	      OSD_suc++;
		  if (suc_order==0){OSD_0_suc++;}
		  else if (suc_order==1){OSD_1_suc++;}
		  else if (suc_order==2){OSD_1_suc++;}
		  else if (suc_order>2){OSD_higher_order_suc++;}
	      return true;
	    }	 

      
       if ((debug/4)%2==1)
	 {
	   cout<<"output e is \n"<<endl;
	   err_pos2(output_e);
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<LR<<endl;
	 }
      max_fail++;
      return false;
 }
//ordered statistical decoder
void full_OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e,int OSDorder,int& suc_order) //r is the rank of H
{ 
    //get the first permutation that abs_LLR is in descending order
  int n=LR.length();
  if (n!=H.cols()){cout<<"OSD error: LR.length!=H.cols"<<endl;}
  int r=H.rows();
  vec LLR(n); // LR=p0/p1
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  LLR(i)=log(LR(i));
	}
      else
	{
	  cout<<"negative LR!"<<endl;
	  cout<<LR(i)<<endl;
	}
    }
  
  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1
  ivec perm1=sort_index(LLR); //the sort function gives a ascending  order
  //for (int i=0;i<n;i++){perm1(i)=i;}
  GF2mat H1=H;
  H1.permute_cols(perm1,0);
  GF2mat perm1_mat=col_permutation_matrix(perm1);
  GF2mat perm1_mat_inv=perm1_mat.inverse();
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*perm1_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*perm1_inv*e:

  int rankH=H1.T_fact(T,U,perm2);
  if (rankH==0){return;}
  GF2mat perm2_mat=col_permutation_matrix(perm2);
  GF2mat perm2_mat_inv=perm2_mat.inverse();
      
  GF2mat syndrome1=T*syndrome;
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);
  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
  
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  GF2mat H_T(1,1);
   
 if (H2.cols()>rankH){H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);}
  GF2mat HS_inv=H_S.inverse();
  
  GF2mat e_S=HS_inv*syndrome2;  
  GF2mat e_T(1,1);
  if(n-rankH>0){ GF2mat e_T_temp(n-rankH,1); e_T=e_T_temp;}
  
  GF2mat new_e_S;
  int min_wt=s_weight(e_S);
  //cout<<"original wt "<<wt<<endl;
  int temp1;
  int temp2=-1;
  int temp3=-1;
 // int lambda=min(100,n-rankH);

  if (n-rankH>0){
		GF2mat new_e_T=e_T;
		int index=0;
		GPT_dfs_For_OSD(index, e_T.rows(), OSDorder, 0, suc_order, new_e_T, e_T, e_S, min_wt, HS_inv, H_T,syndrome2);
		e_S=HS_inv*syndrome2+HS_inv*H_T*e_T;
    }       
		for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
		for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}
		output_e=perm1_mat*perm2_mat*output_e;
   	  	  
		if (H*output_e==syndrome){}
		else
		{
			cout<<"OSD error H*output_e!=syndrom"<<endl;
			cout<<"suc_order="<<suc_order<<endl;
		}  
   
}			 

void GPT_dfs_For_OSD(int index, int n, int OSDorder, int non_zero_count, int& suc_order,GF2mat& new_e_T, GF2mat& e_T,const GF2mat& e_S, int& min_wt, const GF2mat& HS_inv, const GF2mat&H_T,const GF2mat& syndrome2) {
    if (non_zero_count > OSDorder) {
        return;
    }

    if (index == n) {
        
        GF2mat new_e_S=HS_inv*syndrome2+HS_inv*H_T*new_e_T;
		int new_wt=non_zero_count+s_weight(new_e_S);
		if (new_wt<min_wt)
		{
			//cout<<"new e_T"<<endl;
			min_wt=new_wt;
			e_T=new_e_T;
			suc_order=non_zero_count;
		}
       return;
    }

    for (int i = 0; i < 2; i++) {
        new_e_T.set(index,0,i);
        int new_non_zero_count = non_zero_count + (i == 1 ? 1 : 0);

        if (new_non_zero_count <= OSDorder) {
            GPT_dfs_For_OSD(index + 1, n, OSDorder,  new_non_zero_count,suc_order,new_e_T, e_T,e_S, min_wt,HS_inv,H_T,syndrome2);
        }
    }
}

void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e,int& suc_order) //r is the rank of H
{ 
    //get the first permutation that abs_LLR is in descending order
  int n=LR.length();
  if (n!=H.cols()){cout<<"OSD error: LR.length!=H.cols"<<endl;}
  int r=H.rows();
  vec LLR(n); //the sort function gives a ascending  order, LR=p0/p1
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  LLR(i)=log(LR(i));
	}
      else
	{
	  cout<<"negative LR!"<<endl;
	  cout<<LR(i)<<endl;
	}
    }
  
  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1

	ivec perm1=sort_index(LLR);
	global_sort_index=perm1;
	//for (int i=0;i<n;i++){perm1(i)=i;}
	GF2mat H1=H;
	H1.permute_cols(perm1,0);
	GF2mat perm1_mat=col_permutation_matrix(perm1);
	GF2mat perm1_mat_inv=perm1_mat.inverse();
      
	GF2mat T;
	ivec perm2;
	GF2mat U;
	// now we have T*H1*perm2*perm2_inv*perm1_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*perm1_inv*e:
	//There is no T if H.cols()>H.rows()

	int rankH=H1.T_fact(T,U,perm2);

    if (rankH==0){return;}
   //if (rankH<r-1) {cout<<"rankH<H.rows()-1, I have not solved this case yet, stop OSD"<<endl;return;}
   
 //GF2mat H2=H1.get_submatrix(0,0,rankH-1,n-1);
 GF2mat H2=H1;
  GF2mat syndrome2=syndrome;
  GF2mat  perm2_mat;
 if(r<n)
 {

	perm2_mat=column_gaussian(H2,rankH) ;
	H2=H2*perm2_mat;
 
 // global_perm2=perm2; 
  
GF2mat row_perm=row_gaussian(H2,rankH);
H2=row_perm*H2;
  syndrome2=row_perm*syndrome2;
  
  syndrome2=syndrome2.get_submatrix(0,0,rankH-1,0);

 }
 else
 {

	 H2=U;
	 syndrome2=(T*syndrome).get_submatrix(0,0,rankH-1,0);	 
	 perm2_mat=col_permutation_matrix(perm2);
	
 }

  GF2mat perm2_mat_inv=perm2_mat.inverse();
  
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
 
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  
  
  ivec temp_perm;
  int HS_rank=H_S.T_fact(T,U,temp_perm);
  if (HS_rank!=rankH)
  {
	  cout<<"error:HS_rank!=rankH, H2 is\n"<<H2<<"\nHS is \n"<<H_S<<"\n H1 is\n"<<H1<<endl;
  }
  
  GF2mat H_T(1,1);
   
 if (H2.cols()>rankH){H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);}
  GF2mat HS_inv=H_S.inverse();

  GF2mat e_S=HS_inv*syndrome2;  
  global_e_S=e_S;
  GF2mat e_T(1,1);
  if(n-rankH>0){ GF2mat e_T_temp(n-rankH,1); e_T=e_T_temp;}
 
  GF2mat new_e_S;
  int wt=s_weight(e_S);
  //cout<<"original wt "<<wt<<endl;
  int temp1;
  int temp2=-1;
  int temp3=-1;
  int  lambda=min(100,n-rankH);

   if(n-rankH>0) 
   {
	   for (int i=0;i<n-rankH;i++)
		{
			//cout<<i<<": wt is "<<wt<<endl;
			GF2mat new_e_T=e_T;
			new_e_T.set(i,0,1);
			new_e_S=HS_inv*syndrome2+HS_inv*H_T*new_e_T;
			temp1=s_weight(new_e_S);
			if (temp1+1<wt)
			{
				wt=temp1+1;
				temp2=i;
			}
		}


		for (int i=0;i<lambda;i++)
		{
			for (int j=i+1;j<lambda;j++)
			{
				GF2mat new_e_T=e_T;
				new_e_T.set(i,0,1);
				new_e_T.set(j,0,1);		  
				new_e_S=HS_inv*syndrome2+HS_inv*H_T*new_e_T;
				temp1=s_weight(new_e_S);

				if (temp1+2<wt)
				{
					// cout<<"another e_S"<<endl;
					wt=temp1+2;
					// cout<<"new wt is "<<wt<<endl;
					temp2=i;
					temp3=j;
				}
			}
		}
 
		if (temp2!=-1&&temp3==-1)
		{
			e_T.set(temp2,0,1);
			suc_order=1;
			e_S=HS_inv*syndrome2+HS_inv*H_T*e_T;
			global_new_wt=1+s_weight(e_S);
		}
		else if (temp2!=-1&&temp3!=-1)
		{
			e_T.set(temp2,0,1);
			e_T.set(temp2,0,1);
			suc_order=2;
			e_S=HS_inv*syndrome2+HS_inv*H_T*e_T;
			global_new_wt=2+s_weight(e_S);
		}
   }
         
		for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
		for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}
		output_e=perm1_mat*perm2_mat*output_e;
   	  	  
		if (H*output_e==syndrome){}
		else
		{
			cout<<"OSD error H*output_e!=syndrom"<<endl;
		}  
   
}			 


	
//get a column permutation matrix
GF2mat col_permutation_matrix(ivec & perm)
{
  int n=perm.length();
  GF2mat p(n,n);
  for ( int i=0;i<n;i++)
    {
      p.set(perm(i),i,1);
    }
  return p;
}

/*
GF2mat col_permutation_matrix_s(ivec & perm)
{
  int n=perm.length();
  GF2mat p(n,n);
  for ( int i=0;i<n;i++)
    {
      p.set(perm(i),i,1);
    }
  return p;
}
*/

void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha)
{  
  double ipr;

      //fix vj, update all c to vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=(1-pv[j])/pv[j];
   
   for (int i=0;i<vj_degree;i++)
     {
       int cnode=(errors[j].neighbors)(i);
       if( (debug/16)%2==1){ update_ci_to_vj_MS( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));}
       else {update_ci_to_vj( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));}
       final_pr=final_pr*pow(mcv(cnode,j),1.0/alpha);
      }   
   LR(j)=final_pr;
   output_e.set(j,0,final_pr<1? 1:0);   
    }

  
  // fix ci, update all v to ci messages
  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;      
      //update all v-to-ci messages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        ipr=(1-pv[vnode])/pv[vnode];
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr,alpha);	
      }
    }

}

// serial V update:
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha)
{  
    double ipr;
    double final_pr;
   for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 ipr=(1-pv[j])/pv[j];
	 final_pr=ipr;
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
		
	     if( (debug/16)%2==1){ update_ci_to_vj_MS( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));}
	     else {update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));}	 
	     final_pr=final_pr*mcv(ci,j);      	
	   } 
	 //ci is the ith neibor of vj:
	 for (int i=0;i<vj_degree;i++)
	   {
	     //  update the  v_j to c_i massage:      
	     int ci=(errors[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);     		
	   }
	 //  cout<<j<<"   "<<final_pr<<endl;
	 LR(j)=final_pr;
	
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
}

// serial C update:
void quan_s_C_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha)
{  
    double ipr;
    double final_pr;
    
   for (int j=0;j<c;j++)
      {
	 int cj_degree=checks[j].degree;

	 //vi is the ith neibor of cj:
	 for (int i=0;i<cj_degree;i++)
	   {
	     //  update the  v_i to c_j massage:      
	     int vi=(checks[j].neighbors)(i);
	     ipr=(1-pv[vi])/pv[vi];
	     update_vj_to_ci(checks, errors,mcv, mvc,vi,j, ipr);          
	   }	 
	 for (int i=0;i<cj_degree;i++)
	   {     
	     int vi=(checks[j].neighbors)(i);
	     if( (debug/16)%2==1){ update_ci_to_vj_MS( checks, errors,mcv, mvc,j,vi,syndrome(j,0));}
	     else {update_ci_to_vj( checks, errors,mcv, mvc,j,vi,syndrome(j,0));}     
     
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
   
      }
   
   //calculate LR:
   for (int j=0;j<v;j++)
     {
       ipr=(1-pv[j])/pv[j];
       final_pr=ipr;
       int vj_degree=errors[j].degree;
       for (int i=0;i<vj_degree;i++)
	 {
   
	   int ci=(errors[j].neighbors)(i);
	   final_pr=final_pr*mcv(ci,j);
	 }
       LR(j)=final_pr;
       output_e.set(j,0,final_pr<1? 1:0);   
     }
}

GF2mat get_gen(const GF2mat &H){

  GF2mat HT=H.transpose();
  GF2mat T,U;
  ivec P;
  int Hrank= HT.T_fact(T,U,P);
  int r=H.rows();
  int n=H.cols();
  GF2mat G=T.get_submatrix(Hrank,0,n-1,n-1);
  return G;
}

GF2mat get_gen(const GF2mat &Hx, const GF2mat &Hz){

   GF2mat H_star=get_gen(Hx);
   GF2mat L;

	get_logical(Hx, Hz,H_star,L);
	GF2mat G=merge_mat_vert(Hz,L);
 
	return G;
}
/*
void dense_to_sparse(GF2mat &G,GF2mat& G_s)
{
  
  int r=G.rows();
  int c=G.cols();
  for (int i=0;i<r;i++)
    {
      for (int j=0;j<c;j++)
	{
	  if (G(i,j)==1){G_s.set(i,j,1);}	
	}
    }
}
*/
void err_pos1(const GF2mat &error)
{
  cout<<"wt="<<s_weight(error)<<endl;
    int n=error.rows();
    for (int i=0;i<n;i++)
      {
	if (error(i,0)==1)
	  {
	    cout<<i<<" ";
	  }
      }
    cout<<endl;
}
void err_pos2(const GF2mat &error)
{
  // cout<<"\n"<<error<<endl;
  cout<<"\n";
  int n=error.rows();
  int d=sqrt(n);
  if (d*d==n)
    {
      for (int r=0;r<d;r++)
	{
	  cout<<r<<"th row: ";
	  for (int c=0;c<d;c++)	    
	    {
	      if (error(r*d+c,0)==1){cout<<1;}
	      else {cout<<".";}
	    }
	  cout<<endl;
	}
    }
    cout<<"\n";
}

int GF2mat_rank(const GF2mat& H){

  GF2mat T,U;
  ivec P;
  return H.T_fact(T,U,P);	
}

bmat GF2mat_to_bmat(const GF2mat& H)
{
  bmat Hb(H.rows(),H.cols());
  for (int i=0;i<H.rows();i++)
    {
      for (int j=0;j<H.cols();j++)
	{
	  if (H(i,j)==1){Hb.set(i,j,1);}
	  else {Hb.set(i,j,0);}
	}
    }
  return Hb;
}

GF2mat merge_mat_vert(const GF2mat &up,const GF2mat &bottom)
{
  if (up.cols()!=bottom.cols())
    {
      cout<<"the 2 matrices cannot merge vertically"<<endl;
      GF2mat error(1,1); 
      return error;
    }

  else
    {
      int c=up.cols();
      int r=up.rows()+bottom.rows();
      int r1=up.rows();
      int r2=bottom.rows();
      GF2mat m(r,c);

      for (int i=0;i<c;i++)
	{
	  for (int j1=0;j1<r1;j1++)
	    {
	      m.set(j1,i,up(j1,i));
	    }
	  for (int j2=0;j2<r2;j2++)
	    {
	      m.set(j2+r1,i,bottom(j2,i));
	    }	  
	}
      //if (debug==1){cout<<"up:\n"<<up<<"\n  bottom:\n"<<bottom<<"\n result: \n"<<m<<endl;}
      return m;
    }
}

void get_logical(const GF2mat& H, const GF2mat& D,const GF2mat& H_star, GF2mat& L)
{
  vector<GF2mat> Lo;
  GF2mat T;
  ivec perm;
  GF2mat U;
  int prerank=D.T_fact(T,U,perm);
  int temprank=0;
  int first_lo=0;
  GF2mat tempmat2=D;

  for (int i=0;i<H_star.rows();i++)
    {
      GF2mat tempmat1=H_star.get_row(i).transpose();
      tempmat2=merge_mat_vert(tempmat2,tempmat1);
      temprank=tempmat2.T_fact(T,U,perm);
      if (temprank==prerank+1){Lo.push_back(tempmat1);prerank=temprank;}
      else if(temprank!=prerank) {cout<<"sth wrong with the get_locial operator"<<endl;return;}
    }
  L=Lo[0];

  for (int j=1;j<Lo.size();j++)
    {
      L=merge_mat_vert(L,Lo[j]);
    }
  
    
}


bool  decoder_T(GF2mat &H, GF2mat &HT, GF2mat &G, GF2mat &new_H,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,double& num_iter, int lmax,int &wt,int debug, vec &LR,int rankH)
{
  int v=H.cols();
  int c=H.rows();
  double alpha=1;
  int wt_real_e=0;
  GF2mat real_e(v,1); 
//cout<<11<<endl;
  if (wt==0)  
    {
      wt_real_e=error_channel(real_e, pv);
      if (wt_real_e==0)
	{
	  num_no_error++;
	  return true;
	}
    }
  else
    {
      error_channel2(real_e,wt);
    }

  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat tempmat333=G*real_e;
  GF2mat zero_rvec2(tempmat333.rows(),1);
  
  GF2mat syndrome=H*real_e;
 

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  

      mat mcv(c,c);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,c);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, new_H); //initialize to all-1 matrix

      GF2mat output_e(v,1);
	  GF2mat new_e(c,1);
      
	  GF2mat new_H_inv=new_H.inverse();
 GF2mat new_e2=new_H_inv*syndrome;
 GF2mat output_e2=HT*new_e2;
 
   if (H*output_e2==syndrome)
	{		  
		  if(G*(output_e2+real_e)==zero_rvec2)
		    {
		      cout<<"suc"<<endl;
		      return true;
		    }
		  else
		    {
		     
				syn_fail++;		
				if (debug&4)
				{
					cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
					cout<<"sync failed:"<<endl;
					cout<<"output e2 is \n"<<endl;
					
					cout<<"output e2 is \n"<<endl;
					err_pos2(output_e2);
					cout<<"new e is \n"<<new_e2.transpose()<<endl;
				}
			}
	
	}
 
      for (int l=1;l<=lmax;l++)
	{
	    //cout<<l<<endl;

//cout<<1<<endl;
		  quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, c,new_e,LR,debug,alpha);
		 // cout<<2<<endl;
			output_e=HT*new_e;
	    
	   if (H*output_e==syndrome)
		{		  
		  if(G*(output_e+real_e)==zero_rvec2)
		    {
		      num_iter= num_iter+l;		  
		      return true;
		    }
		  else
		    {
		     
		      syn_fail++;		
			  if (debug&4)
		{
		 cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		 cout<<"sync failed:"<<endl;
		cout<<"output e is \n"<<endl;
		err_pos2(output_e);
		cout<<"new e is \n"<<new_e.transpose()<<endl;
		
		
		cout<<"output e2 is \n"<<endl;
		err_pos2(output_e2);
		cout<<"new e is \n"<<new_e2.transpose()<<endl;
		
	   
		cout<<"real e is \n"<<endl;
		err_pos2(real_e);
		GF2mat sume=real_e+output_e;
		cout<<"residual e is \n"<<endl;
		err_pos2(sume);
	
		}
		      return false;      	        
		    }	    	  
		}  
	}
      	    




			
		if (debug&4)
	 {
		
		 cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		 cout<<"sync failed:"<<endl;
		cout<<"output e is \n"<<endl;
		err_pos2(output_e);
		cout<<"new e is \n"<<new_e.transpose()<<endl;
		
		
		cout<<"output e2 is \n"<<endl;
		err_pos2(output_e2);
		cout<<"new e2 is \n"<<new_e2.transpose()<<endl;
		
	   
		cout<<"real e is \n"<<endl;
		err_pos2(real_e);
		GF2mat sume=real_e+output_e;
		cout<<"residual e is \n"<<endl;
		err_pos2(sume);
	
	 }
	   
      max_fail++;
      return false;
 }