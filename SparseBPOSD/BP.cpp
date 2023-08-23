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


// Function to generate a batch of random numbers

std::vector<double> generateRandomBatch(int batchSize)
{
    // Get a seed value from the current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    
    // Initialize the random engine with the seed
    std::default_random_engine generator(seed);
    
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
	
	//cout<<"rightOSD: after col_gaussian H is "<<H2<<endl;
 
 // global_perm2=perm2; 
  
GF2mat row_perm=row_gaussian(H2,rankH);

H2=row_perm*H2;
//cout<<"rightOSD: after row_gaussian H is "<<H2<<endl;

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


