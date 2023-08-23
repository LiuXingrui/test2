#include "Sparse.h"
#include"BP.h"
#include"modified_BP.h"

Sparse_GF2::Sparse_GF2() : num_rows(0), num_cols(0) {}
	
Sparse_GF2::Sparse_GF2(int a, int b) : num_rows(a), num_cols(b) {
        row_mat.resize(a);  // Each row is empty, representing all zeros
        col_mat.resize(b);  // Each column is empty, representing all zeros
    }

    // Constructor that converts a GF2mat to a SparseBinaryMatrix
  Sparse_GF2::Sparse_GF2(const itpp::GF2mat& mat) {
        num_rows = mat.rows();
        num_cols = mat.cols();

        for (int i = 0; i < num_rows; i++) {
            std::vector<int> temp;
            for (int j = 0; j < num_cols; j++) {
                if (mat.get(i, j) == 1) {
                    temp.push_back(j);
                }
            }
            row_mat.push_back(temp);
        }
        row_to_col();       
    }

    // Getter methods
    int Sparse_GF2::rows() const {
        return num_rows;
    }

    int Sparse_GF2::cols() const {
        return num_cols;
    }
	
	std::vector<int> Sparse_GF2::get_row(const int& i) const {
        if (i < 0 || i >= num_rows) {
            throw std::out_of_range("Row index out of bounds.");
        }
        return row_mat[i];
    }

    // Return the jth column of the matrix
    std::vector<int> Sparse_GF2::get_col(const int& j) const {
        if (j < 0 || j >= num_cols) {
            throw std::out_of_range("Column index out of bounds.");
        }
        return col_mat[j];
	}

Sparse_GF2 Sparse_GF2::operator*(const Sparse_GF2& other) const {
    if (num_cols != other.num_rows) {
        throw std::runtime_error("Matrix dimensions mismatch for multiplication.");
    }
    Sparse_GF2 result;
    result.num_rows = num_rows;
    result.num_cols = other.num_cols;

    // Initialize result's rows with all zeros
    result.row_mat.resize(num_rows);

    for (int i = 0; i < num_rows; i++) {
        // This will store the number of times each column is marked
        std::vector<int> col_marks(result.num_cols, 0);
       
        for (int k : row_mat[i]) {
            for (int l : other.row_mat[k]) {
                col_marks[l]++;
            }
        }
        
        // Construct the resulting row based on col_marks
        for (int j = 0; j < result.num_cols; j++) {
            if (col_marks[j] % 2 == 1) {
                result.row_mat[i].push_back(j);
            }
        }
    }

    // Update the col_mat of the result
	result.row_to_col();
    return result;
}

    // Equality check
    bool Sparse_GF2::operator==(const Sparse_GF2& other) const {
        if (num_rows != other.num_rows || num_cols != other.num_cols) {
            return false;
        }
        for (int i = 0; i < num_rows; i++) {
            if (row_mat[i] != other.row_mat[i]) {
                return false;
            }
        }
        return true;
    }

    bool Sparse_GF2::operator!=(const Sparse_GF2& other) const {
        return !(*this == other);
    }

    Sparse_GF2 Sparse_GF2::operator+(const Sparse_GF2& other) const {
        if (num_rows != other.num_rows || num_cols != other.num_cols) {
            // Throw an error or handle the dimension mismatch
            throw std::runtime_error("Matrix dimensions mismatch for addition.");
        }

        Sparse_GF2 result;
        result.num_rows = num_rows;
        result.num_cols = num_cols;
		 // Reserve memory
    result.row_mat.reserve(num_rows);

        for (int i = 0; i < num_rows; i++) {
            std::vector<int> row_result;
            
            // XOR the rows using set_symmetric_difference
            std::set_symmetric_difference(
                row_mat[i].begin(), row_mat[i].end(),
                other.row_mat[i].begin(), other.row_mat[i].end(),
                std::back_inserter(row_result)
            );
            //result.row_mat.push_back(row_result);
			 result.row_mat.emplace_back(row_result);
        }

        // Update the col_mat of the result
		result.row_to_col();
        return result;
    }


void Sparse_GF2::row_to_col() {
        col_mat.clear();
        col_mat.resize(num_cols);

        for (int i = 0; i < num_rows; i++) {
            for (int col_idx : row_mat[i]) {
                col_mat[col_idx].push_back(i);
            }
        }
    }

    // Rebuild row_mat from col_mat
    void Sparse_GF2::col_to_row() {
        row_mat.clear();
        row_mat.resize(num_rows);

        for (int i = 0; i < num_cols; i++) {
            for (int row_idx : col_mat[i]) {
                row_mat[row_idx].push_back(i);
            }
        }
    }

    // Swap two rows of the matrix
    void Sparse_GF2::swap_rows(const int& i, const int& j,const bool& update_col) {
        if (i < 0 || i >= num_rows || j < 0 || j >= num_rows) {
            // Handle index out of bounds
            throw std::out_of_range("Row indices out of bounds.");
        }
        std::swap(row_mat[i], row_mat[j]);

        // Adjust column representations
		if (update_col==true){this->row_to_col();}
       
    }

    // Swap two columns of the matrix
    void Sparse_GF2::swap_cols(const int& i, const int& j,const bool&update_row) {
        if (i < 0 || i >= num_cols || j < 0 || j >= num_cols) {
            // Handle index out of bounds
            throw std::out_of_range("Column indices out of bounds.");
        }
        std::swap(col_mat[i], col_mat[j]);

        // Adjust row representations
        if (update_row==true){this->col_to_row();}
       
    }

    itpp::GF2mat  Sparse_GF2::to_GF2mat() const {
        itpp::GF2mat mat(num_rows, num_cols);
        
        for (int i = 0; i < num_rows; i++) {
            for (int col_idx : row_mat[i]) {
                mat.set(i, col_idx, 1);
            }
        }
        return mat;
    }


    void Sparse_GF2::add_rows(const int& i,const int& j,const bool& update_col) {
    
        if (i < 0 || i >= num_rows || j < 0 || j >= num_rows) {
            // Handle index out of bounds
            throw std::out_of_range("Row indices out of bounds.");
        }
 
        std::vector<int> row_result;
         	
        std::set_symmetric_difference(
            row_mat[i].begin(), row_mat[i].end(),
            row_mat[j].begin(), row_mat[j].end(),
            std::back_inserter(row_result)
        );
        
        row_mat[j] = row_result;
        if (update_col==true){this->row_to_col();}
    }

    // Add the ith column to the jth column
    void Sparse_GF2::add_cols(const int& i, const int& j,const bool& update_row) {
        if (i < 0 || i >= num_cols || j < 0 || j >= num_cols) {
            // Handle index out of bounds
            throw std::out_of_range("Column indices out of bounds.");
        }
        std::vector<int> col_result;
        std::set_symmetric_difference(
            col_mat[i].begin(), col_mat[i].end(),
            col_mat[j].begin(), col_mat[j].end(),
            std::back_inserter(col_result)
        );
      col_mat[j] = col_result;
      if (update_row==true){this->col_to_row();}
    }


Sparse_GF2 Sparse_GF2::row_gaussian(int& rankH) const {
		
        Sparse_GF2 H = *this;
        Sparse_GF2 permutation;
        permutation.num_rows = num_rows;
        permutation.num_cols = num_rows;
        int r=num_rows;

        // Initialize the permutation matrix as an identity permutation
        for (int i = 0; i < num_rows; ++i) {
            permutation.row_mat.push_back({i});
            permutation.col_mat.push_back({i});
        }

        int pivot_row = 0;
        int pivot_col=-1;
        while (pivot_row<r)
        {
        	if (H.row_mat[pivot_row].size()!=0)
			{
				pivot_col=H.row_mat[pivot_row][0];
				
				for (int i=0;i<r;i++)
					{	if(i!=pivot_row)
						{	
							if (find(H.row_mat[i].begin(), H.row_mat[i].end(), pivot_col) != H.row_mat[i].end()) 
							{
								//false means donot update column representation here
								H.add_rows(pivot_row,i,false);
								permutation.add_rows(pivot_row,false);
							}
						}
					}	
        	}       
        pivot_row++;       
    }
	
    rankH=0;
    for (int i=0;i<num_rows;++i){if(H.row_mat[i].size()!=0){rankH++;}}
    int current_row=0;
    int swaped_row=rankH;

     if (rankH<r)
    {
   	 while (current_row<rankH and swaped_row<r)
		{
		if (H.row_mat[current_row].size()==0)
			{
				H.swap_rows(current_row, swaped_row,false);
				permutation.swap_rows(current_row, swaped_row,false);
				swaped_row++;
			}	

			else {current_row++;}
		}
    }
	//update column representation at the end
	permutation.row_to_col();
	  return permutation;
}
	
Sparse_GF2 Sparse_GF2::inverse() const {
    
	if (num_rows!= num_cols) {
            throw std::runtime_error("inverse: need to be a square matrix");
        }
		
        Sparse_GF2 H = *this;
        Sparse_GF2 permutation;
        permutation.num_rows = num_rows;
        permutation.num_cols = num_rows;
        int r=num_rows;

        // Initialize the permutation matrix as an identity permutation
        for (int i = 0; i < num_rows; ++i) {
            permutation.row_mat.push_back({i});
            permutation.col_mat.push_back({i});
        }

        int pivot_row = 0;
        int pivot_col=-1;
        while (pivot_row<r)
        {
        	if (H.row_mat[pivot_row].size()!=0){pivot_col=H.row_mat[pivot_row][0];}
			else{		
            throw std::runtime_error("inverse: there are empty rows, this is not an invertible matrix");
			}
        
        	for (int i=0;i<r;i++)
        		{	if(i!=pivot_row)
        			{	
        				if (find(H.row_mat[i].begin(), H.row_mat[i].end(), pivot_col) != H.row_mat[i].end()) 
        				{
							//I updated the swap functions, false means donot update column representation here
        					H.add_rows(pivot_row,i,false);
        					permutation.add_rows(pivot_row,i,false);
        				}
        			}
        		}
        
        pivot_row++;    
    }
	
    int rankH=0;
	
    for (int i=0;i<num_rows;++i){if(H.row_mat[i].size()!=0){rankH++;}}
	if (rankH!= num_cols) {
            throw std::runtime_error("inverse: this is not an invertible matrix:");
        }
    int current_row=0;
 
   	 while (current_row<rankH)
		{
		if (find(H.row_mat[current_row].begin(), H.row_mat[current_row].end(), current_row)==H.row_mat[current_row].end())
			{
				for (int i=current_row+1;i<rankH;++i)
				{
					if (find(H.row_mat[i].begin(), H.row_mat[i].end(), current_row)!=H.row_mat[i].end())
					{
						H.swap_rows(current_row, i,false);
						permutation.swap_rows(current_row, i,false);
						break;
					}
				}	
			}	
			current_row++;
		}

	//update column representation at the end
	permutation.row_to_col();
  return permutation;
}

 Sparse_GF2 Sparse_GF2::col_gaussian(int& rankH) const {
    
        Sparse_GF2 H = *this;
        Sparse_GF2 permutation;
        permutation.num_rows = num_cols;
        permutation.num_cols = num_cols;
        int c=num_cols;

        // Initialize the permutation matrix as an identity permutation
        for (int i = 0; i < num_cols; ++i) {
            permutation.row_mat.push_back({i});
            permutation.col_mat.push_back({i});
        }
		//Hwrong=permutation.to_GF2mat();

        int pivot_row = -1;
        int pivot_col=0;
        while (pivot_col<c)
        {
        	if (H.col_mat[pivot_col].size()!=0)
			{
				 //pivot_row=*std::min_element((H.col_mat[pivot_col]).begin(), (H.col_mat[pivot_col]).end());
				pivot_row=H.col_mat[pivot_col][0];
					
        	for (int i=0;i<c;i++)
        		{	if(i!=pivot_col)
        			{	
        				if (find(H.col_mat[i].begin(), H.col_mat[i].end(), pivot_row) != H.col_mat[i].end()) 
        				{
        					H.add_cols(pivot_col,i,false);
        					permutation.add_cols(pivot_col,i,false);
        				}
        			}
        		}    
			}				
        pivot_col++;     
    }
	
    rankH=0;
    for (int i=0;i<num_cols;++i){if(H.col_mat[i].size()!=0){rankH++;}}
    int current_col=0;
    int swaped_col=rankH;
    
    if (rankH<c)
    {
   	while (current_col<rankH and swaped_col<c)
		{
			if (H.col_mat[current_col].size()==0)
			{
				H.swap_cols(current_col, swaped_col,false);
				permutation.swap_cols(current_col, swaped_col,false);
				swaped_col++;
			}
			else {current_col++;}
		}
	}

permutation.col_to_row();
  return permutation;
}

Sparse_GF2 Sparse_GF2::permute_cols(const itpp::ivec& Perm) {
        if (Perm.length() != num_cols) {
            throw std::runtime_error("Permutation vector size mismatch.");
        }
		Sparse_GF2 permutation;
        permutation.num_rows = num_cols;
        permutation.num_cols = num_cols;
		
		 for (int i = 0; i < num_cols; ++i) {
            permutation.row_mat.push_back({i});
            permutation.col_mat.push_back({i});
        }
		
        // Create a temporary col_mat based on the permutation
        std::vector<std::vector<int>> new_col_mat(num_cols);
		std::vector<std::vector<int>> new_perm_col_mat(num_cols);
        for (int i = 0; i < num_cols; i++) {
            new_col_mat[i] = col_mat[Perm(i)];
			new_perm_col_mat[i] = permutation.col_mat[Perm(i)];
        }

        // Update col_mat with the permuted version
        col_mat = new_col_mat;
		permutation.col_mat=new_perm_col_mat;

        // Rebuild row_mat based on the updated col_mat
        col_to_row();
		permutation.col_to_row();
		return permutation;
    }
Sparse_GF2 Sparse_GF2::col_permutation_matrix(const itpp::ivec& perm) {
        int size = perm.length();
        Sparse_GF2 perm_matrix;
        perm_matrix.num_rows = size;
        perm_matrix.num_cols = size;

        for (int i = 0; i < size; i++) {
            perm_matrix.row_mat.push_back({perm(i)});
        }

        // Construct the col_mat representation
        perm_matrix.row_to_col();
        return perm_matrix;
    }
	
Sparse_GF2 Sparse_GF2::transpose() const {
        Sparse_GF2 transposed;
        transposed.num_rows = num_cols;
        transposed.num_cols = num_rows;
        transposed.row_mat = col_mat;
        transposed.col_mat = row_mat;
        return transposed;
    }
	
	void Sparse_GF2::del_row(const int&  row_idx) {
        if (row_idx < 0 || row_idx >= num_rows) {
            throw std::out_of_range("Row index out of bounds.");
        }

        // Remove the row from row_mat
        row_mat.erase(row_mat.begin() + row_idx);

        // Adjust col_mat using the row_to_col() method
        row_to_col();

        // Update the number of rows
        --num_rows;
    }

    // Delete the specified column from the matrix
    void Sparse_GF2::del_col(const int&  col_idx) {
        if (col_idx < 0 || col_idx >= num_cols) {
            throw std::out_of_range("Column index out of bounds.");
        }

        // Remove the column from col_mat
        col_mat.erase(col_mat.begin() + col_idx);

        // Adjust row_mat using the col_to_row() method
        col_to_row();

        // Update the number of columns
        --num_cols;
    }
	
	Sparse_GF2 Sparse_GF2::get_submatrix(const int&  begin_row, const int&  begin_col, const int&  end_row, const int&  end_col) {
    if (begin_row < 0 || begin_row >= num_rows || end_row < 0 || end_row >= num_rows || 
        begin_row > end_row) {
        throw std::out_of_range("Row indices out of bounds or inconsistent.");
    }

    if (begin_col < 0 || begin_col >= num_cols || end_col < 0 || end_col >= num_cols || 
        begin_col > end_col) {
        throw std::out_of_range("Column indices out of bounds or inconsistent.");
    }

    Sparse_GF2 S;
    S.num_rows = end_row - begin_row + 1;
    S.num_cols = end_col - begin_col + 1;

    // Construct the row representation of the submatrix directly
    for (int i = begin_row; i <= end_row; i++) {
        std::vector<int> sub_row;
        for (int col_idx : row_mat[i]) {
            if (col_idx >= begin_col && col_idx <= end_col) {
                sub_row.push_back(col_idx - begin_col);  // adjust the column indices
            }
        }
        S.row_mat.push_back(sub_row);
    }

    // Construct the column representation using row_to_col
    S.row_to_col();

    return S;
}

   int Sparse_GF2::col_weight(const int&  i) const {
        if (i < 0 || i >= num_cols) {
            throw std::out_of_range("Column index out of bounds.");
        }
        return col_mat[i].size();
    }

    // Return the number of non-zero entries in the ith row
    int Sparse_GF2::row_weight(const int&  i) const {
        if (i < 0 || i >= num_rows) {
            throw std::out_of_range("Row index out of bounds.");
        }
        return row_mat[i].size();
    }
 void Sparse_GF2::set(const int&  i, const int&  j, const int&  value) {
        // Boundary checks
        if (i < 0 || i >= num_rows || j < 0 || j >= num_cols) {
            throw std::out_of_range("Set:Row or column index out of bounds.");
        }
        if (value != 0 && value != 1) {
            throw std::invalid_argument("Set:Value should be either 0 or 1.");
        }

        if (value == 1) {
            // If the element doesn't already exist, add it
            if (std::find(row_mat[i].begin(), row_mat[i].end(), j) == row_mat[i].end()) {
                row_mat[i].push_back(j);
                std::sort(row_mat[i].begin(), row_mat[i].end());
            }
            if (std::find(col_mat[j].begin(), col_mat[j].end(), i) == col_mat[j].end()) {
                col_mat[j].push_back(i);
                std::sort(col_mat[j].begin(), col_mat[j].end());
            }
        } else {  // value == 0
            // If the element exists, remove it
            row_mat[i].erase(std::remove(row_mat[i].begin(), row_mat[i].end(), j), row_mat[i].end());
            col_mat[j].erase(std::remove(col_mat[j].begin(), col_mat[j].end(), i), col_mat[j].end());
        }
    }

void Sparse_OSD(const Sparse_GF2& H,vec& LR,const GF2mat& denseH,const GF2mat& syndrome,Sparse_GF2& Sparse_output_e,int& suc_order,int OSD_order) //r is the rank of H
{ 
    //get the first permutation that abs_LLR is in descending order
	
  int n=LR.length();
  if (n!=H.cols()){cout<<"OSD error: LR.length!=H.cols"<<endl;}
  vec LLR(n); //the sort function gives a ascending  order, LR=p0/p1
  
  for (int i=0;i<n;i++)
    {
		if (LR(i)>=0) {LLR(i)=log(LR(i));}
      else
		{
		cout<<"negative LR!"<<endl;
		cout<<LR(i)<<endl;
		}
    }
  


	ivec perm1=sort_index(LLR);
	Sparse_GF2 H1=H;	
	Sparse_GF2 perm1_mat=H1.permute_cols(perm1);
	  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1
	 
	// cout<<"Sparse OSD: after perm1 H is "<<H1.to_GF2mat()<<endl;
	 
	int rankH=0;
	//cout<<"before H1 is\n"<<H1.to_GF2mat()<<endl;
	Sparse_GF2 col_perm=H1.col_gaussian(rankH);
	H1=H1*col_perm;
	
	
	 //cout<<"Sparse OSD: after col_gaussian H is "<<H1.to_GF2mat()<<endl;
	Sparse_GF2 row_perm=H1.row_gaussian(rankH);
	
	H1=row_perm*H1;
	//cout<<"Sparse OSD: after row_gaussian H is "<<H1.to_GF2mat()<<endl;
	

		//cout<<"\n after H1 is\n"<<H1.to_GF2mat()<<endl;
	Sparse_GF2  sparse_s(syndrome);
	sparse_s=row_perm*sparse_s;
	
	sparse_s=sparse_s.get_submatrix(0,0,rankH-1,0);
	
    if (rankH==0){return;}
   
  Sparse_GF2 H_S=H1.get_submatrix(0,0,rankH-1,rankH-1);
  //cout<<"H_S is \n"<<H_S.to_GF2mat()<<endl;
  
  Sparse_GF2 H_T(1,1);
   
 if (H1.cols()>rankH){H_T=H1.get_submatrix(0,rankH,rankH-1,n-1);}
  Sparse_GF2 HS_inv=H_S.inverse();
 

  Sparse_GF2 e_S=HS_inv*sparse_s;  

 Sparse_GF2 e_T(1,1);
  if(n-rankH>0){ Sparse_GF2 e_T_temp(n-rankH,1); e_T=e_T_temp;}
 
  Sparse_GF2 new_e_S;
  int wt=e_S.col_weight(0);
  //cout<<"original wt "<<wt<<endl;
  int temp1;
  int temp2=-1;
  int temp3=-1;
  int  lambda=min(OSD_order,n-rankH);
  
  int temp_count=0;
  Sparse_GF2 new_e_T=e_T;
 
 if (OSD_order>=0)
 {
	 Sparse_GF2 tempmat=HS_inv*sparse_s;
	 Sparse_GF2 tempmat2=HS_inv*H_T;
   if(n-rankH>0) 
   {
	    
	   for (int i=0;i<n-rankH;i++)
		{
			 
			new_e_T.set(i,0,1);
			new_e_S=tempmat+tempmat2*new_e_T;
			temp1=new_e_S.col_weight(0);
			if (temp1+1<wt)
			{
				wt=temp1+1;
				temp2=i;
			}
			new_e_T.set(i,0,0);
		}
   }

		for (int i=0;i<lambda;i++)
		{
			for (int j=i+1;j<lambda;j++)
			{	
				new_e_T.set(i,0,1);
				new_e_T.set(j,0,1);		  
				new_e_S=tempmat+tempmat2*new_e_T;
				temp1=new_e_S.col_weight(0);

				if (temp1+2<wt)
				{
					// cout<<"another e_S"<<endl;
					wt=temp1+2;
					// cout<<"new wt is "<<wt<<endl;
					temp2=i;
					temp3=j;
				}
				new_e_T.set(i,0,0);
				new_e_T.set(j,0,0);
			}
		}

		if (temp2!=-1&&temp3==-1)
		{
			e_T.set(temp2,0,1);
			suc_order=1;
			e_S=tempmat+tempmat2*e_T;	
		}
		else if (temp2!=-1&&temp3!=-1)
		{
			e_T.set(temp2,0,1);
			e_T.set(temp2,0,1);
			suc_order=2;
			e_S=tempmat+tempmat2*e_T;
		}
   }
         vector<int> vec_e_S=e_S.get_col(0);
		  vector<int> vec_e_T=e_T.get_col(0);
	
		  
		for (int i=0;i<vec_e_S.size();i++){Sparse_output_e.set(vec_e_S[i],0,1);}
		  
		for (int i=0;i<vec_e_T.size();i++){Sparse_output_e.set(vec_e_T[i]+vec_e_S.size(),0,1);}
		Sparse_output_e=perm1_mat*col_perm*Sparse_output_e;	
		
}
/*
void initialize_checks (const Sparse_GF2mat &H, nodes checks[]){
  int r=H.rows();
  
  
  for (int i=0; i<r;i++)
    {
    checks[i].degree=H.get_row(i).size();
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
*/

bool Fcircuit_decoder(const Sparse_GF2 &sparseH,const Sparse_GF2 &sparseL,const GF2mat &Hx, const GF2mat& Gx, const GF2mat& G,const GF2mat &real_e,const nodes checks[],const nodes errors[],const vec&pv_dec,double& num_iter, const int& lmax,const int &wt,  const int& debug, const int& OSD_order,Sparse_GF2& Sparse_zero_rvec, GF2mat& zero_mat1)
{
 int c=Hx.rows();
 int v=Hx.cols();
 
 //if no error, break

  Sparse_GF2 Sparse_real_e(real_e);
  Sparse_GF2 Sparse_syndrome=sparseH*Sparse_real_e;
  GF2mat syndrome=Sparse_syndrome.to_GF2mat();
	
  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (sparseL*Sparse_real_e==Sparse_zero_rvec)
	    {
			num_no_error++;
			return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	
	//this part can be optimized:
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
      initialize_massages( mcv,mvc, Hx); //initialize to all-1 matrix	
	  
	  //GF2mat output_e(v,1);	
	  vec LR(v);	
	  bool schedule=debug&2;
      for (int l=1;l<=lmax;l++)
	{		 
			Sparse_GF2 Sparse_output_e(v,1);
			if (schedule)
			{
				quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,Sparse_output_e,LR,debug);				
			}
	      else
			{		
				quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,Sparse_output_e,LR,debug);
			}
	    
			//Sparse_GF2 Sparse_output_e(output_e);
			
			if (sparseH*Sparse_output_e==Sparse_syndrome)		   
			{					
				  if(sparseL*(Sparse_output_e+Sparse_real_e)==Sparse_zero_rvec)
					{
						num_iter= num_iter+l;		
						BP_suc++;
						return true;
					}	    
				else{return false; }	    	  
		}
	}
      						
      if (debug&8)
	{
		int suc_order=0;
			  Sparse_GF2 Sparse_output_eOSD(v,1);
			  
	  Sparse_OSD(sparseH,LR,Hx,syndrome,Sparse_output_eOSD,suc_order,OSD_order);
	  /*
	  GF2mat output_e(v,1);
			  GF2mat Hright;
			  GF2mat Hwrong;
	   OSD(LR, Hx,syndrome,output_e,suc_order,Hright);
		if (Hwrong==Hright){}
		else 
		{
			cout<<"error: Hright is\n"<<Hright<<"\n Hwrongis \n"<<Hwrong<<endl;
		}
   if (output_e==Sparse_output_eOSD.to_GF2mat()){}
   else {
	   cout<<"error : output_e is \n"<<output_e.transpose()<<"\n wrong one:\n"<<(Sparse_output_eOSD.to_GF2mat()).transpose()<<endl;
   }
   */
   
	  if(sparseL*(Sparse_output_eOSD+Sparse_real_e)==Sparse_zero_rvec)
	    {
	      OSD_suc++;
		  if (suc_order==0){OSD_0_suc++;}
		else if(suc_order==1){OSD_1_suc++;}
		else if(suc_order==2){OSD_2_suc++;}
		else {OSD_higher_order_suc++;}
	      return true;
	    }	 
		else 
		{		
			OSD_fail++;
			syn_fail++;			
			return false;
		}
	}
      max_fail++;
      return false;
 }
 
 void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int& c, int& v,  Sparse_GF2& output_e, vec &LR,const int& debug,double alpha)
{  
    double ipr;
   for (int j=0;j<v;j++)
    {
		int vj_degree=errors[j].degree;
		ipr=(1-pv[j])/pv[j];
		LR(j)=ipr;
		for (int i=0;i<vj_degree;i++)
		{     
			int ci=(errors[j].neighbors)(i);
			update_ci_to_vj(checks, errors,mcv, mvc,ci,j,syndrome(ci,0));
			LR(j)=LR(j)*mcv(ci,j);      	
	   } 
		//ci is the ith neibor of vj:
		for (int i=0;i<vj_degree;i++)
		{
			//  update the  v_j to c_i massage:      
			int ci=(errors[j].neighbors)(i);
			update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);     		
		}
		if (LR(j)<1) {output_e.set(j,0,1);}
     }  
}

void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  Sparse_GF2& output_e, vec &LR,int debug,double alpha)
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
   if (final_pr<1) {output_e.set(j,0,1);}  
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
