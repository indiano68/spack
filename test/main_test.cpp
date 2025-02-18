#include<linalg.h>
#include<iostream>



int main()
{
    dense_matrix new_matrix = dense_matrix(10,10);
    for(size_t i = 0 ; i < new_matrix.nrows(); i++)
    for(size_t j = 0 ; j < new_matrix.ncols(); j++)
    {
            new_matrix(i,j) = 10*i+j;
    }
    for(size_t idx = 0 ; idx < new_matrix.size(); idx++)
    {
        std::cout << new_matrix[idx] << " " ;
        if((idx + 1) % new_matrix.ncols() == 0) std::cout <<  std::endl;
        
    }
}
