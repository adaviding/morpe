#include "morpe.h"

#include <algorithm>
#include <ranges>

using namespace morpe::numerics;

namespace morpe
{
    //---------------------------------------
    // constructors
    //---------------------------------------

    polynomial::polynomial(int32_t num_dims, int32_t rank)
    {
        ThrowIf(num_dims < 1);
        ThrowIf(rank < 1)

        this->num_dims = num_dims;
        this->rank = rank;

        // Allocate the vector of polynomial coefficients.
        int32_t n = num_coeffs(num_dims, rank);
        this->coeffs.resize(n);

        // Calculate the number of polynomial coefficients for all ranks up to the given rank, and for the
        // homogeneous and inhomogeneous cases.
        this->num_coeffs_each_rank.resize(rank);      // inhomogeneous
        this->num_coeffs_each_rank_homo.resize(rank); // homogeneous

        int32_t num_prior = 0;
        int32_t num_current = 0;
        for (int32_t i_rank = 0; i_rank < rank; i_rank++)
        {
            num_prior = num_current;
            num_current = num_coeffs(num_dims, i_rank+1);
            this->num_coeffs_each_rank[i_rank] = num_current;
            this->num_coeffs_each_rank_homo[i_rank] = num_current - num_prior;
        }

        std::vector<int32_t> row = { -1 };

        // For each coefficient
        for (int i_coeff = 0; i_coeff < n; i_coeff++)
        {
            int i,j;

            // Increment the smallest digit and begin to handle "carry-over" arithmetic.
            for (i = 0; i < row.size(); i++)
            {
                int val = ++row[i];
                if (val < num_dims)
                {
                    break;
                }
                else if (i == row.size() - 1)
                {
                    i++;
                    row.push_back(0);
                }
            }

            // Finish the "carry-over" by ensuring that any leftward maxed digits have been reset.
            for( j=i-1; j>=0; j-- )
            {
                if( row[j]==num_dims )
                    row[j]=row[j+1];
            }

            // Our final answer is just the reverse of the row.
            this->coeffs[i_coeff] = row;
            std::reverse(this->coeffs[i_coeff].begin(), this->coeffs[i_coeff].end());
        }
    }

    //---------------------------------------
    // public member functions
    //---------------------------------------

    Eigen::ArrayXf polynomial::expand(
            _In_ const Eigen::VectorXf &input)
    {
        Eigen::VectorXf output(this->coeffs.size());
        output.resize(this->coeffs.size());
        this->expand(input, output);
        return output;
    }


    void polynomial::expand(
            _In_    const Eigen::VectorXf& input,
            _Inout_ Eigen::VectorXf& output)
    {
        ChkTrue(this->num_dims == input.size());
        ChkTrue(this->coeffs.size() == output.size());

        double y;
        for (int i_coeff = 0; i_coeff < this->coeffs.size(); i_coeff++)
        {
            y = 1.0;

            std::vector<int32_t>* coeff = &this->coeffs[i_coeff];

            for (int i_term = 0; i_term < coeff->size(); i_term++)
            {
                int32_t j = (*coeff)[i_term];
                y *= input[j];
            }

            output[i_coeff] = (float)y;
        }
    }

    //---------------------------------------
    // public static functions
    //---------------------------------------

    std::vector<int32_t> polynomial::map_subspace_to_fullspace(
            _In_ const polynomial& fullpoly,
            _In_ const polynomial& subpoly,
            _In_ const std::vector<int32_t>& subdims)
    {
        ChkTrue(subpoly.num_dims < fullpoly.num_dims)
        ChkTrue(subpoly.rank <= fullpoly.rank);
        ChkTrue(subdims.size() == subpoly.num_dims);
        ChkTrueMsg(subdims.end() == std::ranges::adjacent_find(subdims, std::ranges::greater_equal()),
                   "The list of subspace dimensions must be non-decreasing.");
        ChkTrue(0 <= subdims.front());
        ChkTrue(subdims.back() < fullpoly.num_dims);

        std::vector<int32_t> output(subpoly.coeffs.size());
        std::vector<int32_t> mapping(subpoly.rank);

        int j_coeff = 0;
        for (int i_coeff = 0; i_coeff < output.size(); i_coeff++)
        {
            const std::vector<int32_t>* sub = &subpoly.coeffs[i_coeff];
            mapping.resize(sub->size());

            for (int i_rank = 0; i_rank < sub->size(); i_rank++)
            {
                mapping[i_rank] = subdims[sub->at(i_rank)];
            }

            while (!std::ranges::equal(fullpoly.coeffs[j_coeff], mapping))
            {
                j_coeff++;
            }
            output[i_coeff] = j_coeff;
            j_coeff++;
        }

        return output;
    }

    int32_t polynomial::num_coeffs(int32_t num_dims, int32_t rank)
    {
        // This is the number of coefficients including the 0-th order term.
        int32_t output = I2::pascal(num_dims, rank);

        //                 0    1    2    3     4     5     6     7     8     9
        //                            Rank of Tensor
        //    0            1    1    1    1     1     1     1     1     1     1
        //    1            1    2    3    4     5     6     7     8     9     10
        //    2    Ndims   1    3    6    10    15    21    28    36    45    55
        //    3    of      1    4    10   20    35    56    84    120   165   220
        //    4    Space   1    5    15   35    70    126   210   330   495   715
        //    5            1    6    21   56    126   252   462   924   1716

        // We subtract 1 to remove the 0-th order term.  (MoRPE does not use this term.)
        output--;

        return output;
    }

    int32_t polynomial::num_coeffs_homo(int32_t num_dims, int32_t rank)
    {
        int32_t output = I2::pascal(num_dims - 1, rank);
        return output;
    }
}
