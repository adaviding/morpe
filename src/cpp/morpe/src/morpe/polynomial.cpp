#include "morpe.h"

#include <algorithm>
#include <ranges>

using namespace morpe::numerics;

namespace morpe
{
    /// Constructs the coefficients for a multivariate polynomial having the given rank and spatial dimensionality.
    /// @param num_dims The number of spatial dimensions, see #num_dims.
    /// @param rank The rank of the polynomial, see #rank.
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

    /// Computes 'output' as the polynomial expansion of 'input'.
    /// @input  The feature vector (not expanded).
    /// @return The expanded feature vector.
    Eigen::ArrayXf polynomial::expand(
            _In_ Eigen::VectorXf &input)
    {
        Eigen::VectorXf output(this->coeffs.size());
        output.resize(this->coeffs.size());
        this->expand(input, output);
        return output;
    }


    /// Computes 'output' as the polynomial expansion of 'input'.
    /// @input  The feature vector (not expanded).
    /// @output The expanded feature vector, pre-allocated to be the correct length.
    void polynomial::expand(
            _In_    Eigen::VectorXf& input,
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

    /// STATIC
    /// Computes a mapping from the subspace polynomial coefficients to the fullspace coefficients.
    /// @param fullpoly The fullspace polynomial
    /// @param subpoly  The subspace polynomial
    /// @param subdims  For each spatial dimension of the subspace polynomial, this gives the index of the
    /// corresponding spatial dimension in the fullspace polynomial.
    /// @return For each coefficient of the subspace polynomial, this gives the corresponding index of the coefficient
    /// in the fullspace polynomial.
    std::vector<int32_t> polynomial::map_subspace_to_fullspace(
            _In_ polynomial &fullpoly,
            _In_ polynomial &subpoly,
            _In_ std::vector<int32_t> subdims)
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
            std::vector<int32_t>* sub = &subpoly.coeffs[i_coeff];
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

    /// STATIC
    /// Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank, including
    /// all coefficients of lesser rank.  However, it does not include the 0-th order term which MoRPE does not use.
    /// @param num_dims #num_dims  The spatial dimensionality.
    /// @param rank     #rank      The rank.
    int32_t polynomial::num_coeffs(int32_t num_dims, int32_t rank)
    {
        // This is the number of coefficients including the 0-th order term.
        int32_t output = IX::pascal(num_dims, rank);

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

    /// STATIC
    /// Returns the number of homogeneous polynomial coefficients for a given dimensionality and rank.
    ///
    /// This does NOT including coefficients of lesser rank (as per the term "homogeneous").
    /// @param num_dims #num_dims  The spatial dimensionality.
    /// @param rank     #rank      The rank.
    int32_t polynomial::num_coeffs_homo(int32_t num_dims, int32_t rank)
    {
        int32_t output = IX::pascal(num_dims - 1, rank);
        return output;
    }
}
