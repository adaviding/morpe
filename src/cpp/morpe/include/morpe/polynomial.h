#pragma once

#include "_internal.h"
#include <vector>

#include "Eigen/Eigen"

namespace morpe
{
    /// Represents the coefficients of an inhomogenous polynomial.
    class polynomial
    {
    public:
        //-------------------
        // Member fields
        //-------------------

        /// The inhomogeneous polynomial coefficients.  Each row is non-decreasing and defines a coefficient.
        ///
        ///    For example...
        ///        A row of {0, 0, 2, 5} would correspond to the following polynomial term with a rank of 4.
        ///            x[0] * x[0] * x[2] * x[5]
        ///        A row of {0, 1, 2} would correspond to the following polynomial term with a rank of 3.
        ///            x[0] * x[1] * x[2]
        ///
        /// The number of values in each row is equal to the rank of that coefficient.  The maximum rank is #rank.
        std::vector<std::vector<int32_t>> coeffs;

        /// The number of coefficients for each rank, from 1 to #rank, for the inhomogeneous polynomial.
        ///
        /// Since the array index is zero based, ...
        /// * A rank of 1 corresponds to an index of 0
        /// * A rank of 2 corresponds to an index of 1
        ///   ...
        /// * A rank of R corresponds to an index of R-1
        std::vector<int32_t> num_coeffs_each_rank;

        /// This is like #num_coeffs_each_rank, but for the homogeneous polynomial.
        std::vector<int32_t> num_coeffs_each_rank_homo;

        /// The number of spatial dimensions over which the polynomial is defined.  (i.e. The number of coordinate axes.)
        int32_t num_dims;

        /// The maximum rank of the polynomial coefficients.
        ///
        /// For example:
        ///   * If rank = 1, then the polynomial is linear.
        ///   * If rank = 2, then the polynomial is quadratic.
        ///   * If rank = 3, then the polynomial is cubic.
        ///   * ... and so on
        int32_t rank;

        //-------------------
        // Member functions
        //-------------------

        polynomial(int32_t num_dims, int32_t rank);

        Eigen::ArrayXf expand(
                _In_ Eigen::VectorXf& input);

        void expand(
                _In_    Eigen::VectorXf& input,
                _Inout_ Eigen::VectorXf& output);

        //-------------------
        // Static functions
        //-------------------

        static std::vector<int32_t> map_subspace_to_fullspace(
                _In_ polynomial& fullPoly,
                _In_ polynomial& subPoly,
                _In_ std::vector<int32_t> subdims);

        static int32_t num_coeffs(int32_t num_dims, int32_t rank);

        static int32_t num_coeffs_homo(int32_t num_dims, int32_t rank);
    };
}
