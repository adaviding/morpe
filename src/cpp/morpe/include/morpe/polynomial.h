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
        //---------------------------------------
        // Member fields
        //---------------------------------------

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

        //---------------------------------------
        // Member functions
        //---------------------------------------

        /// Constructs the coefficients for a multivariate polynomial having the given rank and spatial dimensionality.
        /// @param num_dims The number of spatial dimensions, see #num_dims.
        /// @param rank The rank of the polynomial, see #rank.
        polynomial(int32_t num_dims, int32_t rank);

        /// Computes 'output' as the polynomial expansion of 'input'.
        /// @input  The feature vector (not expanded).
        /// @return The expanded feature vector.
        Eigen::ArrayXf expand(
                _In_ const Eigen::VectorXf& input);

        /// Computes 'output' as the polynomial expansion of 'input'.
        /// @input  The feature vector (not expanded).
        /// @output The expanded feature vector, pre-allocated to be the correct length.
        void expand(
                _In_    const Eigen::VectorXf& input,
                _Inout_ Eigen::VectorXf& output);

        //---------------------------------------
        // Static functions
        //---------------------------------------

        /// Computes a mapping from the subspace polynomial coefficients to the fullspace coefficients.
        /// @param fullpoly The fullspace polynomial
        /// @param subpoly  The subspace polynomial
        /// @param subdims  For each spatial dimension of the subspace polynomial, this gives the index of the
        /// corresponding spatial dimension in the fullspace polynomial.
        /// @return For each coefficient of the subspace polynomial, this gives the corresponding index of the coefficient
        /// in the fullspace polynomial.
        static std::vector<int32_t> map_subspace_to_fullspace(
                _In_ const polynomial& fullpoly,
                _In_ const polynomial& subpoly,
                _In_ const std::vector<int32_t>& subdims);

        /// Returns the number of inhomogeneous polynomial coefficients for a given dimensionality and rank, including
        /// all coefficients of lesser rank.  However, it does not include the 0-th order term which MoRPE does not use.
        /// @param num_dims #num_dims  The spatial dimensionality.
        /// @param rank     #rank      The rank.
        static int32_t num_coeffs(int32_t num_dims, int32_t rank);

        /// Returns the number of homogeneous polynomial coefficients for a given dimensionality and rank.
        ///
        /// This does NOT including coefficients of lesser rank (as per the term "homogeneous").
        /// @param num_dims #num_dims  The spatial dimensionality.
        /// @param rank     #rank      The rank.
        static int32_t num_coeffs_homo(int32_t num_dims, int32_t rank);
    };
}
