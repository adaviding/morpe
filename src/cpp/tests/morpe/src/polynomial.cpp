#include "test_morpe.h"

#include <unordered_set>

TEST(polynomial_tests, num_coeffs)
{
    EXPECT_EQ( 1, morpe::polynomial::num_coeffs(1, 1));
    EXPECT_EQ( 2, morpe::polynomial::num_coeffs(2, 1));
    EXPECT_EQ( 5, morpe::polynomial::num_coeffs(2, 2));
    EXPECT_EQ( 3, morpe::polynomial::num_coeffs(3, 1));
    EXPECT_EQ( 9, morpe::polynomial::num_coeffs(3, 2));
    EXPECT_EQ(19, morpe::polynomial::num_coeffs(3, 3));
    EXPECT_EQ(34, morpe::polynomial::num_coeffs(3, 4));
    EXPECT_EQ( 4, morpe::polynomial::num_coeffs(4, 1));
    EXPECT_EQ(14, morpe::polynomial::num_coeffs(4, 2));

    EXPECT_EQ( 1, morpe::polynomial::num_coeffs_homo(1, 1));
    EXPECT_EQ( 2, morpe::polynomial::num_coeffs_homo(2, 1));
    EXPECT_EQ( 3, morpe::polynomial::num_coeffs_homo(2, 2));
    EXPECT_EQ( 3, morpe::polynomial::num_coeffs_homo(3, 1));
    EXPECT_EQ( 6, morpe::polynomial::num_coeffs_homo(3, 2));
    EXPECT_EQ(10, morpe::polynomial::num_coeffs_homo(3, 3));
    EXPECT_EQ(15, morpe::polynomial::num_coeffs_homo(3, 4));
    EXPECT_EQ( 4, morpe::polynomial::num_coeffs_homo(4, 1));
    EXPECT_EQ(10, morpe::polynomial::num_coeffs_homo(4, 2));
}

void coeffs_impl(int32_t num_dims, int32_t rank)
{
    morpe::polynomial poly(num_dims, rank);

    EXPECT_EQ(num_dims, poly.num_dims);
    EXPECT_EQ(rank, poly.rank);

    // We should have the correct number of coefficients.
    EXPECT_EQ(poly.coeffs.size(), morpe::polynomial::num_coeffs(num_dims, rank));

    // The number of inhomogeneous coefficients should increase as a function of rank
    EXPECT_TRUE(poly.num_coeffs_each_rank.end() == std::ranges::adjacent_find(poly.num_coeffs_each_rank, std::ranges::greater()));

    // The number of homogeneous coefficients should increase as a function of rank
    EXPECT_TRUE(poly.num_coeffs_each_rank_homo.end() == std::ranges::adjacent_find(poly.num_coeffs_each_rank_homo, std::ranges::greater()));

    // We are going to build a unique string to represent each coefficient, to ensure there are no duplicates.
    std::unordered_set<std::string> unique_coeff_strings;

    size_t len=0, len_prior=0;

    for (size_t i = 0; i < poly.coeffs.size(); i++)
    {
        // Here we build the unique string for the coefficient.
        std::stringstream unique_coeff_string;

        // The rank of each coefficient is its size().  The rank should make sense.
        len = poly.coeffs[i].size();
        EXPECT_TRUE(len <= rank);
        EXPECT_TRUE(len >= len_prior);
        len_prior = len;

        // The prior value of 'current' (defined below).
        int32_t prior = -1;  // dummy value

        for (int32_t j = 0; j < len; j++)
        {
            int32_t current = poly.coeffs[i][j];

            // All terms must be non-negative and less than the spatial dimensionality.
            EXPECT_TRUE(current >= 0);
            EXPECT_TRUE(current < num_dims);

            if (i > 0)
            {
                EXPECT_LE(prior, current);  // The terms should be non-decreasing.
                unique_coeff_string << "_";
            }

            unique_coeff_string << current;

            prior = current;
        }

        // Ensure that each coefficient is unique.
        EXPECT_TRUE(unique_coeff_strings.emplace(unique_coeff_string.str()).second);
    }
}

TEST(polynomial_tests, coeffs)
{
    coeffs_impl(1,1);
    coeffs_impl(2,1);
    coeffs_impl(3,1);
    coeffs_impl(4,1);
    coeffs_impl(2,2);
    coeffs_impl(3,2);
    coeffs_impl(4,2);
    coeffs_impl(2,3);
    coeffs_impl(3,3);
    coeffs_impl(4,3);
    coeffs_impl(1,4);
    coeffs_impl(2,4);
    coeffs_impl(4,4);
}
