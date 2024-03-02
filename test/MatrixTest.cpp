/**
 * @file MatrixTests.cpp
 * @author ayano
 * @date 2/14/24
 * @brief
*/

#include <gtest/gtest.h>
#include "AppleMath/Matrix.hpp"
#include "AppleMath/Vector.hpp"
#include "AppleMath/configuration.hpp"

TEST(MatrixTests, MatrixConstructorWithArrayProducesCorrectResult) {
    std::array<std::array<float, 2>, 2> data = {{ {1.0, 2.0}, {3.0, 4.0} }};
    AppleMath::Matrix<2, 2> matrix(data);
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
}

TEST(MatrixTests, MatrixConstructorWithInitializerListProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix({ {1.0, 2.0}, {3.0, 4.0} });
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
}

TEST(MatrixTests, MatrixAdditionAssignmentProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix1({ {1.0, 2.0}, {3.0, 4.0} });
    AppleMath::Matrix<2, 2> matrix2({ {5.0, 6.0}, {7.0, 8.0} });
    matrix1 += matrix2;
    EXPECT_DOUBLE_EQ(matrix1(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(matrix1(0, 1), 10.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 0), 8.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 1), 12.0);
}

TEST(MatrixTests, MatrixSubtractionAssignmentProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix1({ {1.0, 2.0}, {3.0, 4.0} });
    AppleMath::Matrix<2, 2> matrix2({ {5.0, 6.0}, {7.0, 8.0} });
    matrix1 -= matrix2;
    EXPECT_DOUBLE_EQ(matrix1(0, 0), -4.0);
    EXPECT_DOUBLE_EQ(matrix1(0, 1), -4.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 0), -4.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 1), -4.0);
}

TEST(MatrixTests, MatrixMultiplicationAssignmentProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix({ {1.0, 2.0}, {3.0, 4.0} });
    matrix *= 2.0;
    EXPECT_DOUBLE_EQ(matrix(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 6.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 8.0);
}

TEST(MatrixTests, MatrixDivisionAssignmentProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix({ {2.0, 4.0}, {6.0, 8.0} });
    matrix /= 2.0;
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
}

TEST(MatrixTests, MatrixMultiplicationAssignmentWithMatrixProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix1({ {1.0, 2.0}, {3.0, 4.0} });
    AppleMath::Matrix<2, 2> matrix2({ {5.0, 6.0}, {7.0, 8.0} });
    matrix1 *= matrix2;
    EXPECT_DOUBLE_EQ(matrix1(0, 0), 23.0);
    EXPECT_DOUBLE_EQ(matrix1(0, 1), 31.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 0), 34.0);
    EXPECT_DOUBLE_EQ(matrix1(1, 1), 46.0);
}

TEST(MatrixTests, makeRotationMatrixR2ProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix = AppleMath::makeRotationMatrixR2(M_PI / 2);
    EXPECT_NEAR(matrix(0, 0), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(matrix(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 1.0);
    EXPECT_NEAR(matrix(1, 1), 0.0, AppleMath::EPS);
}

TEST(MatrixTest, TestMakeEulerRotationMatrixR3) {
    float psi = M_PI / 2; // 90 degrees
    float theta = M_PI; // 180 degrees
    float phi = M_PI / 2; // 90 degrees

    auto result = AppleMath::makeEulerRotationMatrixR3(psi, theta, phi);

    // Expected result for these angles
    AppleMath::Matrix<3, 3> expected {
        {{0, -1, 0},
        {0, 0, -1},
        {1, 0, 0}}
    };
    ASSERT_NEAR((result[0, 0]), (expected[0, 0]), (1e-3));
    ASSERT_NEAR((result[0, 1]), (expected[0, 1]), (1e-3));
    ASSERT_NEAR((result[0, 2]), (expected[0, 2]), (1e-3));
    ASSERT_NEAR((result[1, 0]), (expected[1, 0]), (1e-3));
    ASSERT_NEAR((result[1, 1]), (expected[1, 1]), (1e-3));
    ASSERT_NEAR((result[1, 2]), (expected[1, 2]), (1e-3));
    ASSERT_NEAR((result[2, 0]), (expected[2, 0]), (1e-3));
    ASSERT_NEAR((result[2, 1]), (expected[2, 1]), (1e-3));
    ASSERT_NEAR((result[2, 2]), (expected[2, 2]), (1e-3));
}

TEST(MatrixTests, MatrixVectorMultiplication) {
    AppleMath::Matrix<2, 2> matrix({ {1.0, 2.0}, {3.0, 4.0} });
    AppleMath::Vector<2> vector{1.0, 2.0};
    auto result = matrix * vector;
    EXPECT_DOUBLE_EQ(result[0], 7.0);
    EXPECT_DOUBLE_EQ(result[1], 10.0);
}


TEST(MatrixTests, IdentityMatrix2x2) {
    auto result = AppleMath::makeIdentity<2>();
    auto cmp = AppleMath::Matrix<2, 2>(matrix_identity_float2x2);
    EXPECT_EQ(result, cmp);
}

TEST(MatrixTests, IdentityMatrix3x3) {
    auto result = AppleMath::makeIdentity<3>();
    auto cmp = AppleMath::Matrix<3, 3>(matrix_identity_float3x3);
    EXPECT_EQ(result, cmp);
}

TEST(MatrixTests, IdentityMatrix4x4) {
    auto result = AppleMath::makeIdentity<4>();
    auto cmp = AppleMath::Matrix<4, 4>(matrix_identity_float4x4);
    EXPECT_EQ(result, cmp);
}

TEST(MatrixTests, RotationMatrixR2) {
    double angle = M_PI / 4.0; // 45 degrees
    auto result = AppleMath::makeRotationMatrixR2(angle);
    EXPECT_NEAR(result(0, 0), std::cos(angle), 1e-3);
    EXPECT_NEAR(result(0, 1), -std::sin(angle), 1e-3);
    EXPECT_NEAR(result(1, 0), std::sin(angle), 1e-3);
    EXPECT_NEAR(result(1, 1), std::cos(angle), 1e-3);
}

TEST(MatrixTests, RotationMatrixR3) {
    double angle_phi = M_PI / 4.0; // 45 degrees
    double angle_theta = M_PI / 3.0; // 60 degrees
    double angle_psi = M_PI / 2.0; // 90 degrees
    auto result = AppleMath::makeEulerRotationMatrixR3(angle_psi, angle_theta, angle_phi);
    // Expected result for these angles
    AppleMath::Matrix<3, 3> expected {
        {{0.354, 0.354, -0.866},
        {0.612, 0.612, 0.5},
        {0.707, -0.707, 0}}
    };
    ASSERT_NEAR((result[0, 0]), (expected[0, 0]), (1e-3));
    ASSERT_NEAR((result[0, 1]), (expected[0, 1]), (1e-3));
    ASSERT_NEAR((result[0, 2]), (expected[0, 2]), (1e-3));
    ASSERT_NEAR((result[1, 0]), (expected[1, 0]), (1e-3));
    ASSERT_NEAR((result[1, 1]), (expected[1, 1]), (1e-3));
    ASSERT_NEAR((result[1, 2]), (expected[1, 2]), (1e-3));
    ASSERT_NEAR((result[2, 0]), (expected[2, 0]), (1e-3));
    ASSERT_NEAR((result[2, 1]), (expected[2, 1]), (1e-3));
    ASSERT_NEAR((result[2, 2]), (expected[2, 2]), (1e-3));
}

TEST(MatrixTests, ScaleMatrix2x2) {
    double x = 2.0;
    double y = 3.0;
    auto result = AppleMath::makeScaleMatrix<2>(x, y);
    EXPECT_EQ(result(0, 0), x);
    EXPECT_EQ(result(1, 1), y);
}

TEST(MatrixTests, ScaleMatrix3x3) {
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    auto result = AppleMath::makeScaleMatrix<3>(x, y, z);
    EXPECT_EQ(result(0, 0), x);
    EXPECT_EQ(result(1, 1), y);
    EXPECT_EQ(result(2, 2), z);
}

TEST(MatrixTests, TranslationMatrixR2) {
    double x = 2.0;
    double y = 3.0;
    auto result = AppleMath::makeTranslationMatrixR2(x, y);
    EXPECT_EQ(result(0, 2), x);
    EXPECT_EQ(result(1, 2), y);
}

TEST(MatrixTests, TranslationMatrixR3) {
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    auto result = AppleMath::makeTranslationMatrixR3(x, y, z);
    EXPECT_EQ(result(0, 3), x);
    EXPECT_EQ(result(1, 3), y);
    EXPECT_EQ(result(2, 3), z);
}
