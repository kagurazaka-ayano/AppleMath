/**
 * @file MatrixTest.cpp
 * @author ayano
 * @date 2/14/24
 * @brief
*/

#include <gtest/gtest.h>
#include "Math/Matrix.hpp"
#include "Math/configuration.hpp"

TEST(MatrixTests, MatrixConstructorWithArrayProducesCorrectResult) {
    std::array<std::array<double, 2>, 2> data = {{ {1.0, 2.0}, {3.0, 4.0} }};
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

TEST(MatrixTests, MakeRotationMatrixR2ThrowsForInvalidAngle) {
    EXPECT_THROW(AppleMath::makeRotationMatrixR2(M_PI * 2.1), std::invalid_argument);
}

TEST(MatrixTests, MakeRotationMatrixR3ThrowsForInvalidAngles) {
    EXPECT_THROW(AppleMath::makeRotationMatrixR3(M_PI * 2.1, M_PI / 2, M_PI / 2), std::invalid_argument);
    EXPECT_THROW(AppleMath::makeRotationMatrixR3(M_PI / 2, M_PI * 2.1, M_PI / 2), std::invalid_argument);
    EXPECT_THROW(AppleMath::makeRotationMatrixR3(M_PI / 2, M_PI / 2, M_PI * 2.1), std::invalid_argument);
}

TEST(MatrixTests, MakeRotationMatrixR2ProducesCorrectResult) {
    AppleMath::Matrix<2, 2> matrix = AppleMath::makeRotationMatrixR2(M_PI / 2);
    EXPECT_NEAR(matrix(0, 0), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(matrix(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 1.0);
    EXPECT_NEAR(matrix(1, 1), 0.0, AppleMath::EPS);
}

TEST(MatrixTests, MakeRotationMatrixR3ProducesCorrectResult) {
    auto matrix = AppleMath::makeRotationMatrixR3(M_PI / 2, M_PI / 2, M_PI / 2);
    auto x_matrix = matrix.at("x");
    auto y_matrix = matrix.at("y");
    auto z_matrix = matrix.at("z");
    EXPECT_DOUBLE_EQ(x_matrix(0, 0), 1.0);
    EXPECT_NEAR(x_matrix(0, 1), 0.0, AppleMath::EPS);
    EXPECT_NEAR(x_matrix(0, 2), 0.0, AppleMath::EPS);
    EXPECT_NEAR(x_matrix(1, 0), 0.0, AppleMath::EPS);
    EXPECT_NEAR(x_matrix(1, 1), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(x_matrix(1, 2), -1.0);
    EXPECT_NEAR(x_matrix(2, 0), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(x_matrix(2, 1), 1.0);
    EXPECT_NEAR(x_matrix(2, 2), 0.0, AppleMath::EPS);

    EXPECT_NEAR(y_matrix(0, 0), 0.0, AppleMath::EPS);
    EXPECT_NEAR(y_matrix(0, 1), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(y_matrix(0, 2), 1.0);
    EXPECT_NEAR(y_matrix(1, 0), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(y_matrix(1, 1), 1.0);
    EXPECT_NEAR(y_matrix(1, 2), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(y_matrix(2, 0), -1.0);
    EXPECT_NEAR(y_matrix(2, 1), 0.0, AppleMath::EPS);
    EXPECT_NEAR(y_matrix(2, 2), 0.0, AppleMath::EPS);

    EXPECT_NEAR(z_matrix(0, 0), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(z_matrix(0, 1), -1.0);
    EXPECT_NEAR(z_matrix(0, 2), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(z_matrix(1, 0), 1.0);
    EXPECT_NEAR(z_matrix(1, 1), 0.0, AppleMath::EPS);
    EXPECT_NEAR(z_matrix(1, 2), 0.0, AppleMath::EPS);
    EXPECT_NEAR(z_matrix(2, 0), 0.0, AppleMath::EPS);
    EXPECT_NEAR(z_matrix(2, 1), 0.0, AppleMath::EPS);
    EXPECT_DOUBLE_EQ(z_matrix(2, 2), 1.0);
}

TEST(MatrixTests, MatrixVectorMultiplication) {
    AppleMath::Matrix<2, 2> matrix({ {1.0, 2.0}, {3.0, 4.0} });
    AppleMath::Vector<2> vector{1.0, 2.0};
    auto result = matrix * vector;
    EXPECT_DOUBLE_EQ(result[0], 7.0);
    EXPECT_DOUBLE_EQ(result[1], 10.0);
}
