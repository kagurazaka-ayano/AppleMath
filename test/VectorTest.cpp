#include <gtest/gtest.h>
#include "Math/Vector.h"

using namespace Math;

TEST(Vector2Tests, InitializationWithTwoParameters) {
    Vector2 v(1.0, 2.0);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
}

TEST(Vector2Tests, CopyConstructor) {
    Vector2 v1(1.0, 2.0);
    const Vector2& v2(v1);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
}

TEST(Vector2Tests, AdditionOperator) {
    Vector2 v1(1.0, 2.0);
    Vector2 v2(3.0, 4.0);
    Vector2 v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3[0], 4.0);
    EXPECT_DOUBLE_EQ(v3[1], 6.0);
}

TEST(Vector2Tests, SubtractionOperator) {
    Vector2 v1(1.0, 2.0);
    Vector2 v2(3.0, 4.0);
    Vector2 v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3[0], -2.0);
    EXPECT_DOUBLE_EQ(v3[1], -2.0);
}

TEST(Vector2Tests, MultiplicationOperator) {
    Vector2 v1(1.0, 2.0);
    Vector2 v2 = v1 * 3.0;
    EXPECT_DOUBLE_EQ(v2[0], 3.0);
    EXPECT_DOUBLE_EQ(v2[1], 6.0);
}

TEST(Vector2Tests, DivisionOperator) {
    Vector2 v1(1.0, 2.0);
    Vector2 v2 = v1 / 2.0;
    EXPECT_DOUBLE_EQ(v2[0], 0.5);
    EXPECT_DOUBLE_EQ(v2[1], 1.0);
}

TEST(Vector2Tests, Length) {
    Vector2 v(3.0, 4.0);
    EXPECT_DOUBLE_EQ(v.length(), 5.0);
}

TEST(Vector2Tests, Normalization) {
    Vector2 v(3.0, 4.0);
    Vector2 normalized = v.normalized();
    EXPECT_NEAR(normalized.length(), 1.0, 1e-4);
}

TEST(Vector2Tests, DotProduct) {
    Vector2 v1(1.0, 2.0);
    Vector2 v2(3.0, 4.0);
    double dotProduct = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dotProduct, 11.0);
}

TEST(Vector3Tests, InitializationWithThreeParameters) {
    Vector3 v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}

TEST(Vector3Tests, CopyConstructor) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(v1);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

TEST(Vector3Tests, AdditionOperator) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    Vector3 v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3[0], 5.0);
    EXPECT_DOUBLE_EQ(v3[1], 7.0);
    EXPECT_DOUBLE_EQ(v3[2], 9.0);
}

TEST(Vector3Tests, SubtractionOperator) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    Vector3 v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3[0], -3.0);
    EXPECT_DOUBLE_EQ(v3[1], -3.0);
    EXPECT_DOUBLE_EQ(v3[2], -3.0);
}

TEST(Vector3Tests, MultiplicationOperator) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2 = v1 * 3.0;
    EXPECT_DOUBLE_EQ(v2[0], 3.0);
    EXPECT_DOUBLE_EQ(v2[1], 6.0);
    EXPECT_DOUBLE_EQ(v2[2], 9.0);
}

TEST(Vector3Tests, DivisionOperator) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2 = v1 / 2.0;
    EXPECT_DOUBLE_EQ(v2[0], 0.5);
    EXPECT_DOUBLE_EQ(v2[1], 1.0);
    EXPECT_DOUBLE_EQ(v2[2], 1.5);
}

TEST(Vector3Tests, Length) {
    Vector3 v(1.0, 2.0, 2.0);
    EXPECT_DOUBLE_EQ(v.length(), 3.0);
}

TEST(Vector3Tests, Normalization) {
    Vector3 v(1.0, 2.0, 2.0);
    Vector3 normalized = v.normalized();
    EXPECT_NEAR(normalized.length(), 1.0, 1e-4);
}

TEST(Vector3Tests, DotProduct) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    double dotProduct = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dotProduct, 32.0);
}

TEST(Vector3Tests, ComponentProduct) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    Vector3 product = v1.componentProd(v2);
    EXPECT_DOUBLE_EQ(product[0], 4.0);
    EXPECT_DOUBLE_EQ(product[1], 10.0);
    EXPECT_DOUBLE_EQ(product[2], 18.0);
}

TEST(Vector3Tests, CrossProductWithOrthogonalVectors) {
    Vector3 v1(1.0, 0.0, 0.0);
    Vector3 v2(0.0, 1.0, 0.0);
    Vector3 crossProduct = v1.cross(v2);
    EXPECT_DOUBLE_EQ(crossProduct[0], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[1], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[2], 1.0);
}

TEST(Vector3Tests, CrossProductWithParallelVectors) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(1.0, 2.0, 3.0);
    Vector3 crossProduct = v1.cross(v2);
    EXPECT_DOUBLE_EQ(crossProduct[0], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[1], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[2], 0.0);
}

TEST(Vector3Tests, CrossProductWithZeroVector) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(0.0, 0.0, 0.0);
    Vector3 crossProduct = v1.cross(v2);
    EXPECT_DOUBLE_EQ(crossProduct[0], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[1], 0.0);
    EXPECT_DOUBLE_EQ(crossProduct[2], 0.0);
}

TEST(Vector4Tests, InitializationWithFourParameters) {
    Vector4 v(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
    EXPECT_DOUBLE_EQ(v[3], 4.0);
}

TEST(Vector4Tests, CopyConstructor) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2(v1);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
    EXPECT_DOUBLE_EQ(v2[3], 4.0);
}

TEST(Vector4Tests, AdditionOperator) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2(5.0, 6.0, 7.0, 8.0);
    Vector4 v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3[0], 6.0);
    EXPECT_DOUBLE_EQ(v3[1], 8.0);
    EXPECT_DOUBLE_EQ(v3[2], 10.0);
    EXPECT_DOUBLE_EQ(v3[3], 12.0);
}

TEST(Vector4Tests, SubtractionOperator) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2(5.0, 6.0, 7.0, 8.0);
    Vector4 v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3[0], -4.0);
    EXPECT_DOUBLE_EQ(v3[1], -4.0);
    EXPECT_DOUBLE_EQ(v3[2], -4.0);
    EXPECT_DOUBLE_EQ(v3[3], -4.0);
}

TEST(Vector4Tests, MultiplicationOperator) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2 = v1 * 3.0;
    EXPECT_DOUBLE_EQ(v2[0], 3.0);
    EXPECT_DOUBLE_EQ(v2[1], 6.0);
    EXPECT_DOUBLE_EQ(v2[2], 9.0);
    EXPECT_DOUBLE_EQ(v2[3], 12.0);
}

TEST(Vector4Tests, DivisionOperator) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2 = v1 / 2.0;
    EXPECT_DOUBLE_EQ(v2[0], 0.5);
    EXPECT_DOUBLE_EQ(v2[1], 1.0);
    EXPECT_DOUBLE_EQ(v2[2], 1.5);
    EXPECT_DOUBLE_EQ(v2[3], 2.0);
}

TEST(Vector4Tests, Length) {
    Vector4 v(1.0, 2.0, 2.0, 2.0);
    EXPECT_NEAR(v.length(), 3.6055, 1e-4);
}

TEST(Vector4Tests, Normalization) {
    Vector4 v(1.0, 2.0, 2.0, 2.0);
    Vector4 normalized = v.normalized();
    EXPECT_NEAR(normalized.length(), 1.0, 1e-4);
}

TEST(Vector4Tests, DotProduct) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2(5.0, 6.0, 7.0, 8.0);
    double dotProduct = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dotProduct, 70.0);
}

TEST(Vector4Tests, ComponentProduct) {
    Vector4 v1(1.0, 2.0, 3.0, 4.0);
    Vector4 v2(5.0, 6.0, 7.0, 8.0);
    Vector4 product = v1.componentProd(v2);
    EXPECT_DOUBLE_EQ(product[0], 5.0);
    EXPECT_DOUBLE_EQ(product[1], 12.0);
    EXPECT_DOUBLE_EQ(product[2], 21.0);
    EXPECT_DOUBLE_EQ(product[3], 32.0);
}