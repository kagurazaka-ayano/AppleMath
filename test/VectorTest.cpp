#include <gtest/gtest.h>
#include "AppleMath/Vector.hpp"
#include "AppleMath/Matrix.hpp"

using namespace AppleMath;

TEST(VectorTest, ConstructorFromPVector) {
    simd::double3 data = {1.0, 2.0, 3.0};
    Vector<3> v(data);
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, ConstructorFromArray) {
    std::array<double, 3> data = {1.0, 2.0, 3.0};
    Vector<3> v(data);
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, ConstructorFromInitializerList) {
    Vector<3> v({1.0, 2.0, 3.0});
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, CopyConstructor) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2(v1);
    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);
}

TEST(VectorTest, AssignmentOperator) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2 = v1;
    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);
}

TEST(VectorTest, AdditionOperator) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({4.0, 5.0, 6.0});
    Vector<3> v3 = v1 + v2;
    EXPECT_EQ(v3[0], 5.0);
    EXPECT_EQ(v3[1], 7.0);
    EXPECT_EQ(v3[2], 9.0);
}

TEST(VectorTest, SubtractionOperator) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({4.0, 5.0, 6.0});
    Vector<3> v3 = v1 - v2;
    EXPECT_EQ(v3[0], -3.0);
    EXPECT_EQ(v3[1], -3.0);
    EXPECT_EQ(v3[2], -3.0);
}

TEST(VectorTest, MultiplicationOperator) {
    Vector<3> v({1.0, 2.0, 3.0});
    Vector<3> v2 = v * 2.0;
    EXPECT_EQ(v2[0], 2.0);
    EXPECT_EQ(v2[1], 4.0);
    EXPECT_EQ(v2[2], 6.0);
}

TEST(VectorTest, DivisionOperator) {
    Vector<3> v({1.0, 2.0, 3.0});
    Vector<3> v2 = v / 2.0;
    EXPECT_EQ(v2[0], 0.5);
    EXPECT_EQ(v2[1], 1.0);
    EXPECT_EQ(v2[2], 1.5);
}

TEST(VectorTest, Length) {
    Vector<3> v({1.0, 2.0, 2.0});
    EXPECT_EQ(v.length(), 3.0);
}

TEST(VectorTest, Normalized) {
    Vector<3> v({1.0, 2.0, 2.0});
    Vector<3> v2 = v.normalized();
    EXPECT_NEAR(v2.length(), 1.0, 1e-9);
}

TEST(VectorTest, DotProduct) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({4.0, 5.0, 6.0});
    EXPECT_EQ(v1.dot(v2), 32.0);
}

TEST(VectorTest, CrossProduct) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({4.0, 5.0, 6.0});
    Vector<3> v3 = v1.cross(v2);
    EXPECT_EQ(v3[0], -3.0);
    EXPECT_EQ(v3[1], 6.0);
    EXPECT_EQ(v3[2], -3.0);
}

TEST(VectorTest, EqualityOperator) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({1.0, 2.0, 3.0});
    EXPECT_TRUE(v1 == v2);
}

TEST(VectorTest, InequalityOperator) {
    Vector<3> v1({1.0, 2.0, 3.0});
    Vector<3> v2({4.0, 5.0, 6.0});
    EXPECT_TRUE(v1 != v2);
}

TEST(VectorTest, OutOfRangeIndex) {
    Vector<3> v({1.0, 2.0, 3.0});
    EXPECT_THROW(v[3], std::out_of_range);
}

TEST(VectorTest, ScalarMultiplication) {
    Vector3 v{1.0, 2.0, 3.0};
    auto result = 2.0 * v;
    auto cmp = Vector3{2.0, 4.0, 6.0};
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, ScalarMultiplicationZero) {
    Vector3 v{1.0, 2.0, 3.0};
    auto result = 0.0 * v;
    auto cmp = Vector3{0.0, 0.0, 0.0};
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, MakeHomoCoord2D) {
    Vector2 v{1.0, 2.0};
    auto result = makeHomoCoord(v);
    auto cmp = Vector3{1.0, 2.0, 1.0};
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, MakeHomoCoord3D) {
    Vector3 v{1.0, 2.0, 3.0};
    auto result = makeHomoCoord(v);
    auto cmp = Vector4{1.0, 2.0, 3.0, 1.0};
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, ApplyTrans2D) {
    Vector2 v{1.0, 2.0};
    Matrix<2, 2> m{{1.0, 0.0}, {0.0, 1.0}};
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, ApplyTrans3D) {
    Vector3 v{1.0, 2.0, 3.0};
    Matrix<3, 3> m{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, ApplyTrans4D) {
    Vector4 v{1.0, 2.0, 3.0, 4.0};
    Matrix<4, 4> m{{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, Rotation2D) {
    Vector2 v{1.0, 0.0};
    auto rotation = M_PI / 2;
    auto mat = makeRotationMatrixR2(rotation);
    auto result = applyTrans(v, mat);
    EXPECT_NEAR(result[0], 0, EPS);
    EXPECT_EQ(result[1], 1);
}

TEST(VectorTest, Rotation3D) {
    Vector3 v{1, 0, 0};
    auto phi = M_PI / 2;
    auto theta = M_PI / 2;
    auto psi = M_PI / 2;
    auto mat = makeRotationMatrixR3(phi, theta, psi);
    auto result = applyTrans(v, mat["z"]);
    EXPECT_NEAR(result[0], 0, EPS);
    EXPECT_NEAR(result[1], 1, EPS);
    EXPECT_NEAR(result[2], 0, EPS);
    result = applyTrans(result, mat["x"]);
    EXPECT_NEAR(result[0], 0, EPS);
    EXPECT_NEAR(result[1], 0, EPS);
    EXPECT_NEAR(result[2], 1, EPS);
    result = applyTrans(result, mat["y"]);
    EXPECT_NEAR(result[0], 1, EPS);
    EXPECT_NEAR(result[1], 0, EPS);
    EXPECT_NEAR(result[2], 0, EPS);
}
