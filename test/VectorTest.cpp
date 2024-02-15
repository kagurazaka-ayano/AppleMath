#include <gtest/gtest.h>
#include "Math/Vector.hpp"

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