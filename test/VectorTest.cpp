#include "AppleMath/Vector.hpp"
#include "AppleMath/Matrix.hpp"
#include <gtest/gtest.h>

using namespace AppleMath;

TEST(VectorTest, ConstructorFromPVector)
{
    simd::float3 data = { 1.0, 2.0, 3.0 };
    Vector<3> v(data);
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, ConstructorFromArray)
{
    std::array<float, 3> data = { 1.0, 2.0, 3.0 };
    Vector<3> v(data);
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, ConstructorFromInitializerList)
{
    Vector<3> v({ 1.0, 2.0, 3.0 });
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
}

TEST(VectorTest, CopyConstructor)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2(v1);
    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);
}

TEST(VectorTest, AssignmentOperator)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2 = v1;
    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);
}

TEST(VectorTest, AdditionOperator)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 4.0, 5.0, 6.0 });
    Vector<3> v3 = v1 + v2;
    EXPECT_EQ(v3[0], 5.0);
    EXPECT_EQ(v3[1], 7.0);
    EXPECT_EQ(v3[2], 9.0);
}

TEST(VectorTest, SubtractionOperator)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 4.0, 5.0, 6.0 });
    Vector<3> v3 = v1 - v2;
    EXPECT_EQ(v3[0], -3.0);
    EXPECT_EQ(v3[1], -3.0);
    EXPECT_EQ(v3[2], -3.0);
}

TEST(VectorTest, MultiplicationOperator)
{
    Vector<3> v({ 1.0, 2.0, 3.0 });
    Vector<3> v2 = v * 2.0;
    EXPECT_EQ(v2[0], 2.0);
    EXPECT_EQ(v2[1], 4.0);
    EXPECT_EQ(v2[2], 6.0);
}

TEST(VectorTest, DivisionOperator)
{
    Vector<3> v({ 1.0, 2.0, 3.0 });
    Vector<3> v2 = v / 2.0;
    EXPECT_EQ(v2[0], 0.5);
    EXPECT_EQ(v2[1], 1.0);
    EXPECT_EQ(v2[2], 1.5);
}

TEST(VectorTest, Length)
{
    Vector<3> v({ 1.0, 2.0, 2.0 });
    EXPECT_EQ(v.length(), 3.0);
}

TEST(VectorTest, Normalized)
{
    Vector<3> v({ 1.0, 2.0, 2.0 });
    Vector<3> v2 = v.normalized();
    EXPECT_NEAR(v2.length(), 1.0, 1e-4);
}

TEST(VectorTest, DotProduct)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 4.0, 5.0, 6.0 });
    EXPECT_EQ(v1.dot(v2), 32.0);
}

TEST(VectorTest, CrossProduct)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 4.0, 5.0, 6.0 });
    Vector<3> v3 = v1.cross(v2);
    EXPECT_EQ(v3[0], -3.0);
    EXPECT_EQ(v3[1], 6.0);
    EXPECT_EQ(v3[2], -3.0);
}

TEST(VectorTest, EqualityOperator)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 1.0, 2.0, 3.0 });
    EXPECT_TRUE(v1 == v2);
}

TEST(VectorTest, InequalityOperator)
{
    Vector<3> v1({ 1.0, 2.0, 3.0 });
    Vector<3> v2({ 4.0, 5.0, 6.0 });
    EXPECT_TRUE(v1 != v2);
}

TEST(VectorTest, OutOfRangeIndex)
{
    Vector<3> v({ 1.0, 2.0, 3.0 });
    EXPECT_THROW(v[3], std::out_of_range);
}

TEST(VectorTest, ScalarMultiplication)
{
    Vector3 v { 1.0, 2.0, 3.0 };
    auto result = 2.0 * v;
    auto cmp = Vector3 { 2.0, 4.0, 6.0 };
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, ScalarMultiplicationZero)
{
    Vector3 v { 1.0, 2.0, 3.0 };
    auto result = 0.0 * v;
    auto cmp = Vector3 { 0.0, 0.0, 0.0 };
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, MakeHomoCoord2D)
{
    Vector2 v { 1.0, 2.0 };
    auto result = makeHomoCoord(v);
    auto cmp = Vector3 { 1.0, 2.0, 1.0 };
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, MakeHomoCoord3D)
{
    Vector3 v { 1.0, 2.0, 3.0 };
    auto result = makeHomoCoord(v);
    auto cmp = Vector4 { 1.0, 2.0, 3.0, 1.0 };
    EXPECT_EQ(result, cmp);
}

TEST(VectorTest, ApplyTrans2D)
{
    Vector2 v { 1.0, 2.0 };
    Matrix<2, 2> m { { 1.0, 0.0 }, { 0.0, 1.0 } };
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, ApplyTrans3D)
{
    Vector3 v { 1.0, 2.0, 3.0 };
    Matrix<3, 3> m { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, ApplyTrans4D)
{
    Vector4 v { 1.0, 2.0, 3.0, 4.0 };
    Matrix<4, 4> m { { 1.0, 0.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0, 0.0 },
        { 0.0, 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 0.0, 1.0 } };
    auto result = applyTrans(v, m);
    EXPECT_EQ(result, v);
}

TEST(VectorTest, ApplyTrans_Rotation2D) {
    Vector<2> v{1.0, 0.0};
    float rotation = M_PI / 2;
    auto mat = makeRotationMatrixR2(rotation);
    Vector<2> result = applyTrans(v, mat);
    EXPECT_NEAR(result[0], 0.0, 1e-4);
    EXPECT_NEAR(result[1], 1.0, 1e-4);
}

TEST(VectorTest, ApplyTrans_Rotation3D) {
    Vector<3> v{1.0, 0.0, 0.0};
    float phi = M_PI / 2;
    float theta = M_PI / 2;
    float psi = M_PI / 2;
    auto mat = makeEulerRotationMatrixR3(phi, theta, psi);
    Vector<3> result = applyTrans(v, mat);
    EXPECT_NEAR(result[0], 0.0, 1e-4);
    EXPECT_NEAR(result[1], 0.0, 1e-4);
    EXPECT_NEAR(result[2], -1.0, 1e-4);
}


TEST(VectorTest, AngleBetweenR3_SameDirection) {
    Vector<3> v1{1.0, 0.0, 0.0};
    Vector<3> v2{1.0, 0.0, 0.0};
    auto angles = getAngleBetweenR3(v1, v2);
    EXPECT_DOUBLE_EQ(angles["psi"], 0);
    EXPECT_DOUBLE_EQ(angles["theta"], 0);
    EXPECT_DOUBLE_EQ(angles["phi"], 0);
}

TEST(VectorTest, AngleBetweenR3_OppositeDirection) {
    Vector<3> v1{1.0, 1, 1};
    Vector<3> v2{-1.0, -1, -1};
    auto angles = getAngleBetweenR3(v1, v2);
    v1 = applyTrans(v1, makeEulerRotationMatrixR3(angles["psi"], angles["theta"], angles["phi"]));
    EXPECT_NEAR(v1[0], v2[0], EPS);
    EXPECT_NEAR(v1[1], v2[1], EPS);
    EXPECT_NEAR(v1[2], v2[2], EPS);
}

TEST(VectorTest, AngleBetweenR3_DifferentDirection) {
    Vector<3> v1{1.0, 0.0, 0.0};
    Vector<3> v2{0.0, 1.0, 0.0};
    auto angles = getAngleBetweenR3(v1, v2);
    v1 = applyTrans(v1, makeEulerRotationMatrixR3(angles["psi"], angles["theta"], angles["phi"]));
    EXPECT_NEAR(v1[0], v2[0], EPS);
    EXPECT_NEAR(v1[1], v2[1], EPS);
    EXPECT_NEAR(v1[2], v2[2], EPS);
}

TEST(VectorTest, AngleBetweenR2)
{
    Vector<2> v1 { 1.0, 0.0 };
    Vector<2> v2 { 0.0, 1.0 };
    float angle = getAngleBetweenR2(v1, v2);
    // Check if the function correctly calculates the angle between the vectors
    EXPECT_NEAR(angle, M_PI_2, 1e-4);
}

TEST(VectorTest, AngleBetweenR2_SpecialAngles)
{
    Vector<2> v1 { 1.0, 0.0 };
    Vector<2> v2 { 1.0, 0.0 };
    float angle = getAngleBetweenR2(v1, v2);
    // Check if the function correctly calculates the angle between the vectors
    // when they are the same
    EXPECT_NEAR(angle, 0.0, 1e-4);
}
