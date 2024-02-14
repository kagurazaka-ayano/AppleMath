/**
 * @file Vector.h
 * @author ayano
 * @date 13/02/24
 * @brief A vector class
 * @remark If the host is apple, apple simd is used.
 */

#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP

#include <simd/simd.h>
#include <string>

#ifdef __APPLE__
namespace Math {

    class Vector2;

    class Vector3;

    class Vector4;

    class Vector2 {
    public:
        Vector2();

        Vector2(double x, double y);

        explicit Vector2(const simd::double2& data);

        Vector2(const Vector2& other);

        Vector2& operator=(const Vector2& other);

        Vector2 operator+(const Vector2& rhs) const;

        Vector2 operator-(const Vector2& rhs) const;

        Vector2 operator*(double rhs) const;

        Vector2 operator/(double rhs) const;

        Vector2 operator+=(const Vector2& rhs);

        Vector2 operator-=(const Vector2& rhs);

        Vector2 operator*=(double rhs);

        Vector2 operator/=(double rhs);

        double operator[](int index) const;

        bool operator==(const Vector2& rhs) const;

        bool operator!=(const Vector2& rhs) const;

        explicit operator std::string() const;

        [[nodiscard]] double lengthSq() const;

        [[nodiscard]] double length() const;

        [[nodiscard]] Vector2 normalized() const;

        [[nodiscard]] double dot(const Vector2& rhs) const;

        [[nodiscard]] Vector3 cross(const Vector2& rhs) const;
        
        [[nodiscard]] Vector2 componentProd(const Vector2& rhs) const;

        [[nodiscard]] Vector2 unit() const;

    private:
        simd::double2 data;
    };
    inline Vector2 operator*(double lhs, const Vector2& rhs) {
        return rhs * lhs;
    }

    inline std::ostream& operator<<(std::ostream& os, const Vector2& v) {
        os << static_cast<std::string>(v);
        return os;
    }

    class Vector3 {
    public:
        Vector3();

        Vector3(double x, double y, double z);

        Vector3(const Vector2& other, double z);

        explicit Vector3(const simd::double3& data);

        Vector3(const Vector3& other);

        Vector3& operator=(const Vector3& other);

        Vector3 operator+(const Vector3& rhs) const;

        Vector3 operator-(const Vector3& rhs) const;

        Vector3 operator*(double rhs) const;

        Vector3 operator/(double rhs) const;

        Vector3 operator+=(const Vector3& rhs);

        Vector3 operator-=(const Vector3& rhs);

        Vector3 operator*=(double rhs);

        Vector3 operator/=(double rhs);

        double operator[](int index) const;

        bool operator==(const Vector3& rhs) const;

        bool operator!=(const Vector3& rhs) const;

        explicit operator std::string() const;

        [[nodiscard]] double lengthSq() const;

        [[nodiscard]] double length() const;

        [[nodiscard]] Vector3 normalized() const;

        [[nodiscard]] Vector3 cross(const Vector3& rhs) const;

        [[nodiscard]] double dot(const Vector3& rhs) const;

        [[nodiscard]] Vector3 componentProd(const Vector3& rhs) const;

        [[nodiscard]] Vector3 unit() const;
    private:
        simd::double3 data;
    };

    inline Vector3 operator*(double lhs, const Vector3& rhs) {
        return rhs * lhs;
    }

    inline std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        os << static_cast<std::string>(v);
        return os;
    }

    class Vector4 {
    public:
        Vector4();

        Vector4(double x, double y, double z, double w);

        explicit Vector4(const simd::double4& data);

        Vector4(const Vector4& other);

        Vector4& operator=(const Vector4& other);

        Vector4 operator+(const Vector4& rhs) const;

        Vector4 operator-(const Vector4& rhs) const;

        Vector4 operator*(double rhs) const;

        Vector4 operator/(double rhs) const;

        Vector4 operator+=(const Vector4& rhs);

        Vector4 operator-=(const Vector4& rhs);

        Vector4 operator*=(double rhs);

        Vector4 operator/=(double rhs);

        double operator[](int index) const;

        bool operator==(const Vector4& rhs) const;

        bool operator!=(const Vector4& rhs) const;

        explicit operator std::string() const;

        [[nodiscard]] double lengthSq() const;

        [[nodiscard]] double length() const;

        [[nodiscard]] Vector4 normalized() const;

        [[nodiscard]] double dot(const Vector4& rhs) const;

        [[nodiscard]] Vector4 componentProd(const Vector4& rhs) const;

        [[nodiscard]] Vector4 unit() const;
    private:
        simd::double4 data;
    };

    inline Vector4 operator*(double lhs, const Vector4& rhs) {
        return rhs * lhs;
    }

    inline std::ostream& operator<<(std::ostream& os, const Vector4& v) {
        os << static_cast<std::string>(v);
        return os;
    }

}
#endif

#endif // MATH_VECTOR_HPP