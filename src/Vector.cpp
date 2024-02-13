/**
 * @file Vector.cpp
 * @author ayano
 * @date 13/02/24
 */

#include "Math/Vector.h"

#ifdef __APPLE__
namespace Math {
    Vector2::Vector2(const simd::double2& data) : data(data) {}
    Vector2::Vector2() : Vector2(simd::double2{0.0, 0.0}) {}
    Vector2::Vector2(double x, double y) : Vector2(simd::double2{x, y}) {};
    Vector2::Vector2(const Vector2& other) {
        data = other.data;
    }
    Vector2& Vector2::operator=(const Vector2& other) {
        if (this != &other) {
            data = other.data;
        }
        return *this;
    }
    Vector2 Vector2::operator+(const Vector2& rhs) const {
        return Vector2(data + rhs.data);
    }
    Vector2 Vector2::operator-(const Vector2& rhs) const {
        return Vector2(data - rhs.data);
    }
    Vector2 Vector2::operator*(double rhs) const {
        return Vector2(data * rhs);
    }
    Vector2 Vector2::operator/(double rhs) const {
        return Vector2(data / rhs);
    }
    Vector2 Vector2::operator+=(const Vector2& rhs) {
        data += rhs.data;
        return *this;
    }
    Vector2 Vector2::operator-=(const Vector2& rhs) {
        data -= rhs.data;
        return *this;
    }
    Vector2 Vector2::operator*=(double rhs) {
        data *= rhs;
        return *this;
    }
    Vector2 Vector2::operator/=(double rhs) {
        data /= rhs;
        return *this;
    }
    double Vector2::operator[](int index) const {
        switch(index) {
            case 0:
                return data[0];
            case 1:
                return data[1];
            default:
                throw std::out_of_range("Index out of range");
        }
    }
    bool Vector2::operator==(const Vector2& rhs) const {
        return simd::all(data == rhs.data);
    }
    bool Vector2::operator!=(const Vector2& rhs) const {
        return !(*this == rhs);
    }
    double Vector2::lengthSq() const {
        return simd::length_squared(data);
    }
    double Vector2::length() const {
        return simd::length(data);
    }
    Vector2 Vector2::normalized() const {
        return Vector2(simd::normalize(data));
    }
    double Vector2::dot(const Vector2& rhs) const {
        return simd::dot(data, rhs.data);
    }
    Vector2 Vector2::componentProd(const Vector2& rhs) const {
        return Vector2(data * rhs.data);
    }
    Vector2::operator std::string() const {
        return std::string("(" + std::to_string(data[0]) + ", " + std::to_string(data[1]) + ")");
    }

    Vector3::Vector3() : Vector3(simd::double3{0.0, 0.0, 0.0}) {

    }

    Vector3::Vector3(double x, double y, double z): Vector3(simd::double3{x, y, z}) {

    }

    Vector3::Vector3(const simd::double3 &data) : data(data) {

    }

    Vector3::Vector3(const Vector3 &other) : data(other.data){

    }

    Vector3 &Vector3::operator=(const Vector3 &other) {
        if(this != &other) {
            data = other.data;
        }
        return *this;
    }

    Vector3 Vector3::operator+(const Vector3 &rhs) const {
        return Vector3(data + rhs.data);
    }

    Vector3 Vector3::operator-(const Vector3 &rhs) const {
        return Vector3(data - rhs.data);
    }

    Vector3 Vector3::operator*(double rhs) const {
        return Vector3(data * rhs);
    }

    Vector3 Vector3::operator/(double rhs) const {
        return Vector3(data / rhs);
    }

    Vector3 Vector3::operator+=(const Vector3 &rhs) {
        data += rhs.data;
        return *this;
    }

    Vector3 Vector3::operator-=(const Vector3 &rhs) {
        data -= rhs.data;
        return *this;
    }

    Vector3 Vector3::operator*=(double rhs) {
        data *= rhs;
        return *this;
    }

    Vector3 Vector3::operator/=(double rhs) {
        data /= rhs;
        return *this;
    }

    double Vector3::operator[](int index) const {
        switch(index) {
            case 0:
                return data[0];
            case 1:
                return data[1];
            case 2:
                return data[2];
            default:
                throw std::out_of_range("Index out of range");
        }
    }

    bool Vector3::operator==(const Vector3 &rhs) const {
        return simd::all(data == rhs.data);
    }

    bool Vector3::operator!=(const Vector3 &rhs) const {
        return !(*this == rhs);
    }

    Vector3::operator std::string() const {
        return std::string("(" + std::to_string(data[0]) + ", " + std::to_string(data[1]) + ", " + std::to_string(data[2]) + ")");
    }

    double Vector3::lengthSq() const {
        return simd::length_squared(data);
    }

    double Vector3::length() const {
        return simd::length(data);
    }

    Vector3 Vector3::normalized() const {
        return Vector3(simd::normalize(data));
    }

    double Vector3::dot(const Vector3 &rhs) const {
        return simd::dot(data, rhs.data);
    }

    Vector3 Vector3::componentProd(const Vector3 &rhs) const {
        return Vector3(data * rhs.data);
    }

    Vector3 Vector3::cross(const Vector3 &rhs) const {
        return Vector3(simd::cross(data, rhs.data));
    }

    Vector4::Vector4() : Vector4(simd::double4{0.0, 0.0, 0.0, 0.0}) {

    }

    Vector4::Vector4(double x, double y, double z, double w) : Vector4(simd::double4{x, y, z, w}) {

    }

    Vector4::Vector4(const simd::double4 &data) : data(data) {

    }

    Vector4::Vector4(const Vector4 &other) {
        data = other.data;
    }

    Vector4 &Vector4::operator=(const Vector4 &other) {
        if(this != &other) {
            data = other.data;
        }
        return *this;
    }

    Vector4 Vector4::operator+(const Vector4 &rhs) const {
        return Vector4(data + rhs.data);
    }

    Vector4 Vector4::operator-(const Vector4 &rhs) const {
        return Vector4(data - rhs.data);
    }

    Vector4 Vector4::operator*(double rhs) const {
        return Vector4(data * rhs);
    }

    Vector4 Vector4::operator/(double rhs) const {
        return Vector4(data / rhs);
    }

    Vector4 Vector4::operator+=(const Vector4 &rhs) {
        data += rhs.data;
        return *this;
    }

    Vector4 Vector4::operator-=(const Vector4 &rhs) {
        data -= rhs.data;
        return *this;
    }

    Vector4 Vector4::operator*=(double rhs) {
        data *= rhs;
        return *this;
    }

    Vector4 Vector4::operator/=(double rhs) {
        data /= rhs;
        return *this;
    }

    double Vector4::operator[](int index) const {
        switch(index) {
            case 0:
                return data[0];
            case 1:
                return data[1];
            case 2:
                return data[2];
            case 3:
                return data[3];
            default:
                throw std::out_of_range("Index out of range");
        }
    }

    bool Vector4::operator==(const Vector4 &rhs) const {
        return simd::all(data == rhs.data);
    }

    bool Vector4::operator!=(const Vector4 &rhs) const {
        return !(*this == rhs);
    }

    Vector4::operator std::string() const {
        return std::string("(" + std::to_string(data[0]) + ", " + std::to_string(data[1]) + ", " + std::to_string(data[2]) + ", " + std::to_string(data[3]) + ")");
    }

    double Vector4::lengthSq() const {
        return simd::length_squared(data);
    }

    double Vector4::length() const {
        return simd::length(data);
    }

    Vector4 Vector4::normalized() const {
        return Vector4(simd::normalize(data));
    }

    double Vector4::dot(const Vector4 &rhs) const {
        return simd::dot(data, rhs.data);
    }

    Vector4 Vector4::componentProd(const Vector4 &rhs) const {
        return Vector4(data * rhs.data);
    }
}

#endif
