/**
 * @file Vector.hpp
 * @author ayano
 * @date 2/14/24
 * @brief
*/

#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP
#include <simd/simd.h>
#include <vector>
namespace AppleMath {
    template<std::size_t N>
    requires (N == 2 || N == 3 || N == 4)
    class Vector {
        using p_vector = std::conditional_t<
        N == 2, simd::double2, std::conditional_t<
        N == 3, simd::double3, std::conditional_t<
        N == 4, simd::double4, void>>>;
    public:
        explicit Vector(p_vector data) : data(data) {

        }
        explicit Vector(std::array<double, N> data) {
            for (std::size_t i = 0; i < N; i++) {
                this->data[i] = data[i];
            }
        }
        explicit Vector(std::initializer_list<double> li) {
            if(li.size() != N) {
                throw std::invalid_argument("Initializer list size does not match vector size");
            }
            std::size_t i = 0;
            for (double it : li) {
                data[i] = it;
                i++;
            }
        }
        Vector(const Vector& other) {
            data = other.data;
        }
        Vector& operator=(const Vector& other) {
            if (this != &other) {
                data = other.data;
            }
            return *this;
        }
        Vector operator+(const Vector& rhs) const {
            return Vector(data + rhs.data);
        }
        Vector operator-(const Vector& rhs) const {
            return Vector(data - rhs.data);
        }
        Vector operator-(void) const {
            return Vector(-data);
        }
        Vector operator*(double rhs) const {
            return Vector(data * rhs);
        }
        Vector operator/(double rhs) const {
            return Vector(data / rhs);
        }
        Vector operator+=(const Vector& rhs) {
            data += rhs.data;
            return *this;
        }
        Vector operator-=(const Vector& rhs) {
            data -= rhs.data;
            return *this;
        }
        Vector operator*=(double rhs) {
            data *= rhs;
            return *this;
        }
        Vector operator/=(double rhs) {
            data /= rhs;
            return *this;
        }
        double operator[](int index) const {
            if (index >= N) {
                throw std::out_of_range("Index out of range");
            }
            return data[index];
        }

        double lengthSq() {
            return simd::length_squared(data);
        }

        double length() const {
            return simd::length(data);
        }
        Vector normalized() const {
            return Vector(simd::normalize(data));
        }
        double dot(const Vector& rhs) const {
            return simd::dot(data, rhs.data);
        }
        Vector cross(const Vector& rhs) const {
            if (N == 4) {
                throw std::invalid_argument("Cross product is not defined for 4D vectors");
            }
            return Vector(simd::cross(data, rhs.data));
        }
        Vector componentProd(const Vector& rhs) const {
            return Vector(data * rhs.data);
        }
        bool operator==(const Vector& rhs) const {
            return simd_equal(data, rhs.data);
        }
        bool operator!=(const Vector& rhs) const {
            return !simd_equal(data, rhs.data);
        }
        explicit operator std::string() const {
            std::string str = "(";
            for (std::size_t i = 0; i < N; i++) {
                str += std::to_string(data[i]);
                if (i != N - 1) {
                    str += ", ";
                }
            }
            str += ")";
            return str;
        }
        std::ostream& operator<<(std::ostream& os) const {
            os << static_cast<std::string>(*this);
            return os;
        }
        p_vector getData() const {
            return data;
        }
    private:
        p_vector data;
    };
    inline Vector<2> operator*(double lhs, const Vector<2>& rhs) {
        return rhs * lhs;
    }
    inline Vector<3> operator*(double lhs, const Vector<3>& rhs) {
        return rhs * lhs;
    }
    inline Vector<4> operator*(double lhs, const Vector<4>& rhs) {
        return rhs * lhs;
    }
    using Vector2 = Vector<2>;
    using Vector3 = Vector<3>;
    using Vector4 = Vector<4>;
};
#endif //MATH_VECTOR_HPP