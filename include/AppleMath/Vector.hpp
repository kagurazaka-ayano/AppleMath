/**
 * @file Vector.hpp
 * @author ayano
 * @date 2/14/24
 * @brief
 */

#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP
#include "AppleMath/Matrix.hpp"
#include "AppleMath/configuration.hpp"
#include <cmath>
#include <iostream>
#include <simd/conversion.h>
#include <simd/simd.h>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace AppleMath {

template <std::size_t R, std::size_t C>
    requires(R == 2 || R == 3 || R == 4) && (C == 2 || C == 3 || C == 4)
class Matrix;

template <std::size_t N>
    requires(N == 2 || N == 3 || N == 4)
class Vector {
    using p_vector = std::conditional_t<
        N == 2, simd::double2,
        std::conditional_t<N == 3, simd::double3,
            std::conditional_t<N == 4, simd::double4, void>>>;

public:
    /**
     * @brief default constructor, constructor a vector with all element be 0
     */
    explicit Vector() {
        if constexpr (N == 2) {
            data = simd::double2 { 0.0, 0.0 };
        } else if constexpr (N == 3) {
            data = simd::double3 { 0.0, 0.0, 0.0 };
        } else if constexpr (N == 4) {
            data = simd::double4 { 0.0, 0.0, 0.0, 0.0 };
        }
    }

    /**
     * @brief constructor
     * @param data another raw apple simd vector(must have 2 ~ 4 elements)
     */
    explicit Vector(const p_vector& data)
        : data(data) {
    }

    /**
     * constructor
     * @param data array with elements
     */
    explicit Vector(const std::array<double, N>& data) {
        for (std::size_t i = 0; i < N; i++) {
            this->data[i] = data[i];
        }
    }

    /**
     * @brief constructor
     * @param li initializer list with elements
     */
    Vector(const std::initializer_list<double>& li) {
        if (li.size() != N) {
            throw std::invalid_argument(
                "Initializer list size does not match vector size");
        }
        std::size_t i = 0;
        for (double it : li) {
            data[i] = it;
            i++;
        }
    }

    /**
     * @brief copy constructor
     * @param other another vector
     */
    Vector(const Vector<N>& other)
        : data(other.data) { }

    /**
     * @brief equality
     * @param other another vector
     * @return reference to this vector
     */
    Vector& operator=(const Vector<N>& other) {
        if (this != &other) {
            data = other.data;
        }
        return *this;
    }

    /**
     * @brief add
     * @param rhs another vector
     * @return result
     */
    Vector operator+(const Vector<N>& rhs) const {
        return Vector(data + rhs.data);
    }

    /**
     * @brief minus
     * @param rhs another vector
     * @return result
     */
    Vector operator-(const Vector<N>& rhs) const {
        return Vector(data - rhs.data);
    }

    /**
     * @brief negation
     * @return this vector with every element negated
     */
    Vector operator-() const { return Vector(-data); }

    /**
     * @brief scalar multiplication
     * @param rhs scalar
     * @return result
     */
    Vector operator*(double rhs) const { return Vector(data * rhs); }

    /**
     * @brief scalar division
     * @param rhs scalar
     * @return result
     */
    Vector operator/(double rhs) const { return Vector(data / rhs); }

    /**
     * @brief self increment
     * @param rhs another vector
     * @return reference to this vector
     */
    Vector& operator+=(const Vector<N>& rhs) {
        data += rhs.data;
        return *this;
    }

    /**
     * @brief self decrement
     * @param rhs another vector
     * @return reference to this vector
     */
    Vector& operator-=(const Vector<N>& rhs) {
        data -= rhs.data;
        return *this;
    }

    /**
     * @brief self scalar multiplication
     * @param rhs scalar
     * @return reference to this vector
     */
    Vector& operator*=(double rhs) {
        data *= rhs;
        return *this;
    }

    /**
     * @brief self scalar division
     * @param rhs scalar
     * @return reference to this vector
     */
    Vector& operator/=(double rhs) {
        data /= rhs;
        return *this;
    }

    /**
     * @brief element getter
     * @param index index
     * @return element at index position
     */
    double operator[](std::size_t index) const {
        if (index >= N) {
            throw std::out_of_range("Index out of range");
        }
        return data[index];
    }

    /**
     * @brief element getter
     * @param idx index
     * @return element at index position
     */
    double operator()(std::size_t idx) const {
        if (idx >= N) {
            throw std::out_of_range("Index out of range");
        }
        return data[idx];
    }

    /**
     * @brief self multiplication with matrix
     *
     * @param mat matrix, must have N columns
     * @return the result vector
     */
    template <std::size_t C>
        requires(C >= 2) && (C <= 4)
    Vector<N> operator*=(const Matrix<C, N>& mat) {
        return simd_mul(mat, this->data);
    }

    /**
     * @brief x component
     *
     * @return x component
     */
    double x() const {
        return *this[0];
    }

    /**
     * @brief y component
     *
     * @return y component
     */
    double y() const {
        return *this[1];
    }

    /**
     * @brief z component
     *
     * @return z component
     */
    double z() const {
        if constexpr (N < 3) {
            throw std::out_of_range("vector of size " + std::to_string(N) + " doesn't have z component");
        }
        return *this[2];
    }

    /**
     * @brief w component
     *
     * @return w component
     */
    double w() const {
        if constexpr (N < 4) {
            throw std::out_of_range("vector of size " + std::to_string(N) + " doesn't have z component");
        }
        return *this[3];
    }

    /**
     * @brief element setter
     *
     * @param idx index
     * @param val value you want to set
     */
    void setElement(std::size_t idx, double val) {
        if (N <= idx)
            throw std::out_of_range("index out of range");
        data[idx] = val;
    }

    /**
     * @brief vector length squared
     * @return vector length squared
     */
    double lengthSq() const { return simd::length_squared(data); }

    /**
     * @brief vector length
     * @return vector length
     */
    double length() const { return simd::length(data); }

    /**
     * @brief return normalized (length = 1) version of this vector
     * @return this vector, normalized
     */
    Vector normalized() const { return Vector(simd::normalize(data)); }

    /**
     * @brief dot product
     * @param rhs another vector with same size
     * @return result
     */
    double dot(const Vector<N>& rhs) const { return simd::dot(data, rhs.data); }

    /**
     * @brief cross product
     * @param rhs another vector with same size and less than 4
     * @return result vector
     * @remark if vector is in R2, 0 will be added as the last element
     */
    Vector cross(const Vector<N>& rhs) const {
        static_assert(N != 4, "Cross product is not defined for 4D vectors");
        return Vector(simd::cross(data, rhs.data));
    }

    /**
     * @brief component product
     * @param rhs another vector with same size
     * @return a vector with each element multiplied with another
     */
    Vector componentProd(const Vector<N>& rhs) const {
        return Vector(data * rhs.data);
    }

    /**
     * @brief equality
     * @param rhs another vector
     * @return true if equal
     */
    bool operator==(const Vector<N>& rhs) const {
        return simd_equal(data, rhs.data);
    }

    /**
     * @brief inequality
     * @param rhs another vector
     * @return false if equal
     */
    bool operator!=(const Vector<N>& rhs) const {
        return !simd_equal(data, rhs.data);
    }

    /**
     * @brief toString
     * @return string value of this vector
     */
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

    /**
     * @brief get size of the vector
     * @return size
     */
    std::size_t size() const { return N; }

    /**
     * @brief get the raw apple simd vector version of this vector
     * @return raw apple simd vector version of this vector
     */
    p_vector getData() const { return data; }

private:
    p_vector data;
};

/**
 * @brief scalar multiplication
 * @param lhs scalar
 * @param rhs vector
 * @tparam N dimension
 * @return result
 */
template <std::size_t N>
    requires(N == 2 || N == 3 || N == 4)
inline Vector<N> operator*(double lhs, const Vector<N>& rhs) {
    return rhs * lhs;
}

/**
 * @brief make a vector to its P^N correspondence
 * @tparam N dimension
 */
template <std::size_t N>
    requires(N == 2) || (N == 3)
Vector<N + 1> makeHomoCoord(const Vector<N>& vec) {
    Vector<N + 1> ret;
    for (int i = 0; i < N; i++) {
        ret.setElement(i, vec[i]);
    }
    ret.setElement(N, 1);
    return ret;
}

/**
 * @brief extract the R3/R2 value from the homogeneous coordinate
 * @tparam N
 */
template <std::size_t N>
    requires(N == 3) || (N == 4)
Vector<N - 1> deHomonize(const Vector<N>& vec) {

    Vector<N - 1> ret;
    for (int i = 0; i < N - 1; i++) {
        ret.setElement(i, vec[i]);
    }
    return ret;
}

/**
 * @brief apply a matrix transformation to the given vector
 * @tparam N dimension
 * @remark This is the recommended way to make matrix transformation since it checks the dimension
 */
template <std::size_t N>
    requires(N == 2) || (N == 3) || (N == 4)
inline Vector<N> applyTrans(const Vector<N>& lhs, const Matrix<N, N>& rhs) {
    return rhs * lhs;
}

/**
 * @brief stream operator
 * @param os stream
 * @return stream
 */
template <std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const Vector<N>& v) {
    os << std::string(v);
    return os;
}

/**
 * @brief Get the Angle Between 2 R3 vectors
 *
 * @param from one vector
 * @param to another vector
 * @return a key-value mapping, key is the rotation angle name, value is the rotation in radians
 */
inline std::unordered_map<std::string, double> getAngleBetweenR3(const Vector<3>& from, const Vector<3>& to) {
    auto from_n = from.normalized();
    auto to_n = to.normalized();
    Vector<3> axis;
    double angle;
    // same direction
    if (simd::fabs(from_n.dot(to_n) - 1) < EPS) {
        return {
            { "psi", 0 },
            { "theta", 0 },
            { "phi", 0 }
        };
    } else if (simd::fabs(from_n.dot(to_n) + 1) < EPS) {
        angle = M_PI;
        // kernal of the linear system used to find the vecthr such that the dot product is 0
        axis = { -from_n[1] * 1 / from_n[0] - from_n[2] * 1 / from_n[0], 1, 1 };
        axis = axis.normalized();
    } else {
        // orthogonal axis
        axis = from.cross(to).normalized();
        angle = simd::acos(from.normalized().dot(to.normalized()));
    }
    // by Rodrigues' rotation formula
    auto K = Matrix<3, 3> {
        { 0, axis[2], -axis[1] },
        { -axis[2], 0, axis[0] },
        { axis[1], -axis[0], 0 }
    };
    auto R = makeIdentity<3>() + simd::sin(angle) * K + (1 - simd::cos(angle)) * K * K;
    double psi, theta, phi;
    // rotation matrix decomposition
    // if there is a lock
    if (R[2, 0] == 1 || R[2, 0] == -1) {
        phi = 0;
        if (R[2, 0] == 1) {
            theta = M_PI_2;
            psi = phi + simd::atan2(R[0, 1], R[0, 2]);
        } else {
            theta = -M_PI_2;
            psi = -phi + simd::atan2(-R[0, 1], -R[0, 2]);
        }
    } else {
        theta = -simd::asin(R[2, 0]);
        double c_theta = simd::cos(theta);
        psi = simd::atan2(R[2, 1] / c_theta, R[2, 2] / c_theta);
        phi = simd::atan2(R[1, 0] / c_theta, R[0, 0] / c_theta);
    }
    // check for very small value
    if (simd::fabs(psi) < EPS)
        psi = 0;
    if (simd::fabs(theta) < EPS)
        theta = 0;
    if (simd::fabs(phi) < EPS)
        phi = 0;
    return {
        { "psi", psi },
        { "theta", theta },
        { "phi", phi }
    };
}

/**
 * @brief Get the Angle Between 2 R2 vectors
 *
 * @param from one vector
 * @param to another vector
 * @return the rotation in radians
 */
inline double getAngleBetweenR2(Vector<2> from, Vector<2> to) {
    double dot = from.dot(to);
    return simd::acos(dot / from.length() / to.length());
}

using Vector2 = Vector<2>;
using Vector3 = Vector<3>;
using Vector4 = Vector<4>;
}; // namespace AppleMath
#endif // MATH_VECTOR_HPP