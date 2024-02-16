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

    template<std::size_t R, std::size_t C>
    requires (R == 2 || R == 3 || R == 4) && (C == 2 || C == 3 || C == 4)
    class Matrix;

    template<std::size_t N>
    requires (N == 2 || N == 3 || N == 4)
    class Vector {
        using p_vector = std::conditional_t<
        N == 2, simd::double2, std::conditional_t<
        N == 3, simd::double3, std::conditional_t<
        N == 4, simd::double4, void>>>;
    public:
        /**
         * default constructor, constructor a vector with all element be 0
         */
        explicit Vector() {
            if constexpr (N == 2) {
                data = simd::double2{0.0, 0.0};
            } else if constexpr (N == 3) {
                data = simd::double3{0.0, 0.0, 0.0};
            } else if constexpr (N == 4) {
                data = simd::double4{0.0, 0.0, 0.0, 0.0};
            }
        }

        /**
         * constructor
         * @param data another raw apple simd vector(must have 2 ~ 4 elements)
         */
        explicit Vector(p_vector data) : data(data) {

        }

        /**
         * constructor
         * @param data array with elements
         */
        explicit Vector(std::array<double, N> data) {
            for (std::size_t i = 0; i < N; i++) {
                this->data[i] = data[i];
            }
        }

        /**
         * constructor
         * @param li initializer list with elements
         */
        Vector(std::initializer_list<double> li) {
            if(li.size() != N) {
                throw std::invalid_argument("Initializer list size does not match vector size");
            }
            std::size_t i = 0;
            for (double it : li) {
                data[i] = it;
                i++;
            }
        }

        /**
         * copy constructor
         * @param other another vector
         */
        Vector(const Vector<N>& other) {
            data = other.data;
        }

        /**
         * equality
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
         * add
         * @param rhs another vector
         * @return result
         */
        Vector operator+(const Vector<N>& rhs) const {
            return Vector(data + rhs.data);
        }

        /**
         * minus
         * @param rhs another vector
         * @return result
         */
        Vector operator-(const Vector<N>& rhs) const {
            return Vector(data - rhs.data);
        }

        /**
         * negation
         * @return this vector with every element negated
         */
        Vector operator-() const {
            return Vector(-data);
        }

        /**
         * scalar multiplication
         * @param rhs scalar
         * @return result
         */
        Vector operator*(double rhs) const {
            return Vector(data * rhs);
        }

        /**
         * scalar division
         * @param rhs scalar
         * @return result
         */
        Vector operator/(double rhs) const {
            return Vector(data / rhs);
        }

        /**
         * self increment
         * @param rhs another vector
         * @return reference to this vector
         */
        Vector& operator+=(const Vector<N>& rhs) {
            data += rhs.data;
            return *this;
        }

        /**
         * self decrement
         * @param rhs another vector
         * @return reference to this vector
         */
        Vector& operator-=(const Vector<N>& rhs) {
            data -= rhs.data;
            return *this;
        }

        /**
         * self scalar multiplication
         * @param rhs scalar
         * @return reference to this vector
         */
        Vector& operator*=(double rhs) {
            data *= rhs;
            return *this;
        }

        /**
         * self scalar division
         * @param rhs scalar
         * @return reference to this vector
         */
        Vector& operator/=(double rhs) {
            data /= rhs;
            return *this;
        }

        /**
         * element getter
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
         * element getter
         * @param index index
         * @return element at index position
         */
        double operator()(std::size_t index) const {
            if (index >= N) {
                throw std::out_of_range("Index out of range");
            }
            return data[index];
        }

        void setElement(std::size_t idx, double val) {
            if (N <= idx) throw std::out_of_range("index out of range");
            data[idx] = val;
        }


        /**
         * vector length squared
         * @return vector length squared
         */
        double lengthSq() {
            return simd::length_squared(data);
        }

        /**
         * vector length
         * @return vector length
         */
        double length() const {
            return simd::length(data);
        }

        /**
         * return normalized (length = 1) version of this vector
         * @return this vector, normalized
         */
        Vector normalized() const {
            return Vector(simd::normalize(data));
        }

        /**
         * dot product
         * @param rhs another vector with same size
         * @return result
         */
        double dot(const Vector<N>& rhs) const {
            return simd::dot(data, rhs.data);
        }

        /**
         * cross product
         * @param rhs another vector with same size and less than 4
         * @return result vector
         * @remark if vector is in R2, 0 will be added as the last element
         */
        Vector cross(const Vector<N>& rhs) const {
            static_assert(N != 4, "Cross product is not defined for 4D vectors");
            return Vector(simd::cross(data, rhs.data));
        }

        /**
         * component product
         * @param rhs another vector with same size
         * @return a vector with each element multiplied with another
         */
        Vector componentProd(const Vector<N>& rhs) const {
            return Vector(data * rhs.data);
        }

        /**
         * equality
         * @param rhs another vector
         * @return true if equal
         */
        bool operator==(const Vector<N>& rhs) const {
            return simd_equal(data, rhs.data);
        }

        /**
         * inequality
         * @param rhs another vector
         * @return false if equal
         */
        bool operator!=(const Vector<N>& rhs) const {
            return !simd_equal(data, rhs.data);
        }

        /**
         * toString
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
         * stream operator
         * @param os stream
         * @return stream
         */
        std::ostream& operator<<(std::ostream& os) const {
            os << static_cast<std::string>(*this);
            return os;
        }

        /**
         * get dimension of the vector
         * @return dimension
         */
        std::size_t size() {
            return N;
        }

        /**
         * get the raw apple simd vector version of this vector
         * @return raw apple simd vector version of this vector
         */
        p_vector getData() const {
            return data;
        }
    private:
        p_vector data;
    };

    /**
     * scalar multiplication
     * @param lhs scalar
     * @param rhs vector
     * @tparam N dimension
     * @return result
     */
    template<std::size_t N>
    requires (N == 2 || N == 3 || N == 4)
    inline Vector<N> operator*(double lhs, const Vector<N>& rhs) {
        return rhs * lhs;
    }

    /**
     * make a vector to its P^N correspondence
     * @tparam N dimension
     */
    template<std::size_t N>
    requires (N == 2) || (N == 3)
    Vector<N + 1> makeHomoCoord(Vector<N> vec) {
        Vector<N + 1> ret;
        for(int i = 0; i < N; i++) {
            ret.setElement(i, vec[i]);
        }
        ret.setElement(N, 1);
        return ret;
    }

    /**
     * extract the R3/R2 value from the homogeneous coordinate
     * @tparam N
     */
    template<std::size_t N>
    requires (N == 3) || (N == 4)
    Vector<N - 1> deHomonize(Vector<N> vec) {
        Vector<N - 1> ret;
        for(int i = 0; i < N - 1; i++) {
            ret.setElement(i, vec[i]);
        }
        return ret;
    }

    /**
     * apply a matrix transformation to the given vector
     * @tparam N dimension
     */
    template<std::size_t N>
    requires (N == 2) || (N == 3) || (N == 4)
    inline Vector<N> applyTrans(Vector<N> lhs, Matrix<N, N> rhs) {
        return rhs * lhs;
    }

    using Vector2 = Vector<2>;
    using Vector3 = Vector<3>;
    using Vector4 = Vector<4>;
};
#endif //MATH_VECTOR_HPP