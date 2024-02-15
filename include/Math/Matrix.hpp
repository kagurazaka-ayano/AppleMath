/**
 * @file Matrix.h
 * @author ayano
 * @date 2/14/24
 * @brief
*/

#ifndef MATH_MATRIX_HPP
#define MATH_MATRIX_HPP

#include <vector>
#include <simd/simd.h>
#include "Math/configuration.hpp"
#include "Math/Vector.hpp"

namespace AppleMath{


    /**
     * @brief A column major class, only support 2x2, 2x3, 2x4, 3x2, 3x3, 3x4, 4x2, 4x3, 4x4 matrices
     * @tparam R row count
     * @tparam C column count
     */
    template<std::size_t R, std::size_t C>
    requires (R == 2 || R == 3 || R == 4) && (C == 2 || C == 3 || C == 4)
    class Matrix {

        /**
         * @brief CheckSize struct, used to determine the size of the matrix
         */
        struct CheckSize{
            static constexpr bool row2 = R == 2;
            static constexpr bool row3 = R == 3;
            static constexpr bool row4 = R == 4;
            static constexpr bool col2 = C == 2;
            static constexpr bool col3 = C == 3;
            static constexpr bool col4 = C == 4;
        };

        static constexpr bool is2x2 = CheckSize::row2 && CheckSize::col2;
        static constexpr bool is2x3 = CheckSize::row2 && CheckSize::col3;
        static constexpr bool is2x4 = CheckSize::row2 && CheckSize::col4;
        static constexpr bool is3x2 = CheckSize::row3 && CheckSize::col2;
        static constexpr bool is3x3 = CheckSize::row3 && CheckSize::col3;
        static constexpr bool is3x4 = CheckSize::row3 && CheckSize::col4;
        static constexpr bool is4x2 = CheckSize::row4 && CheckSize::col2;
        static constexpr bool is4x3 = CheckSize::row4 && CheckSize::col3;
        static constexpr bool is4x4 = CheckSize::row4 && CheckSize::col4;
        using p_matrix = std::conditional_t<
        is2x2, simd::double2x2, std::conditional_t<
        is2x3, simd::double2x3, std::conditional_t<
        is2x4, simd::double2x4, std::conditional_t<
        is3x2, simd::double3x2, std::conditional_t<
        is3x3, simd::double3x3, std::conditional_t<
        is3x4, simd::double3x4, std::conditional_t<
        is4x2, simd::double4x2, std::conditional_t<
        is4x3, simd::double4x3, std::conditional_t<
        is4x4, simd::double4x4, void>>>>>>>>>;
    public:
        explicit Matrix(p_matrix data) : data(data) {

        }
        explicit Matrix(std::array<std::array<double, C>, R> data) {
            for (std::size_t i = 0; i < R; i++) {
                for (std::size_t j = 0; j < C; j++) {
                    this->data.columns[i][j] = data[i][j];
                }
            }
        }
        explicit Matrix(std::initializer_list<std::initializer_list<double>> li) {
            int col_idx = 0;
            if (li.size() != R) throw std::invalid_argument("Invalid size");
            for (auto col : li) {
                if (col.size() != C) throw std::invalid_argument("Invalid size");
                int row_idx = 0;
                for (auto row : col) {
                    data.columns[col_idx][row_idx] = row;
                    row_idx++;
                }
                col_idx++;
            }
        }
        Matrix(const Matrix<R, C>& other) {
            data = other.data;
        }
        Matrix(Matrix<R, C>&& other) noexcept {
            data = other.data;
        }
        Matrix& operator=(const Matrix<R, C>& other) {
            if (this != &other) {
                data = other.data;
            }
            return *this;
        }
        Matrix operator+(const Matrix<R, C>& rhs) const {
            return Matrix(simd_add(data, rhs.getData()));
        }
        Matrix operator-(const Matrix<R, C>& rhs) const {
            return Matrix(simd_sub(data, rhs.getData()));
        }
        Matrix operator*(double rhs) const {
            return Matrix(simd_mul(rhs, data));
        }
        Matrix<R, R> operator*(Matrix<C, R> rhs) const {
            return Matrix<R, R>(simd_mul(data, rhs.getData()));
        }
        Vector<R> operator*(Vector<R> rhs) const {
            return Vector<R>(simd_mul(data, rhs.getData()));
        }
        Matrix operator/(double rhs) const {
            return Matrix(simd_mul(1.0 / rhs, data));
        }
        Matrix operator+=(const Matrix<R, C>& rhs) {
            data = simd_add(data, rhs.getData());
            return *this;
        }
        Matrix operator-=(const Matrix<R, C>& rhs) {
            data = simd_sub(data, rhs.getData());
            return *this;
        }
        Matrix operator*=(double rhs) {
            data = simd_mul(rhs, data);
            return *this;
        }
        Matrix operator/=(double rhs) {
            data = simd_mul(1.0 / rhs, data);
            return *this;
        }
        bool operator==(const Matrix<R, C>& rhs) const {
            return simd_almost_equal_elements(data, rhs.getData(), EPS);
        }
        bool operator!=(const Matrix<R, C>& rhs) const {
            return !simd_almost_equal_elements(data, rhs.getData(), EPS);
        }
        Matrix operator-() const {
            return Matrix(simd_mul(-1, data));
        }
        Matrix operator*(const Matrix<C, R>& rhs) {
            return Matrix(matrix_multiply(data, rhs.getData()));
        }
        Matrix Transpose() const {
            return Matrix(simd_transpose(data));
        }
        std::size_t row() const {
            return R;
        }
        std::size_t col() const {
            return C;
        }
        double operator()(std::size_t i, std::size_t j) const {
            return data.columns[j][i];
        }
        double& operator[](std::size_t i, std::size_t j) {
            return data.columns[j][i];
        }
        p_matrix getData() const {
            return data;
        }
    private:
        p_matrix data;
    };

    template<std::size_t N>
    requires (N == 2 || N == 3 || N == 4)
    Matrix<N, N> operator*=(Matrix<N, N>& lhs, Matrix<N, N>& rhs) {
        lhs = Matrix<N, N>(matrix_multiply(lhs.getData(), rhs.getData()));
        return lhs;
    }


    template<std::size_t R>
    requires (R == 2 || R == 3 || R == 4)
    Matrix<R, R> makeIdentity() {
        if constexpr (R == 2) {
            return Matrix<2, 2>(matrix_identity_double2x2);
        } else if constexpr (R == 3) {
            return Matrix<3, 3>(matrix_identity_double3x3);
        } else if constexpr (R == 4) {
            return Matrix<4, 4>(matrix_identity_double4x4);
        }
        else {
            throw std::invalid_argument("Matrix must be 2x2, 3x3, or 4x4");
        }
    }

    inline Matrix<2, 2> makeRotationMatrixR2(double angle_rad) {
        if (angle_rad > M_PI * 2)
            throw std::invalid_argument("angle should be less than 2 PI");
        double c = std::cos(angle_rad);
        double s = std::sin(angle_rad);
        return Matrix<2, 2>(simd::double2x2{simd_double2{c, s}, simd_double2{-s, c}});
    }

    inline std::unordered_map<std::string, Matrix<3, 3>> makeRotationMatrixR3(double angle_phi, double angle_theta, double angle_psi) {
        if(angle_phi > M_PI * 2 || angle_theta > M_PI * 2 || angle_psi > M_PI * 2)
            throw std::invalid_argument("angle should be less than 2 PI");
        double c_phi = std::cos(angle_phi);
        double s_phi = std::sin(angle_phi);
        double c_theta = std::cos(angle_theta);
        double s_theta = std::sin(angle_theta);
        double c_psi = std::cos(angle_psi);
        double s_psi = std::sin(angle_psi);
        auto x = Matrix<3, 3>(simd::double3x3{simd_double3{1.0, 0.0, 0.0}, simd_double3{0.0, c_phi, s_phi}, simd_double3{0.0, -s_phi, c_phi}});
        auto y = Matrix<3, 3>(simd::double3x3{simd_double3{c_theta, 0, -s_theta}, simd_double3{0, 1, 0}, simd_double3{s_theta, 0, c_theta}});
        auto z = Matrix<3, 3>(simd::double3x3{simd_double3{c_psi, s_psi, 0}, simd_double3{-s_psi, c_psi, 0}, simd_double3{0, 0, 1}});
        return {{"x", x}, {"y", y}, {"z", z}};
    }

    template<std::size_t N>
    requires (N == 2 || N == 3)
    Matrix<N, N> makeScaleMatrix(double x, double y) {
        if (x == 0.0 || y == 0.0) throw std::invalid_argument("Scale factor cannot be zero");
        if constexpr (N == 2) return Matrix<2, 2>(simd::double2x2{simd_double2{x, 0.0}, simd_double2{0.0, y}});
        if constexpr (N == 3) return Matrix<3, 3>(simd::double3x3{simd_double3{x, 0.0, 0.0}, simd_double3{0.0, y, 0.0}, simd_double3{0.0, 0.0, 1.0}});
        else throw std::invalid_argument("Matrix must be 2x2 or 3x3");
    }

    Matrix<3, 3> makeTranslationMatrix(double x, double y) {
        return Matrix<3, 3>(simd::double3x3{simd_double3{1.0, 0.0, 0.0}, simd_double3{0.0, 1.0, 0.0}, simd_double3{x, y, 1.0}});
    }

}


#endif //MATH_MATRIX_HPP
