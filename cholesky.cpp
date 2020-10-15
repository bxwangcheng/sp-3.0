/**
     * Specializing for "dec" member function.
     */

#include "cholesky.h"
namespace splab {

    template<>
    void Cholesky<complex<float> >::dec(const Matrix<complex<float> > &A) {
        float d;
        complex<float> s;
        L = Matrix<complex<float> >(A.cols(), A.cols());

        MAINLOOP;
    }

    template<>
    void Cholesky<complex<double> >::dec(const Matrix<complex<double> > &A) {
        double d;
        complex<double> s;
        L = Matrix<complex<double> >(A.cols(), A.cols());

        MAINLOOP;
    }

    template<>
    void Cholesky<complex<long double> >::dec(const Matrix<complex<long double> > &A) {
        long double d;
        complex<long double> s;
        L = Matrix<complex<long double> >(A.cols(), A.cols());

        MAINLOOP;
    }


/**
 * Specializing for "solve" member function.
 */
    template<>
    Vector<complex<float> > Cholesky<complex<float> >::solve(const Vector<complex<float> > &b) {
        int n = L.rows();
        if (b.dim() != n)
            return Vector<complex<float> >();
        Vector<complex<float> > x = b;

        SOLVE1
    }

    template<>
    Vector<complex<double> > Cholesky<complex<double> >::solve(const Vector<complex<double> > &b) {
        int n = L.rows();
        if (b.dim() != n)
            return Vector<complex<double> >();
        Vector<complex<double> > x = b;

        SOLVE1
    }

    template<>
    Vector<complex<long double> > Cholesky<complex<long double> >::solve(const Vector<complex<long double> > &b) {
        int n = L.rows();
        if (b.dim() != n)
            return Vector<complex<long double> >();
        Vector<complex<long double> > x = b;

        SOLVE1
    }


/**
 * Specializing for "solve" member function.
 */
    template<>
    Matrix<complex<float> > Cholesky<complex<float> >::solve(const Matrix<complex<float> > &B) {
        int n = L.rows();
        if (B.rows() != n)
            return Matrix<complex<float> >();
        Matrix<complex<float> > X = B;

        SOLVE2
    }

    template<>
    Matrix<complex<double> > Cholesky<complex<double> >::solve(const Matrix<complex<double> > &B) {
        int n = L.rows();
        if (B.rows() != n)
            return Matrix<complex<double> >();
        Matrix<complex<double> > X = B;

        SOLVE2
    }

    template<>
    Matrix<complex<long double> > Cholesky<complex<long double> >::solve(const Matrix<complex<long double> > &B) {
        int n = L.rows();
        if (B.rows() != n)
            return Matrix<complex<long double> >();
        Matrix<complex<long double> > X = B;

        SOLVE2
    }
}