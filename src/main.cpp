#include <fmt/format.h>

#include "Matrix.hpp"
#include "Polynomial.hpp"

int main() {

    constexpr cte::mat::Matrix singular{ { 1.0l, 2.0l, 3.0l }, { 4.0l, 5.0l, 6.0l }, { 7.0l, 8.0l, 9.0l } };
    constexpr cte::mat::Matrix non_singular{ { 1.0l, 2.0l, 3.0l }, { 4.0l, 5.0l, 6.0l }, { 7.0l, 8.0l, 8.0l } };

    constexpr cte::mat::Matrix copy = singular;
    constexpr cte::mat::Matrix add = singular + copy;
    constexpr cte::mat::Matrix hadamard = singular ^ copy;
    constexpr cte::mat::Matrix product = singular * copy;
    constexpr cte::mat::Matrix concat = product.concatenateCols(singular);
    constexpr cte::mat::Matrix inverse = non_singular.inverse();

    constexpr cte::mat::Matrix QR_test_matrix = { { 12.0l, -51.0l, 4.0l }, { 6.0l, 167.0l, -68.0l }, { -4.0l, 24.0l, -41.0l } };
    constexpr auto QR = QR_test_matrix.QRDecomposition();
    constexpr auto Q = QR.getOrthogonal();
    constexpr auto R = QR.getUpperTriangular();

    constexpr cte::poly::Polynomial fifth{ 2.0l, 0.0l, 3.0l, -7.0l, 1.0l, -10.0l };
    constexpr auto roots = fifth.roots();

    constexpr cte::poly::Polynomial p1{ 1.0l, 2.0l, 3.0l, 4.0l };
    constexpr cte::poly::Polynomial p2{ 2.0l, 3.0, 4.0l, 1.0l, 10.0l, 7.0l };
    constexpr cte::poly::Polynomial p11{ 1.0l, 1.0l, 2.0l, 1.0l, 6.0l };
    constexpr cte::poly::Polynomial p22{ 4.0l, 1.0l, 2.0l, 2.0l };
    constexpr cte::poly::Polynomial p3 = cte::poly::add<p1, p2>();
    constexpr cte::poly::Polynomial p4 = cte::poly::sub<p1, p2>();
    constexpr cte::poly::Polynomial p5 = cte::poly::mult<p1, p2>();
    constexpr cte::poly::Polynomial p6 = cte::poly::div<p11, p22>().getQuotient();
    constexpr cte::poly::Polynomial p7 = cte::poly::div<p11, p22>().getRemainder();

    constexpr cte::mat::Matrix eigen_test{ { 1.0l, 2.0l, 3.0l, 4.0l }, { 4.0l, 5.0l, 6.0l, 7.0l }, { 7.0l, 8.0l, 8.0l, 9.0l }, {10.0l, 11.0l, 12.0l, 16.0l } };
    constexpr auto eigenvalues = cte::mat::calculateEigenvalues<eigen_test>();
    constexpr auto eigenvectors = cte::mat::calculateEigenvectors<eigen_test>();

    ignore(add, hadamard, product, concat);

    fmt::print("The inverse matrix of\n{}\n", non_singular);
    fmt::print("is the following matrix:\n{}\n", inverse);
    fmt::print("and the value of the determinant is: \n");
    fmt::print("{:.4f}\n\n", non_singular.determinant());

    fmt::print("QR decomposition of the matrix:\n");
    fmt::print("{}\n", QR_test_matrix);
    fmt::print("is the following:\n");
    fmt::print("Q:\n{}\nR:\n{}\n", Q, R);

    fmt::print("Roots of the polynomial {} are:\n", fifth);
    fmt::print("{}\n", roots);

    fmt::print("\nGiven polynomial functions:\n");
    fmt::print("p1: {}\np2: {}\n", p1, p2);
    fmt::print("we can do the following operations:\n");
    fmt::print("Adition: {}\n", p3);
    fmt::print("Subtraction: {}\n", p4);
    fmt::print("Multiplication: {}\n", p5);
    fmt::print("Division: Q: {}\n", p6);
    fmt::print("          R: {}\n\n", p7);

    fmt::print("Eigenvalues of matrix:\n{}\nare:\n", eigen_test);
    fmt::print("{}\n", eigenvalues);
    fmt::print("and eigenvectors are:\n");
    fmt::print("{}\n", eigenvectors);

    return 0;
}

