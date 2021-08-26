#include <iostream>

#include "Matrix.hpp"
#include "Polynomial.hpp"

int main() {

    constexpr cte::mat::Matrix non_singular{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
    constexpr cte::mat::Matrix singular{ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 8 } };

    constexpr cte::mat::Matrix copy = non_singular;
    constexpr cte::mat::Matrix add = non_singular + copy;
    constexpr cte::mat::Matrix hadamard = non_singular ^ copy;
    constexpr cte::mat::Matrix product = non_singular * copy;
    constexpr cte::mat::Matrix concat = product.concatenateCols(non_singular);
    constexpr cte::mat::Matrix inverse = singular.inverse();

    std::cout << "The inverse matrix of\n" << singular;
    std::cout << "is the following matrix: \n" << inverse;
    std::cout << "and the value of the determinant is: \n";
    std::cout << singular.determinant() << '\n';

    ignore(add, hadamard, product, concat);

    std::cout << '\n';

    constexpr cte::poly::Polynomial fifth{ 2, 0, 3, -7, 1, -10 };
    constexpr auto roots = fifth.roots();

    std::cout << "Roots of the polynomial " << fifth << " are:\n";
    for (std::size_t i = 0; i < roots.size(); i++) {
        std::cout << roots[i] << '\n';
    }

    std::cout << '\n';

    constexpr cte::poly::Polynomial p1{1,2,3,4};
    constexpr cte::poly::Polynomial p2{2,3,4,1,10,7};
    constexpr cte::poly::Polynomial p11{1,1,2,1,6};
    constexpr cte::poly::Polynomial p22{4,1,2,2};
    constexpr cte::poly::Polynomial p3 = cte::poly::add<p1, p2>();
    constexpr cte::poly::Polynomial p4 = cte::poly::sub<p1, p2>();
    constexpr cte::poly::Polynomial p5 = cte::poly::mult<p1, p2>();
    constexpr cte::poly::Polynomial p6 = cte::poly::div<p11, p22>().getQuotient();
    constexpr cte::poly::Polynomial p7 = cte::poly::div<p11, p22>().getRemainder();
    std::cout << "Given polynomial functions:\n";
    std::cout << "p1: " << p1 << '\n' << "p2: " << p2 << '\n';
    std::cout << "we can do the following operations:\n";
    std::cout << "Adition: " << p3 << '\n';
    std::cout << "Subtraction: " << p4 << '\n';
    std::cout << "Multiplication: " << p5 << '\n';
    std::cout << "Division: Q: " << p6 << '\n';
    std::cout << "          R: " << p7 << '\n';

    std::cout << '\n';

    constexpr cte::mat::Matrix QR_test_matrix = { { 12, -51, 4 }, {6, 167, -68 }, {-4, 24, -41} };
    constexpr auto QR = QR_test_matrix.QRDecomposition();
    constexpr auto Q = QR.getOrthogonal();
    constexpr auto R = QR.getUpperTriangular();
    std::cout << "QR decomposition of the matrix:\n";
    std::cout << QR_test_matrix;
    std::cout << "is the following:\n";
    std::cout << "Q:\n" << Q << "R:\n" << R << '\n';

    constexpr cte::poly::Polynomial e1{-1.0, 1.0};
    constexpr cte::poly::Polynomial e2{2.0};
    constexpr cte::poly::Polynomial e3{3.0};
    constexpr cte::poly::Polynomial e4{4.0};
    constexpr cte::poly::Polynomial e5{-1.0, 5.0};
    constexpr cte::poly::Polynomial e6{6.0};
    constexpr cte::poly::Polynomial e7{7.0};
    constexpr cte::poly::Polynomial e8{8.0};
    constexpr cte::poly::Polynomial e9{-1.0, 8.0};
    
    std::cout << cte::poly::sub<cte::poly::mult<cte::poly::Polynomial{-1.0, 1.0}, cte::poly::Polynomial{-1.0, 5.0}>(), cte::poly::mult<cte::poly::Polynomial{2.0}, cte::poly::Polynomial{4.0}>()>() << '\n';
    std::cout << cte::poly::sub<cte::poly::mult<cte::poly::Polynomial{-1.0, 1.0}, cte::poly::Polynomial{6.0}>(), cte::poly::mult<cte::poly::Polynomial{3.0}, cte::poly::Polynomial{4.0}>()>() << '\n';
    std::cout << cte::poly::sub<cte::poly::mult<cte::poly::Polynomial{-1.0, 1.0}, cte::poly::Polynomial{8.0}>(), cte::poly::mult<cte::poly::Polynomial{2.0}, cte::poly::Polynomial{7.0}>()>() << '\n';
    std::cout << cte::poly::add<cte::poly::Polynomial{-1.0, 1.0}, cte::poly::Polynomial{-1.0, 8.0}>() << '\n';

    return 0;
}

