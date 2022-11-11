#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "excerpt/inc/excerpt.h"

using namespace std;

template<typename fp_t>
void unified_method(vector<fp_t> coef, vector<fp_t> &roots){
    fp_t a = coef[3];
    fp_t b = coef[2];
    fp_t c = coef[1];
    fp_t d = coef[0];
    cout << b << " " << c << " " << d << endl;
    // Коэффициенты, которые в вещественном поле
    fp_t f_1 = fma(c, b, (fp_t)-9 * d)/fma(b, b, (fp_t)-3 * c);
    fp_t f_2 = fma(c, c, (fp_t)-3 * d * b)/fma(b, b, (fp_t)-3 * c);

    cout << endl << f_1 << endl << f_2 << endl;
    // Коэффициенты, которые в комплексном поле
    complex<fp_t> image = sqrt(fma(f_1, f_1, -4*f_2));

    complex<fp_t> b_0 = f_1/(fp_t)2 + image/(fp_t)2;
    complex<fp_t> c_0 = f_1/(fp_t)2 - image/(fp_t)2;
    cout << endl << b_0 << endl << c_0 << endl;

    complex<fp_t> p = pow(((fp_t)2 * b - (fp_t)3 * f_1 - (fp_t)3 * image)/((fp_t)2 * b - (fp_t)3 * f_1 + (fp_t)3 * image), (fp_t)1/3);

    cout << endl << p << endl;
    // Корни уравнения
    roots[0] = ((p * c_0 - b_0)/((complex<fp_t>)1 - p)).real();

    complex<fp_t> coef_p = (fp_t)2 * b_0 + p * ((fp_t)2 * p * c_0 + f_1);
    complex<fp_t> denominator_p = (fp_t)1 + p * ((fp_t)1 + p);
    complex<fp_t> numerator_p = pow(pow(coef_p, 2) - (fp_t)4 * (denominator_p) * (pow(b_0, 2) + p * (p*pow(c_0,2)+f_2)),
                                    (fp_t)1/2);

    roots[1] = ((-coef_p + numerator_p) / ((fp_t)2 * denominator_p)).real();
    roots[2] = ((-coef_p - numerator_p) / ((fp_t)2 * denominator_p)).real();

    cout << endl << roots[0] << " " << roots[1] << " " << roots[2] << endl;
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t deviation;
    vector<fp_t> roots_computed(roots_count);
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, std::numeric_limits<fp_t>::min(), -1, 1, roots, coefficients);
    int cnt_real_roots = 1;
    unified_method(coefficients, roots_computed);
    if (cnt_real_roots!=0) {
        auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, deviation);
        switch (result) {
            case PR_2_INFINITE_ROOTS:
                cout << "INFINITE ROOTS";
                break;
            case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
                cout << "AT LEAST ONE ROOT IS FAKE";
                break;
            case PR_AT_LEAST_ONE_ROOT_LOST:
                cout << "AT LEAST ONE ROOT LOST";
                break;
            default:
                break;
        }

    }
    else deviation = std::numeric_limits<fp_t>::infinity();
    return deviation;
}

int main() {

    float deviation, max_deviation = 0;
    for (auto i = 0; i < 100; ++i) {
        deviation = testPolynomial<float>(3);
        if (deviation != std::numeric_limits<float>::infinity()) {
            cout << "deviation = " << deviation << endl;
            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }
        else cout << "\t\tComplex roots!"<<endl;
    }
    cout<< endl<<"MAX_deviation = "<< max_deviation<<endl;

}