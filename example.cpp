#include "MTLS.h"

double evaluate(const Eigen::VectorXd& x) {
    /* calculate x**2 */
    return x.squaredNorm();
} // end function evaluate

Eigen::VectorXd derivative(const Eigen::VectorXd& x) {
    return x;
} // end function derivative

double evalRaw(const double* x) {
    /* calculate x**2 */
    double norm = 0.0;
    for (unsigned int i = 0; i < sizeof(x)/sizeof(x[0]); ++i) {
        norm += x[i]*x[i];
    } // end fori
    return norm;
} // end function evalRaw

double* derRaw(const double* x) {
    /* calculate x**2 */
    return (double*)x;
} // end function evalRaw
            
class Dummy {
    /* Dummy class for containing calculation function func and
     * derivative derFunc */
    public:
        Dummy(int size) {
            derivative = Eigen::VectorXd::Zero(size);
        };
        virtual ~Dummy() {};

        Eigen::VectorXd derivative;

        double f(const Eigen::VectorXd& x) {
            derivative = der(x);
            return x.squaredNorm();
        } // end function eva

        Eigen::VectorXd der(const Eigen::VectorXd& x) {
            return x;
        }

        const Eigen::VectorXd& g() const {
            return derivative;
        } // end function der
};

int main() {
    /* main function */
    struct Params {
        /* struct of default parameters */
        double maxIterations = 100; // maximum number of iterations
        double mu = 0.5; // step scaling factor, (0<mu<=1/2<eta)
        double eta = 1.0; // termination parameter (0<eta<1)
        double delta = 4.0; // delta: scaling of step updating ([1.1,4.0])
        double bisectWidth = 0.66; // extrapolation tolerance for bounds
        double bracketTol = 1e-14; // termination tolerance for brackets
        double aMin0 = 0.0; // lower bound for step (aMin>=0 and aMin<aMax)
        double aMax0 = 100.0; // upper bound for step (aMax>aMin)
    } params; // end struct Params
    double step;

    // silly initial values
    Eigen::VectorXd x0 = Eigen::VectorXd::Constant(2, 1.0);
    double f0 = evaluate(x0);
    Eigen::VectorXd p = -derivative(x0);
    
    // simple call with default parameters
    step = MTLS::linesearchMoreThuente(p, x0, f0, &evaluate, &derivative);

    // call with struct
    step = MTLS::linesearchMoreThuente(&params, p, x0, f0, &evaluate,
            &derivative);

    // call with class
    Dummy* d = new Dummy(2);
    step = MTLS::linesearchMoreThuente(p, x0, f0, d, &Dummy::f, &Dummy::g);
    step = MTLS::linesearchMoreThuente(&params, p, x0, f0, d, &Dummy::f,
            &Dummy::g);
    delete d;

    // initialize with plain arrays
    double x0Raw[2];
    x0Raw[0] = 2.0;
    x0Raw[1] = 2.0;
    f0 = evalRaw(x0Raw);
    double pRaw[2];
    memcpy(&pRaw, derRaw(x0Raw), sizeof(pRaw));

    // call with plain arrays and default parameters
    step = MTLS::linesearchMoreThuente(pRaw, x0Raw, f0, &evalRaw, &derRaw);
    
    // call with plain arrays and parameters struct
    step = MTLS::linesearchMoreThuente(&params, pRaw, x0Raw, f0, &evalRaw,
            &derRaw);

    return 0;
} // end main
