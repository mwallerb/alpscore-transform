#include <alps/transform/nonuniform.hpp>

namespace alps { namespace transform {

gaussian_window::gaussian_window(unsigned niw, unsigned sigma,
                                 unsigned half_width, double beta)
    : niw_(niw)
    , sigma_(sigma)
    , half_width_(half_width)
    , beta_(beta)
    , lfact_(2 * half_width)
    , zeroval_()
    , step_()
{
    // Runtime checks
    if (sigma == 0)
        throw std::invalid_argument("sigma must be positive");

    // Follows Staar et al (2012), eq. (C.1)
    unsigned m = half_width;
    unsigned n = 4 * niw;
    double b = 2. * sigma/(2 * sigma - 1) * m/M_PI;

    tnorm_ = 1./std::sqrt(M_PI * b);
    texp_ = -1./b;
    knorm_ = 1./n;
    kexp_ = -b * M_PI * M_PI/(1. * n * n);

    // for -m < l <= m, precompute first term of:
    // a*exp[c(l-delta)**2)] = a*exp[c*l**2] exp[c*delta**2] exp[-2*c*delta]**l
    for (unsigned l = 0; l != m + 1; ++l) {
        lfact_[l] = tnorm_ * std::exp(texp_ * l * l);
    }
}

void gaussian_window::set_tau(double *tau, unsigned ntau)
{
}


}}
