#include <alps/transform/nonuniform.hpp>

#include <cassert>

namespace alps { namespace transform {

conv_gaussian::conv_gaussian(unsigned nfreq, unsigned sigma,
                             unsigned half_width, double beta)
    : nfreq_(nfreq)
    , sigma_(sigma)
    , half_width_(half_width)
    , beta_(beta)
    , lfact_(2 * half_width)
    , precomp_()
{
    // Runtime checks
    if (sigma == 0)
        throw std::invalid_argument("sigma must be positive");

    // Follows Staar et al (2012), eq. (C.1)
    unsigned m = half_width;
    double b = 2. * sigma/(2 * sigma - 1) * m/M_PI;

    tnorm_ = 1./std::sqrt(M_PI * b);
    texp_ = -1./b;

    //unsigned n = 4 * nfreq;
    //knorm_ = 1./n;
    //kexp_ = -b * M_PI * M_PI/(1. * n * n);

    // for -m < l <= m, precompute first term of:
    // a*exp[c(l-delta)**2)] = a*exp[c*l**2] exp[c*delta**2] exp[-2*c*delta]**l
    for (int l = -int(m) + 1; l != int(m) + 1; ++l)
        lfact_[l + m - 1] = tnorm_ * std::exp(texp_ * l * l);
}

void conv_gaussian::set_tau(double *tau, unsigned ntau)
{
    precomp_.resize(ntau);

    for (unsigned ix = 0; ix != ntau; ++ix) {
        double tauval = tau[ix];
        assert(tauval >= -beta_ && tauval < beta_);

        // Map the interval -0.5 .. 0.5 to 0 .. 1 of the grid.  Keiner et al.'s
        // convention is to wrap the negative part over to the positive one.
        // Note that this causes a "disconnect": the positive part comes first.
        if (tauval < 0)
            tauval += beta_;

        // Sort of implements Staar, Alg.3, formula 1.  We do not do the
        // truncation of coefficients like they do, since I think it is wrong
        //
        // For pos being an exact integer we would actually have 2n+1 points, but
        // to avoid an if we conceptually replace t -> t - epsilon, thereby
        // eliminating the additional point (it lies epsilon outside [-m,m]).
        double pos = (nfreq_ * sigma_)/(2 * beta_) * tauval + half_width_;
        unsigned start = pos;
        double delta = pos - start;
        assert(start < conv_size() - 2 * half_width_);

        // precompute second as well third term of:
        // a*exp[c*(l-delta)**2)]
        //   = a*exp[c*l**2] exp[c*delta**2] exp[-2*c*delta]**l
        //   = a*exp[...] exp[c*delta*(delta+2*m-2)] exp[-2*c*delta]**(l+m-1)
        double baseval = std::exp(texp_ * delta * (delta + 2*half_width_ - 2));
        double step = std::exp(-2. * texp_ * delta);
        tau_point curr = {start, baseval, step};
        precomp_[ix] = curr;
    }
}

template <typename T>
void conv_gaussian::operator() (const T *in, T *out) const
{
    // now do: buffer[l] += phi(l-delta) * fx for -m < l <= m
    for (unsigned i = 0; i != precomp_.size(); ++i) {
        T *val = in[i] * precomp_[i].baseval;
        unsigned start = precomp_[i].start;
        for (unsigned p = start; p != start + 2 * half_width_; ++p) {
            out[p] += val;
            val *= precomp_[i].step;
        }
    }
}

}}
