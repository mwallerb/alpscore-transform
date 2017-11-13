#include <alps/transform/model.hpp>

namespace alps { namespace transform {

double euler_p(int n, double x)
{
    switch (n) {
    case 0:
        return 1;
    case 1:
        return x - 1/2.;
    case 2:
        return x*x - x;
    default:
        throw std::runtime_error("Euler P not implemented for this case");
    }
}

// trivial model without contribution
template <typename T>
moments_model<T>::moments_model()
    : mom_()
    , beta_(1.0)
    , stat_(fermionic)
{ }

template <typename T>
moments_model<T>::moments_model(std::vector<T> mom, double beta, statistics stat)
    : mom_(mom)
    , beta_(beta)
    , stat_(stat)
{ }

template <typename T>
typename moments_model<T>::ComplexT moments_model<T>::fiw(int k) const
{
    double freq = M_PI/beta_ * (2 * k + int(stat_));
    if (freq == 0)
        throw std::runtime_error("Moments model undefined for bosonic freq 0.");

    // mom[0]/iw + mom[1]/(iw**2) + mom[2]/(iw**3) + ...
    // = 1/iw (mom[0] + 1/iw (mom[1] + 1/iw * (mom[2] + ... )
    ComplexT one_over_iw = 1.0 / ComplexT(0, freq);
    ComplexT result = 0;
    for (typename std::vector<T>::const_reverse_iterator it = mom_.rbegin();
            it != mom_.rend(); ++it) {
        result += *it;
        result *= one_over_iw;
    }
    return result;
}

template <typename T>
T moments_model<T>::ftau(double tau) const
{
    if (tau < -beta_ || tau > beta_)
        throw std::runtime_error("Other tau in [-beta, beta]");

    if (stat_ != fermionic)
        throw std::runtime_error("Fourier transform does not exist");

    // Abramowitz & Stegun [23.1.17 & 23.1.18]
    // TODO: possibly unstable for large moment orders
    T result = 0;
    T prefactor = 0.5;
    for (unsigned i = 0; i != mom_.size(); ++i) {
        if (i != 0)
            prefactor *= beta_/i;     // beta^n / 2 n!
        result += mom_[i] * euler_p(i, -tau/beta_);
    }
    return result;
}

template class moments_model<double>;
template class moments_model<std::complex<double> >;


}}

