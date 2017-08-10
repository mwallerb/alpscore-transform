#include <alps/transform/fourier.hpp>

namespace alps { namespace transform {

dft::dft(unsigned n, int direction, bool use_fftw)
    : fftw_()
    , n_(n)
    , direction_(direction)
{
    if (use_fftw)
        fftw_ = fftw::wrapper<>(fftw::plan_data(n, direction, FFTW_ESTIMATE));
}

void dft::operator() (std::complex<double> *out, const std::complex<double> *in)
{
    if (use_fftw()) {
        // copy to internal buffers
        std::copy(in, in + n_, fftw().in());
        fftw().execute();
        std::copy(fftw().out(), fftw().out() + n_, out);
    } else {
        naive(out, in);
    }
}

void dft::naive(std::complex<double> *out, const std::complex<double> *in) const
{
    // naive implementation of the formula
    //         f_hat[k] = sum_j exp(+2pi i/N k j) f_conv[j]
    double p2piiIn = direction_*2*M_PI/n_;
    for (unsigned k = 0; k != n_; ++k) {
        std::complex<double> f_hatk = 0;
        for (unsigned j = 0; j != n_; ++j) {
            f_hatk += std::exp(std::complex<double>(0, p2piiIn * k * j))
                        * in[j];
        }
        *(out++) = f_hatk;
    }
}




}}
