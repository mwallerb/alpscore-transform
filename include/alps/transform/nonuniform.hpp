/*
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_NONUNIFORM_HPP
#define ALPS_TRANSFORM_NONUNIFORM_HPP

#include <alps/transform/common.hpp>
#include <alps/transform/fftw.hpp>

#include <vector>
#include <complex>

namespace alps { namespace transform {

class gaussian_window
{
public:
    gaussian_window();

    gaussian_window(unsigned nfreq, unsigned sigma, unsigned half_width,
                    double beta);

    void set_tau(double *tau, unsigned ntau);

    template <typename T>
    void convolve_tau(const T *in, T *out) const;

    void deconvolve_freq(const std::complex<double> *in, std::complex<double> *out);

    unsigned tau_size() const { return precomp_.size(); }

    unsigned conv_size() const { return nfreq_ * sigma_ + 4 * half_width_; }

    unsigned freq_size() const { return nfreq_; }

    template <typename T>
    void convolve_fx_naive(const T *in, T *out);

private:
    unsigned nfreq_, sigma_, half_width_;
    double beta_;

    double tnorm_, texp_, knorm_, kexp_;
    std::vector<double> lfact_;

    struct tau_point {
        unsigned start;
        double baseval;
        double step;
    };
    std::vector<tau_point> precomp_;
};

template <> void gaussian_window::convolve_tau(const double *, double *) const;
template <> void gaussian_window::convolve_tau(const std::complex<double>*, std::complex<double>*) const;

} }

#endif /* ALPS_TRANSFORM_NONUNIFORM_HPP */
