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

class conv_gaussian
{
public:
    conv_gaussian();

    conv_gaussian(unsigned nfreq, unsigned sigma, unsigned half_width,
                  double beta);

    void set_tau(double *tau, unsigned ntau);

    template <typename T>
    void operator() (const T *in, T *out) const;

    unsigned tau_size() const { return precomp_.size(); }

    unsigned conv_size() const { return nfreq_ * sigma_ + 2 * width(); }

    unsigned freq_size() const { return nfreq_; }

    unsigned width() const { return 2 * half_width_; }

private:
    unsigned nfreq_, sigma_, half_width_;
    double beta_;

    double tnorm_, texp_;
    std::vector<double> lfact_;

    struct tau_point {
        unsigned start;
        double baseval;
        double step;
    };
    std::vector<tau_point> precomp_;
};

template <> void conv_gaussian::operator()(const double *, double *) const;
template <> void conv_gaussian::operator()(const std::complex<double>*, std::complex<double>*) const;

} }

#endif /* ALPS_TRANSFORM_NONUNIFORM_HPP */
