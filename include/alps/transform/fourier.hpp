/*
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_FOURIER_HPP
#define ALPS_TRANSFORM_FOURIER_HPP

#include "fftw.hpp"

namespace fftw = alps::fftw;

namespace alps { namespace transform {

/**
 * Transformer representing a discrete Fourier transform.
 *
 * Performs the unnormalized discrete Fourier transform (FFTW convention):
 *
 *     out[k] = sum[i] exp(direction * 2*pi/n * k * i) * in[i]
 *
 * By default, it uses FFTW when the library is available, and a naive
 * computation otherwise.
 */
class dft
{
public:
    dft() : fftw_() { }

    dft(unsigned n, int direction, bool use_fftw=fftw::SUPPORTED);

    void operator() (std::complex<double> *out, const std::complex<double> *in);

    unsigned in_size() const { return n_; }

    unsigned out_size() const { return n_; }

    bool use_fftw() const { return fftw_.is_initialized(); }

    const fftw::wrapper<> &fftw() const { return fftw_; }

    void naive(std::complex<double> *out, const std::complex<double> *in) const;

protected:
    fftw::wrapper<> &fftw() { return fftw_; }

private:
    fftw::wrapper<> fftw_;
    unsigned n_;
    int direction_;
};

enum statistics
{
    bosonic = 0,
    fermionic = 1
};

/**
 * Transformation from matsubara frequencies to imaginary time.
 */
class iw_to_tau_real
{
public:
    iw_to_tau_real() { }

    iw_to_tau_real(unsigned int niw, unsigned int ntau, double beta,
                   statistics stat, bool use_fftw = fftw::SUPPORTED);

    void operator() (const std::complex<double> *in, double *out);

    void naive(const std::complex<double> *in, double *out) const;

private:
    unsigned niw_, ntau_, oversampling_;
    double beta_;
    statistics stat_;
    fftw::wrapper<> fft_;
};



}} /* alps::transform */

#endif

