/*
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_FOURIER_HPP
#define ALPS_TRANSFORM_FOURIER_HPP

#include <alps/transform/common.hpp>
#include <alps/transform/fftw.hpp>
#include <alps/transform/model.hpp>

namespace alps { namespace transform {

class dft;
class iw_to_tau_real;
class tau_to_iw_real;



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

    dft(unsigned n, int direction, bool use_fftw=alps::fftw::SUPPORTED);

    void operator() (const std::complex<double> *in, std::complex<double> *out);

    unsigned in_size() const { return n_; }

    unsigned out_size() const { return n_; }

    bool use_fftw() const { return fftw_.is_initialized(); }

    const alps::fftw::wrapper<> &fftw() const { return fftw_; }

    void naive(const std::complex<double> *in, std::complex<double> *out) const;

protected:
    alps::fftw::wrapper<> &fftw() { return fftw_; }

private:
    alps::fftw::wrapper<> fftw_;
    unsigned n_;
    int direction_;
};

template <>
struct traits<dft>
{
    typedef std::complex<double> in_type;
    typedef std::complex<double> out_type;
};


/**
 * Transformation from positive Matsubara frequencies to imaginary time.
 *
 * Performs the following transformation (s = 0/1 for bosonic/fermionic):
 *
 *      iw[k] = pi/beta * (2 * k + s)
 *      tau[n] = beta/Ntau * n
 *      f[n] = 2/beta * sum(k,0,Niw-1) exp(i * iw[k] * tau[n]) fhat[k]
 */
class iw_to_tau_real
{
public:
    iw_to_tau_real() { }

    iw_to_tau_real(unsigned int niw, unsigned int ntau, double beta,
                   statistics stat, bool use_fftw=alps::fftw::SUPPORTED);

    void operator() (const std::complex<double> *in, double *out);

    unsigned in_size() const { return niw_; }

    unsigned out_size() const { return ntau_; }

    double beta() const { return beta_; }

    statistics stat() const { return stat_; }

    double tau_value(unsigned i) const { return beta_/ntau_ * i; }

    void naive(const std::complex<double> *in, double *out) const;

private:
    unsigned niw_, ntau_, oversampling_;
    double beta_;
    statistics stat_;
    alps::fftw::wrapper<> fft_;
};

template<>
struct traits<iw_to_tau_real>
{
    typedef std::complex<double> in_type;
    typedef double out_type;
};


/**
 * Transformation from Matsubara frequencies to imaginary time using a model
 */
class iw_to_tau_model_real
{
    iw_to_tau_model_real() { }

    iw_to_tau_model_real(unsigned int niw, unsigned int ntau, double beta,
                         statistics stat, const std::vector<double> &moments,
                         bool use_fftw=alps::fftw::SUPPORTED);

    void operator() (const std::complex<double> *in, double *out);

    unsigned in_size() const { return transform_.in_size(); }

    unsigned out_size() const { return transform_.out_size(); }

    double beta() const { return transform_.beta(); }

    statistics stat() const { return transform_.stat(); }

private:
    iw_to_tau_real transform_;
    moments_model<double> model_;
    std::vector< std::complex<double> > in_buffer_;
};

template<>
struct traits<iw_to_tau_model_real>
{
    typedef std::complex<double> in_type;
    typedef double out_type;
};

/**
 * Transformation from imaginary time to positive Matsubara frequencies.
 *
 * Performs the following transformation (s = 0/1 for bosonic/fermionic):
 *
 *     tau[n] = beta/Ntau * n
 *     iw[k] = pi/beta * (2 * k + s)
 *     f[k] = beta/N * sum(k,0,Ntau-1) exp(i * iw[k] * tau[n]) f[n]
 */
class tau_to_iw_real
{
public:
    tau_to_iw_real() { }

    tau_to_iw_real(unsigned int ntau, unsigned int niw, double beta,
                   statistics stat, bool use_fftw=alps::fftw::SUPPORTED);

    void operator() (const double *in, std::complex<double> *out);

    unsigned in_size() const { return ntau_; }

    unsigned out_size() const { return niw_; }

    double beta() const { return beta_; }

    statistics stat() const { return stat_; }

    void naive(const double *in, std::complex<double> *out) const;

private:
    unsigned niw_, ntau_, oversampling_;
    double beta_;
    statistics stat_;
    alps::fftw::wrapper<> fft_;
};

template<>
struct traits<tau_to_iw_real>
{
    typedef double in_type;
    typedef std::complex<double> out_type;
};

}} /* alps::transform */

#endif

