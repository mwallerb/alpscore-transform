%module(package="alps") transform

%include "attribute.i"
%include "std_except.i"
%include "std_complex.i"

%include "numpy.i"
%include "transform.i"

%{
#define SWIG_FILE_WITH_INIT

#include <alps/transform/common.hpp>
#include <alps/transform/fftw.hpp>
#include <alps/transform/fourier.hpp>
#include <alps/transform/model.hpp>
#include <alps/transform/nonuniform.hpp>
%}

/*
 * %init
 * %{
 *     // Ensures that numpy is set up properly
 *     import_array();
 * %}
 */

#pragma SWIG nowarn=320
#define ALPS_NO_WRAPPERS

/* ---------------------------- */

%include <alps/transform/common.hpp>

/* ---------------------------- */

%ignore alps::fftw::alloc;
%ignore alps::fftw::allocator::rebind;
%ignore alps::fftw::wrapper::operator=;
%ignore alps::fftw::swap;

%rename(input) in();
%rename(input) in() const;
%rename(output) out();
%rename(output) out() const;

%include <alps/transform/fftw.hpp>

%template(WrapperCR) alps::fftw::wrapper< std::complex<double>, double >;
%template(WrapperRC) alps::fftw::wrapper< double, std::complex<double> >;
%template(WrapperCC) alps::fftw::wrapper< std::complex<double>, std::complex<double> >;

/* ---------------------------- */

%include <alps/transform/fourier.hpp>

/* ---------------------------- */

%include <alps/transform/nonuniform.hpp>

%template(convolve_tauR) alps::transform::conv_gaussian::operator()<double>;
%template(convolve_tauC) alps::transform::conv_gaussian::operator()< std::complex<double> >;
