%module(package="alps") transform

%include "attribute.i"
%include "std_except.i"
%include "std_complex.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init
%{
    // Ensures that numpy is set up properly
    import_array();
%}

%{
#include <alps/transform/common.hpp>
#include <alps/transform/fftw.hpp>
#include <alps/transform/fourier.hpp>
#include <alps/transform/model.hpp>
#include <alps/transform/nonuniform.hpp>
%}

#pragma SWIG nowarn=320
#define ALPS_NO_WRAPPERS

/* --------------------------- COMMON ---------------------------- */

%include <alps/transform/common.hpp>

%define TRANSFORM(class, in_type, out_type)
    %attribute(alps::transform::class, unsigned, in_size, in_size)
    %attribute(alps::transform::class, unsigned, out_size, out_size)

    %extend alps::transform::class {
        %ignore class();
        %ignore operator() (const in_type *in, out_type *out);

        %apply (in_type *IN_ARRAY1, int DIM1) { (const in_type *in, int nin) };
        %apply (out_type *INPLACE_ARRAY1, int DIM1) { (out_type *out, int nout) };
        %catches(std::runtime_error);
        void operator() (const in_type *in, int nin, out_type *out, int nout)
        {
            if (nin != self->in_size())
                throw std::runtime_error("Invalid input size");
            if (nout != self->out_size())
                throw std::runtime_error("Invalid output size");
            self->operator() (in, out);
        }
    }
%enddef

/* --------------------------- FFTW ---------------------------- */

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

/* --------------------------- FOURIER ---------------------------- */


TRANSFORM(dft, std::complex<double>, std::complex<double>)
TRANSFORM(iw_to_tau_real, std::complex<double>, double)
TRANSFORM(tau_to_iw_real, double, std::complex<double>)

%include <alps/transform/fourier.hpp>

/* --------------------------- NONUNIFORM ---------------------------- */

TRANSFORM(conv_gaussian, double, double)

%include <alps/transform/nonuniform.hpp>


