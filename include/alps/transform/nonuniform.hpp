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

template <typename T>
class ntau_to_gauss
{
public:
    ntau_to_gauss(unsigned niw, unsigned sigma, unsigned width);

    void set_tau(double *tau, unsigned ntau);

    void operator() (const T *in, T *out);

    unsigned in_size() const { return zeroval_.size(); }

    unsigned out_size() const { return 2 * (niw_ * sigma_ + width_); }

private:
    unsigned niw_, sigma_, width_;
    std::vector<double> lfact_;

    std::vector<T> zeroval_;
    std::vector<double> step_;
};

extern template class ntau_to_gauss<double>;
extern template class ntau_to_gauss< std::complex<double> >;

} }

#endif /* ALPS_TRANSFORM_NONUNIFORM_HPP */
