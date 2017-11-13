/*
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_MODEL_HPP
#define ALPS_TRANSFORM_MODEL_HPP

#include <vector>

#include "common.hpp"

namespace alps { namespace transform {

/** Complexifies a type */
template <typename T>
struct make_complex
{
    typedef std::complex<T> type;
    typedef T real_type;
};

template <typename T>
struct make_complex< std::complex<T> >
{
    typedef std::complex<T> type;
    typedef T real_type;
};


template <typename T>
struct model_traits;

template <typename T>
class moments_model
{
public:
    typedef typename make_complex<T>::type ComplexT;

public:
    moments_model();

    moments_model(std::vector<T> mom, double beta, statistics stat);

    ComplexT fiw(int k) const;

    T ftau(double tau) const;

private:
    std::vector<T> mom_;
    double beta_;
    statistics stat_;
};

extern template class moments_model<double>;
extern template class moments_model<std::complex<double> >;


} }

#endif /* ALPS_TRANSFORM_MODEL_HPP */
