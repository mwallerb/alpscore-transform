/*
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_COMMON_HPP
#define ALPS_TRANSFORM_COMMON_HPP

#include <vector>
#include <stdexcept>

namespace alps { namespace transform {

/**
 * Common traits struct.
 *
 * You should specialize this template for every transformer in such a way:
 *
 *     template <> struct transform_traits<my_transformer>
 *     {
 *         typedef double in_type;                  // input type
 *         typedef double out_type;                 // output type
 *     };
 */
template <typename T>
struct transform_traits;

/**
 * Convenience function that transforms vector
 */
template <typename T>
std::vector<typename transform_traits<T>::out_type> apply(
                    T &tf,
                    const std::vector<typename transform_traits<T>::in_type> &in
                    )
{
    if (tf.in_size() != in.size())
        throw std::runtime_error("size mismatch");

    typedef typename transform_traits<T>::out_type out_type;
    std::vector<out_type> out(tf.out_size(), out_type(0));
    tf(&out[0], &in[0]);
    return out;
}



}}

#endif
