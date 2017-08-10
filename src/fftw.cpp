#include <alps/transform/fftw.hpp>

namespace alps { namespace fftw {

template <>
fftw_plan create_plan(plan_data data, std::complex<double> *in, std::complex<double> *out)
{
    return fftw_plan_many_dft(data.n.size(), &data.n[0], 1,
                              (fftw_complex *)in, NULL, 1, 0,
                              (fftw_complex *)out, NULL, 1, 0,
                              data.sign, data.flags);
}

template <>
fftw_plan create_plan(plan_data data, std::complex<double> *in, double *out)
{
    return fftw_plan_many_dft_c2r(data.n.size(), &data.n[0], 1,
                                  (fftw_complex *)in, NULL, 1, 0,
                                  out, NULL, 1, 0,
                                  data.flags);
}

template <>
fftw_plan create_plan(plan_data data, double *in, std::complex<double> *out)
{
    return fftw_plan_many_dft_r2c(data.n.size(), &data.n[0], 1,
                                  in, NULL, 1, 0,
                                  (fftw_complex *)out, NULL, 1, 0,
                                  data.flags);
}

template <typename InT, typename OutT>
wrapper<InT, OutT>::wrapper()
    : plan_(NULL)
{ }

template <typename InT, typename OutT>
wrapper<InT, OutT>::wrapper(const plan_data &data)
    : data_(data)
    , in_(data.total_size())
    , out_(in_.size())
    , plan_(create_plan(data_, &in_[0], &out_[0]))
{ }

template <typename InT, typename OutT>
wrapper<InT, OutT>::wrapper(wrapper &other)
    : data_(other.data_)
    , in_(other.in_)
    , out_(other.out_)
    , plan_(other.plan_ != NULL ? create_plan(data_, &in_[0], &out_[0]) : NULL)
{ }

template <typename InT, typename OutT>
wrapper<InT, OutT> &wrapper<InT, OutT>::operator=(wrapper rhs)
{
    swap(*this, rhs);
    return *this;
}

template <typename InT, typename OutT>
wrapper<InT, OutT>::~wrapper()
{
    if (plan_ != NULL)
        fftw_destroy_plan(plan_);
}

template class wrapper< std::complex<double>, double >;
template class wrapper< double, std::complex<double> >;
template class wrapper< std::complex<double>, std::complex<double> >;

}} /* namespace alps::fftw */
