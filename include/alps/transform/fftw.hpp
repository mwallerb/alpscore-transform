/**
 * Transparently wrap FFTW library, provide dummy fall-back functions.
 *
 * This header transparently wraps the FFTW 3.x interface, providing dummy
 * functions if FFTW is not available.  This allows to replace #ifdef's with
 * regular ifs.  There is no runtime cost, since every single compiler is able
 * to fold the constant and eliminate the dead path, but is still able to
 * perform static analysis on the code.
 *
 * Also provides the following:
 *
 *     alps::fftw::allocator      // STL allocator
 *     alps::fftw::wrapper        // memory-safe wrapper
 *
 * Copyright (C) 1998-2017 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */
#ifndef ALPS_TRANSFORM_FFTW_WRAP_HH_
#define ALPS_TRANSFORM_FFTW_WRAP_HH_

#include <stdexcept>
#include <vector>

struct fftw_not_available : public std::exception
{ };

#if defined(ALPS_HAVE_FFTW) || defined(ALPS_NO_WRAPPERS)

#include <fftw3.h>

#else /* ALPS_HAVE_FFTW */

#include <cstdlib>
#include <complex>

typedef std::pair<double, double> fftw_complex;

// Dummy library for FFTW3
struct fftw_plan_s;
typedef fftw_plan_s *fftw_plan;

// Dummy flags for FFTW3
#define FFTW_MEASURE        1
#define FFTW_ESTIMATE       1
#define FFTW_DESTROY_INPUT  1
#define FFTW_PRESERVE_INPUT 1
#define FFTW_UNALIGNED      2
#define FFTW_FORWARD        -1
#define FFTW_BACKWARD       +1

#define FFTW_DUMMY_ { throw fftw_not_available(); }


static inline fftw_plan fftw_plan_many_dft(int, const int *, int,
                             fftw_complex *, const int *, int, int,
                             fftw_complex *, const int *, int, int,
                             int, unsigned) FFTW_DUMMY_

static inline fftw_plan fftw_plan_many_dft_r2c(int, const int *, int,
                             double *, const int *, int, int,
                             fftw_complex *, const int *, int, int,
                             unsigned) FFTW_DUMMY_

static inline fftw_plan fftw_plan_many_dft_c2r(int, const int *, int,
                             fftw_complex *, const int *, int, int,
                             double *, const int *, int, int,
                             unsigned) FFTW_DUMMY_

static inline void fftw_execute(const fftw_plan) FFTW_DUMMY_

static inline void fftw_execute_dft(const fftw_plan,
                                    fftw_complex *, fftw_complex *) FFTW_DUMMY_

static inline void fftw_execute_dft_r2c(const fftw_plan,
                                        double *, fftw_complex *) FFTW_DUMMY_

static inline void fftw_execute_dft_c2r(const fftw_plan,
                                        fftw_complex *, double *) FFTW_DUMMY_

static inline void fftw_destroy_plan(fftw_plan) FFTW_DUMMY_

static inline double *fftw_alloc_real(size_t) FFTW_DUMMY_

static inline fftw_complex *fftw_alloc_complex(size_t) FFTW_DUMMY_

static inline int fftw_alignment_of(double *) FFTW_DUMMY_

static inline void fftw_free(void *) FFTW_DUMMY_

#undef FFTW_DUMMY_

#endif /* ALPS_HAVE_FFTW */

namespace alps { namespace fftw {

#ifdef ALPS_HAVE_FFTW
static const bool SUPPORTED = true;
#else
static const bool SUPPORTED = false;
#endif

#if defined(ALPS_DEBUG_FFTW) || !defined(NDEBUG)
static const bool DEBUG = true;
#else
static const bool DEBUG = false;
#endif

struct alignment_mismatch : std::exception { };

/** Helper function for memory-aligned allocation from FFTW */
template <typename T>
static inline T *alloc(size_t);

template <>
inline double *alloc(size_t num) { return fftw_alloc_real(num); }

template <>
inline fftw_complex *alloc(size_t num) { return fftw_alloc_complex(num); }

template <>
inline std::complex<double> *alloc(size_t num)
{
    return (std::complex<double> *)fftw_alloc_complex(num);
}

/**
 * Custom STL allocator for FFTW.
 *
 * This allows to use the aligned-memory allocators of FFTW in the STL
 * containers, by, e.g., using `std::vector< double, alps::fftw::allocator<double> >`
 * for a vector.
 */
template <typename T>
struct allocator
    : public std::allocator<T>
{
    typedef T value_type;
    typedef T &reference;
    typedef T *pointer;
    typedef const T &const_reference;
    typedef const T *const_pointer;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    template <typename U> struct rebind { typedef allocator<U> other; };

    allocator() : std::allocator<T>() {}

    allocator(const allocator &other) : std::allocator<T>(other) {}

    template <typename U>
    allocator(const allocator<U> &other) : std::allocator<T>(other) {}

    ~allocator() {}

    /** Aligned allocation from FFTW */
    pointer allocate(size_type num, const void* /*hint*/ = 0)
    {
        return alloc<T>(num);
    }

    /** Deallocation of FFTW aligned array */
    void deallocate(pointer p, size_type /*num*/)
    {
        fftw_free(p);
    }
};

/**
 * Unified structure for holding (parts of) FFTW plan data.
 */
struct plan_data
{
    plan_data() { }

    plan_data(const std::vector<int> &n, int sign, unsigned flags)
        : n(n)
        , sign(sign)
        , flags(flags)
    { }

    plan_data(int n1, int sign, unsigned flags)
        : n(1, n1)
        , sign(sign)
        , flags(flags)
    { }

    int total_size() const
    {
        int size = 1;
        for (unsigned i = 0; i != n.size(); ++i)
            size *= n[i];
        return size;
    }

    std::vector<int> n;
    int sign;
    unsigned flags;
};

/**
 * Wrapper for creating plans
 */
template <typename InT, typename OutT>
fftw_plan create_plan(plan_data data, InT *in, OutT *out);

template <> fftw_plan create_plan(plan_data, std::complex<double> *, std::complex<double> *out);
template <> fftw_plan create_plan(plan_data, std::complex<double> *, double *);
template <> fftw_plan create_plan(plan_data, double *, std::complex<double> *);

/**
 * Wrapper for executing plans
 */
template <typename InT, typename OutT>
inline void execute_plan(const fftw_plan data, InT *in, OutT *out);

template <>
inline void execute_plan(const fftw_plan p, std::complex<double> *in, std::complex<double> *out) {
    fftw_execute_dft(p, (fftw_complex *)in, (fftw_complex *)out);
}

template <>
inline void execute_plan(const fftw_plan p, std::complex<double> *in, double *out) {
    fftw_execute_dft_c2r(p, (fftw_complex *)in, out);
}

template <>
inline void execute_plan(const fftw_plan p, double *in, std::complex<double> *out) {
    fftw_execute_dft_r2c(p, in, (fftw_complex *)out);
}

/**
 * Memory-safe wrapper for a FFTW plan.
 *
 * FFTW makes it fiendishly difficult to acquire memory safety, because there
 * is not API for cloning a plan and there is no portable way to extract the
 * specifics of a plan since its hidden behind a opaque pointer.  Therefore,
 * we must store the properties of a plan and the corresponding arrays together
 * with the plan.
 */
template <typename InT=std::complex<double>, typename OutT=std::complex<double> >
class wrapper
{
public:
    wrapper();

    wrapper(const plan_data &data);

    wrapper(wrapper &other);

#if __cplusplus >= 201103L
    wrapper(wrapper &&other) : wrapper() { swap(*this, other); }
#endif

    wrapper &operator=(wrapper rhs);

    ~wrapper();

    void execute() { fftw_execute(plan_); }

    void execute_on(InT *in, OutT *out);

    const plan_data &data() const { return data_; }

    const fftw_plan &plan() const { return plan_; }

    unsigned size() const { return in_.size(); }

    const InT *in() const { return &in_[0]; }

    InT *in() { return &in_[0]; }

    const OutT *out() const { return &out_[0]; }

    OutT *out() { return &out_[0]; }

    bool is_initialized() const { return plan_ != NULL; }

private:
    friend void swap(wrapper &left, wrapper &right)
    {
        using std::swap;  // ADL
        swap(left.data_, right.data_);
        swap(left.in_, right.in_);
        swap(left.out_, right.out_);
        swap(left.plan_, right.plan_);
    }

    plan_data data_;
    std::vector<InT, allocator<InT> > in_;
    std::vector<OutT, allocator<OutT> > out_;
    fftw_plan plan_;
};

extern template class wrapper< std::complex<double>, double >;
extern template class wrapper< double, std::complex<double> >;
extern template class wrapper< std::complex<double>, std::complex<double> >;

}} /* namespace alps::fftw */

#endif /* ALPS_TRANSFORM_FFTW_WRAP_HH_ */

