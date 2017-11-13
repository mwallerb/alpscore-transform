#include <alps/transform/fourier.hpp>

#include "gtest_util.hpp"

template <typename Tf>
class consistency_case
{
public:
    typedef typename alps::transform::traits<Tf>::in_type in_type;
    typedef typename alps::transform::traits<Tf>::out_type out_type;

    consistency_case(const Tf &tf)
    {
        tf_ = tf;
        in_buffer_.resize(tf.in_size());
        out_buffer_.resize(tf.out_size());
        naive_buffer_.resize(tf.out_size());

        for (size_t i = 0; i != tf.in_size(); ++i)
            in_buffer_[i] = i;
    }

    void test_zero_out()
    {
        // see if same result
        for (size_t i = 0; i != tf_.in_size(); ++i) {
            out_buffer_[i] = 0;
            naive_buffer_[i] = 0;
        }
        tf_(&in_buffer_[0], &out_buffer_[0]);
        tf_.naive(&in_buffer_[0], &out_buffer_[0]);
        check();
    }

    void test_prefilled()
    {
        // pre-filled buffer
        for (size_t i = 0; i != tf_.in_size(); ++i) {
            out_buffer_[i] = 1.0 - i;
            naive_buffer_[i] = 1.0 - i;
        }
        tf_(&in_buffer_[0], &out_buffer_[0]);
        tf_.naive(&in_buffer_[0], &out_buffer_[0]);
        check();
    }

    void check()
    {
        EXPECT_EQ(in_buffer_.size(), tf_.in_size());
        EXPECT_EQ(out_buffer_.size(), tf_.out_size());

        for (size_t i = 0; i != out_buffer_.size(); ++i)
            EXPECT_NEAR(naive_buffer_[i], out_buffer_[i], 1e-6);
    }

private:
    Tf tf_;
    std::vector<in_type> in_buffer_;
    std::vector<out_type> out_buffer_;
    std::vector<out_type> naive_buffer_;
};

template <typename T>
consistency_case<T> make_consistency_case(const T &t)
{
    return consistency_case<T>(t);
}

TEST(consistency, test_dft_zero)
{
    make_consistency_case(alps::transform::dft(40, 1)).test_zero_out();
}

TEST(consistency, test_dft_pre)
{
    make_consistency_case(alps::transform::dft(40, 1)).test_prefilled();
}

TEST(consistency, test_idft_zero)
{
    make_consistency_case(alps::transform::dft(26, -1)).test_zero_out();
}

TEST(consistency, test_idft_pre)
{
    make_consistency_case(alps::transform::dft(26, -1)).test_prefilled();
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
