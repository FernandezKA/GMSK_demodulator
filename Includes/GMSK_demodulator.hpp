#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

template <typename sample_type, typename sample, size_t N> class GMSK {

  const static size_t Fs = 100000;
  const static size_t BAUDRATE = 9600;
  const static size_t sample_per_symbol = 10; 

  size_t size;

  vector<sample_type> stream;
  vector<sample> I, Q;

public:
  GMSK() : size(N) {}

  ~GMSK() {}

  void add_samples(const vector<sample_type> &rhs) {
    stream = rhs;
    separate_complex(stream, I, Q);
  }

  vector<sample> end_diff(const vector<sample> &in) const {
    std::vector<sample> diff;
    diff.resize(in.size());
    size_t size = in.size();
    for (size_t i = 1; i < size; ++i) {
      diff.at(i) = in.at(i) - in.at(i - 1);
    }
    std::cout << "Size is " << diff.size() << " " << in.size() << std::endl;
    assert(diff.size() == in.size());
    return diff;
  }

  vector<sample> vec_sum(const vector<sample> &lhs,
                         const vector<sample> &rhs) const {
    assert(lhs.size() == rhs.size());
    vector<sample> sum(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      sum.at(i) = lhs.at(i) + rhs.at(i);
    }
    return sum;
  }

  vector<sample> vec_sub(const vector<sample> &lhs,
                         const vector<sample> &rhs) const {
    assert(lhs.size() == rhs.size());
    vector<sample> sub(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      sub.at(i) = lhs.at(i) - rhs.at(i);
    }
    return sub;
  }

  vector<sample> vec_mul(const vector<sample> &lhs,
                         const vector<sample> &rhs) const {
    assert(lhs.size() == rhs.size());
    vector<sample> mul(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      mul.at(i) = (lhs.at(i) * rhs.at(i));
    }
    return mul;
  }

  vector<sample> vec_div(const vector<sample> &lhs,
                         const vector<sample> &rhs) const {
    assert(lhs.size() == rhs.size());
    vector<sample> div(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      div.at(i) = lhs.at(i) - rhs.at(i);
    }
    return div;
  }

  vector<sample> vec_pow2(const vector<sample> &lhs) const {
    vector<sample> pow(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      pow.at(i) = (lhs.at(i) * lhs.at(i));
    }
    return pow;
  }

  sample get_average(const vector<sample>& in) const
  {
    int64_t sum = in.at(0); 
    for(const auto& s : in){
      sum+=s;
    }
    sum = sum / static_cast<int64_t>(in.size());
    return static_cast<sample>(sum);
  }

  vector<sample> cum_sum(const vector<sample>& in, size_t win_len) const
  {
    vector<sample> c_sum(in.size()); 
    sample sum = in.at(0); 
    for(size_t i = 0; i < in.size(); ++i){
      if(i % win_len)
        sum += in.at(i);
      else
        sum = in.at(i);
      c_sum.at(i) = sum; 
    }
    return c_sum;
  }

  vector<bool> hard_decision(const vector<sample>& in) const
  {
    vector<bool> bitstream;
    auto cs_freq = cum_sum(in, sample_per_symbol);
    auto avg = get_average(cs_freq);
    std::cout << "Average " << avg << std::endl;
    get_log(cs_freq, "cs_freq.txt");
    for(size_t i = sample_per_symbol/2; i < cs_freq.size(); i += sample_per_symbol){
      if(cs_freq.at(i) > avg){
        bitstream.push_back(1);
      }
      else{
        bitstream.push_back(0);
      }
    }
    return bitstream;
  }

  bool separate_complex(const vector<sample_type> &in, vector<sample> &re,
                        vector<sample> &im) {
    if (in.size() == 0)
      return false;

    re.clear();
    im.clear();

    re.reserve(in.size());
    im.reserve(in.size());
    for (const auto &s : in) {
      re.push_back(s.real());
      im.push_back(s.imag());
    }
    return true;
  }

  vector<bool> demodulate() const {
    auto dt_real = end_diff(I);
    auto dt_imag = end_diff(Q);

    get_log(dt_real, "dt_real.txt");
    get_log(dt_imag, "dt_imag.txt");

    auto I_mul = vec_mul(dt_real, Q);
    auto Q_mul = vec_mul(dt_imag, I);

    get_log(I_mul, "I_mul.txt");
    get_log(Q_mul, "Q_mul.txt");

    auto Sub = vec_sub(Q_mul, I_mul);

    get_log(Sub, "Sub.txt");

    auto pow_re = vec_pow2(I);
    auto pow_im = vec_pow2(Q);

    get_log(pow_re, "pow_re.txt");
    get_log(pow_im, "pow_im.txt");

    auto sum_pow = vec_sum(pow_re, pow_im);

    get_log(sum_pow, "sum_pow.txt");

    auto freq = vec_div(Sub, sum_pow);

    get_log(freq, "freq.txt");

    auto bitstream = hard_decision(freq);

    return bitstream;
  }

  bool get_log(vector<sample> &in, string path) const {
    ofstream dump(path);
    if (dump.is_open()) {
      for (const auto &s : in) {
        dump << s << std::endl;
      }
      dump.close();
    }
    return true;
  }

  bool get_log(vector<sample_type> &in, string path) const {
    ofstream dump(path);
    if (dump.is_open()) {
      for (const auto &s : in) {
        dump << s.real() << ", " << s.imag() << std::endl;
      }
      dump.close();
    }
    return true;
  }
};