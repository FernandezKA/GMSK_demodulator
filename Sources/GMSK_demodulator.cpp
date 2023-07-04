#include "GMSK_demodulator.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iomanip>

namespace Demodulators {

GMSK::GMSK(size_t frame_size, size_t baudrate, size_t sample_rate)
    : frame_size(frame_size), baudrate(baudrate), sample_rate(sample_rate) {
  std::cout << "GMSK demodulator is created \n\r";
  input_stream.reserve(frame_size);
}

GMSK::~GMSK() { std::cout << "GMSK demodulator is removed \n\r"; }

int GMSK::Add_Samples(std::vector<std::complex<int16_t>> &input) {
  input_stream.clear();
  for (const auto &s : input)
    input_stream.push_back(s);
  return (input_stream.size() == input.size()) ? (0) : (-1);
}

size_t GMSK::Normalization(std::vector<std::complex<int16_t>> &samples) {
  int16_t bitmask_shift = 0;

  for (const auto &s : samples) {
    bitmask_shift |= abs(s.real());
    bitmask_shift |= abs(s.imag());
  }

  uint32_t ones = (1u << 31);
  uint8_t shift_val = 0;

  while ((bitmask_shift & ones) == ones) {
    shift_val++;
  }

  for (auto &s : samples) {
    std::complex<int16_t> shifted =
        (s.real() << shift_val, s.imag() << shift_val);
    s = shifted;
  }
  return static_cast<size_t>(shift_val);
}

const std::vector<uint8_t> &GMSK::Get_Bitstream() const {
  return output_stream;
}

std::vector<bool> GMSK::Demodulate(void) {
  std::vector<int> demodulatedBits;
  int16_t phase = 0.0;

  for (int i = 0; i < input_stream.size(); i++) {
    std::complex<int16_t> sample = input_stream.at(i);
    double curr_phase = std::atan(sample.imag() / sample.real());
    int bit = (curr_phase - phase > M_PI) ? 1 : 0;
    demodulatedBits.push_back(bit);
    phase = curr_phase;
  }

  size_t oversampling =
      sample_rate / baudrate; /*Количество отсчетов на один символ*/

  std::vector<bool> output_stream;

  for (auto &s : demodulatedBits) {
    static size_t counter = 0;
    static bool last = demodulatedBits.at(0);

    if (s ==
        1) { /*Если в рамках одного временного интервала символа мы
                обнаружили сдвиг фазы - значит это была посылка с единицей*/
      counter = 0;
      output_stream.push_back(true);
      continue;
    }

    if (counter ==
        oversampling) { /*Если за один символьный интервал скачка
                           фазы не было, то это была посылка с нулем*/
      counter = 0;
      output_stream.push_back(false);
      continue;
    }
    ++counter;
  }

  return output_stream;
}

} // namespace Demodulators
