#include "GMSK_demodulator.hpp"
#include <cstddef>

namespace Demodulators {

GMSK::GMSK(size_t frame_size, size_t baudrate, size_t sample_rate,
           size_t freq_deviation)
    : frame_size(frame_size), baudrate(baudrate), sample_rate(sample_rate),
      freq_deviation(freq_deviation) {
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

size_t GMSK::Normalization(std::vector<std::complex<int16_t>>& samples){
    int16_t bitmask_shift = 0; 

    for(const auto& s : samples){
        bitmask_shift |= s.real; 
        bitmask_shift |= s.imag; 
    }

    for(auto& s : )
}

const std::vector<uint8_t> &GMSK::Get_Bitstream() const {
  return output_stream;
}

} // namespace Demodulators
