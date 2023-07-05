#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace Demodulators {

class GMSK {

  std::vector<std::complex<int16_t>> input_stream;
  std::vector<bool> output_stream;

  std::vector<uint8_t> NRZI_conversion(const std::vector<uint8_t> &) const;
  bool CheckQuadrantSequence(uint8_t last, uint8_t next) const;
  const double phase_threshold = 1;

public:
  GMSK(size_t frame_size = 4096, size_t baudrate = 9600,
       size_t sample_rate = 100000);
  ~GMSK();
  int Add_Samples(std::vector<std::complex<int16_t>> &);
  const std::vector<bool> &Get_Bitstream() const;

  std::vector<bool> Demodulate(void);
  std::vector<std::complex<int16_t>> &LPF(std::vector<std::complex<int16_t>> &);

protected:
  size_t frame_size;
  size_t baudrate;
  size_t sample_rate;
  size_t freq_deviation;
};

}; // namespace Demodulators
