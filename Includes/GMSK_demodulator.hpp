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
  std::vector<bool> RemoveInsertions(const std::vector<bool> &) const;

public:
  GMSK();
  ~GMSK();
  int Add_Samples(std::vector<std::complex<double>> &);
  const std::vector<bool> &Get_Bitstream() const;

  std::vector<bool> Demodulate(void);
  std::vector<std::complex<int16_t>> &LPF(std::vector<std::complex<int16_t>> &);
  std::vector<double> MulVec(std::vector<double> &lhs_1, std::vector<double> &lhs_2);
  std::vector<double> CumSum(std::vector<double> &, size_t)const; /*Кумуллятивная сумма, с индексом сброса окна*/
  std::vector<double> FiniteDimensions(const std::vector<double> &in) const;

protected:
  const double Fs = 100000;
  const double Ts = 10.4e-3;
  const size_t baudrate = 9600;
  const size_t samples_per_symbol = static_cast<size_t>(Ts * Fs);
  const size_t frame_size = 0x1000;
};

}; // namespace Demodulators
