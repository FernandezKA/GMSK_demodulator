#include "GMSK_demodulator.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

int main() {

  std::ifstream dump;
  dump.open("gmsk_1.txt");
  if (dump.is_open()) {
    std::vector<std::complex<int16_t>> iq;
    std::complex<int16_t> sample;
    while (dump >> sample)
      iq.push_back(sample);
    std::cout << "IQ samples size : " << iq.size() << std::endl;
    Demodulators::GMSK gmsk;

    gmsk.Add_Samples(iq);
    auto demodulated = gmsk.Demodulate();
    std::cout << "Demodulated bitstream size is " << demodulated.size()
              << std::endl;
  } else {
    std::cout << "Can't open file with I/Q samples! \n\r";
  }

  return 0;
}