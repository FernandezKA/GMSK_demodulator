#include "GMSK_demodulator.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

int main() {

  std::ifstream dump;
  dump.open("/home/fka/dev_dsp/GMSK_demodulator/gmsk_1.txt");
  if (dump.is_open()) {
    std::vector<std::complex<int16_t>> iq;
    std::complex<int16_t> sample;
    while (dump >> sample)
      iq.push_back(sample);
    std::cout << "IQ samples size : " << iq.size() << std::endl;
    Demodulators::GMSK gmsk(4096, 9600, 100000);

    gmsk.Add_Samples(gmsk.LPF(iq));
    auto demodulated = gmsk.Demodulate();
    std::cout << "Demodulated bitstream size is " << demodulated.size()
              << std::endl;

    for(const auto& s : demodulated)
        std::cout<<s;
    std::cout << std::endl; 

  } else {
    std::cout << "Can't open file with I/Q samples! \n\r";
  }

  return 0;
}