#include "GMSK_demodulator.hpp"
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>
#include <complex>

int main() {
  GMSK<std::complex<float>, float, 2310> gmsk; 

  std::ifstream dump; 
  dump.open("/home/fka/dev_dsp/GMSK_demodulator/gmsk_1.txt");

  if(dump.is_open()){
    std::vector<std::complex<float>> iq;
    std::complex<float> sample;
    while (dump >> sample){
      iq.push_back(sample);
    }
    std::cout << "IQ samples size : " << iq.size() << std::endl;
    gmsk.add_samples(iq); 
    auto result = gmsk.demodulate();

    for(const auto& s : result){
      std::cout << s << " ";
    }
    std::cout << std::endl;

  }
  else{
    std::cout << "Can't open file with IQ /n/r";
  }
 
  return 0;
}