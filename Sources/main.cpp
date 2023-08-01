#include "GMSK_demodulator.hpp"
#include <cstddef>
#include <cstdint>
#include <ios>
#include <iostream>
#include <vector>
#include <complex>

int main()
{
  GMSK<std::complex<int16_t>, std::complex<double>> gmsk(100'000, 9600);

  std::ifstream dump;
  dump.open("/home/fka/dev/GMSK_demodulator/gmsk_1.bin", std::ios::binary);

  if (dump.is_open())
  {
    dump.seekg(0, std::ios::end);
    long unsigned int size = static_cast<long unsigned int>(dump.tellg());
    std::vector<std::complex<int16_t>> iq(size);
    std::cout << "File size " << size << " bytes, " << size / sizeof(std::complex<int16_t>) << " samples. /n/r";
    dump.read(reinterpret_cast<char *>(&iq[0]), static_cast<std::streamsize>(iq.size() * sizeof(std::complex<int16_t>)));
    gmsk.add_samples(iq);
    gmsk.demodulate_without_noise();
    auto result = gmsk.demodulate();

    for (const auto &s : result)
    {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
  else
  {
    std::cout << "Can't open file with IQ /n/r";
  }

  return 0;
}