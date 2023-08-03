#include "GMSK_demodulator.hpp"
#include <cstddef>
#include <cstdint>
#include <ios>
#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>

int main()
{
  GMSK<std::complex<int16_t>, std::complex<double>> gmsk(100'000, 9600);

  std::ifstream dump;
  dump.open("/home/fka/dev_dsp/GMSK_demodulator/Includes/preamble.bin", std::ios::binary);

  if (dump.is_open())
  {
    std::cout << "File has been open successfully\n\r";
    dump.seekg(0, std::ios::end);
    long unsigned int size = static_cast<long unsigned int>(dump.tellg());
    dump.close();
    dump.open("/home/fka/dev_dsp/GMSK_demodulator/Includes/preamble.bin", std::ios::binary);
    std::vector<std::complex<int16_t>> iq(size / sizeof(std::complex<int16_t>));
    std::cout << "File size " << size << " bytes, " << size / sizeof(std::complex<int16_t>) << " samples. \n\r";
    dump.read(reinterpret_cast<char *>(&iq[0]), size);
    auto start = chrono::steady_clock::now();
    gmsk.add_samples(iq);
    auto packet_index = gmsk.try_to_find_packet();
    auto result = gmsk.demodulate(packet_index);
    auto end = chrono::steady_clock::now();
    std::cout << "Elapsed time [us]: "
              << chrono::duration_cast<chrono::microseconds>(end - start).count()
              << " us" << endl;
    std::cout << "Demodulated bitstream: \n\r";
    size_t index = 1;
    for (const auto &s : result)
    {
      std::cout << s;
    }
    std::cout << endl;
  }
  else
  {
    std::cout << "Can't open file with IQ /n/r";
  }

  return 0;
}