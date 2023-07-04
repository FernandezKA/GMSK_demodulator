#pragma once

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string> 
#include <complex>
#include <vector> 

namespace Demodulators {

class GMSK {

std::vector<std::complex<int16_t>> input_stream; 
std::vector<uint8_t> output_stream; 

std::vector<uint8_t> NRZI_conversion(const std::vector<uint8_t>& ) const; 


public: 

GMSK(size_t frame_size = 4096, size_t baudrate = 9600, size_t sample_rate = 100000); 
~GMSK(); 
int Add_Samples(std::vector<std::complex<int16_t>>&);
const std::vector<uint8_t>& Get_Bitstream() const;  

size_t Normalization(std::vector<std::complex<int16_t>>& );
std::vector<bool> Demodulate(void); 

protected:

size_t frame_size; 
size_t baudrate; 
size_t sample_rate; 
size_t freq_deviation; 

};

}; // namespace Demodulators
