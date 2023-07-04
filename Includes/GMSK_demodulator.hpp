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

GMSK(size_t frame_size, size_t baudrate, size_t sample_rate, size_t freq_deviation); 
~GMSK(); 
int Add_Samples(std::vector<std::complex<int16_t>>&);
const std::vector<uint8_t>& Get_Bitstream() const;  

size_t Normalization(std::vector<std::complex<int16_t>>& );

protected:

size_t frame_size; 
size_t baudrate; 
size_t sample_rate; 
size_t freq_deviation; 

};

}; // namespace Demodulators
