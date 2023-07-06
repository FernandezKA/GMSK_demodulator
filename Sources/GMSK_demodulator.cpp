/**
 * @file GMSK_demodulator
 * @author FernandezKA
 * @date 2023-07-05
 * @brief Реализация GMSK демодулятора
 * @details
 *  1. Пропускаем I/Q через ФНЧ, таким образом убирая ВЧ;
 *  2. Находим значение фазы для каждого отсчета;
 *  3. Оцениваем разность фаз в точках семплирования, по этому значению смотрим, был ли скачок фазы.
 *  4. Преобразуем из NRZI в битовый поток.
 * @TODO Нужно обеспечить синхронизацию сообщения, чтобы не декодировать поток
 * без сигнала (шумы).
 */

#include "GMSK_demodulator.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <math.h>

namespace Demodulators {

GMSK::GMSK()
{
  std::cout << "GMSK demodulator is created \n\r";
  input_stream.reserve(frame_size);
}

GMSK::~GMSK() { std::cout << "GMSK demodulator is removed \n\r"; }

int GMSK::Add_Samples(std::vector<std::complex<double>> &input) {
  input_stream.clear();
  for (const auto &s : input)
    input_stream.push_back(s);
  return (input_stream.size() == input.size()) ? (0) : (-1);
}

const std::vector<bool> &GMSK::Get_Bitstream() const { return output_stream; }

/**
 * @brief Демодуляция происходит здесь. Описание приложено в начале файла.
 */
std::vector<bool> GMSK::Demodulate(void) {
  std::vector<double> phase_samples;
  for (auto &s : input_stream) {
    phase_samples.push_back(abs(std::atan2(s.imag(), s.real())));
  }

  auto cum_sum = CumSum(phase_samples, samples_per_symbol);  
  std::ofstream phase("phase.txt");
  for(auto s : cum_sum){
      phase << s << std::endl; 
  }

  for(size_t i = samples_per_symbol - 1; i < phase_samples.size(); i+=samples_per_symbol){
    bool bit = (cum_sum.at(i) > static_cast<double>(samples_per_symbol));  
    output_stream.push_back(bit);
  }
  std::cout << "Phase samples size " << phase_samples.size() << std::endl;

  return output_stream;
}

/*Простая имплементация ФНЧ, которая необходима перед фильтрацией*/
std::vector<std::complex<int16_t>> &
GMSK::LPF(std::vector<std::complex<int16_t>> &in) {
  /*TODO: Need to  write LPF filter! */
  return in; 
}


std::vector<double> GMSK::CumSum(std::vector<double>& in, size_t reset_index) const{
  std::vector<double> cum_sum;
  cum_sum.reserve(in.size()); 
  for(size_t i = 0; i < in.size(); ++i){
    static double sum = 0; 
    if(!(i%reset_index))
      sum = 0; 
    sum+=in.at(i); 
    cum_sum.push_back(sum);
  } 
  return cum_sum; 
}

} // namespace Demodulators
