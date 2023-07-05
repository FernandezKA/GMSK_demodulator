/**
 * @file GMSK_demodulator
 * @author FernandezKA
 * @date 2023-07-05
 * @brief Реализация GMSK демодулятора
 * @details
 *  1. Пропускаем I/Q через ФНЧ, таким образом убирая ВЧ;
 *  2. Находим значение фазы для каждого отсчета;
 *  3. Оцениваем, является ли на символьном периоде фаза возрастающей, или
 * убывающей. В зависимости от этого формируем битовый поток до кодирования;
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

GMSK::GMSK(size_t frame_size, size_t baudrate, size_t sample_rate)
    : frame_size(frame_size), baudrate(baudrate), sample_rate(sample_rate) {
  std::cout << "GMSK demodulator is created \n\r";
  input_stream.reserve(frame_size);
}

GMSK::~GMSK() { std::cout << "GMSK demodulator is removed \n\r"; }

int GMSK::Add_Samples(std::vector<std::complex<int16_t>> &input) {
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
  std::vector<int> demodulatedBits;
  size_t sample_per_symbol =
      sample_rate / baudrate; /*Количество отсчетов на один символ*/
  std::vector<double> phase_samples;
  for (const auto &s : input_stream) {
    double tan_phase = 0; 
  
    if(s.real()!= 0){
      tan_phase = static_cast<double>(s.imag()) / static_cast<double>(s.real());
    }
    double phase = atan(tan_phase);
    phase_samples.push_back(phase);
  }
  std::cout << "Phase samples size " << phase_samples.size() << std::endl;

  size_t sample_index = 0;
  double last_phase = 0.0, curr_phase = 0.0;
  while (sample_index < phase_samples.size()) {
    curr_phase = phase_samples.at(sample_index);

    double deltaPhase = curr_phase - last_phase; 
    if(deltaPhase < 0){
      deltaPhase = -deltaPhase; 
    }
    bool bit = (deltaPhase > M_PI * 0.25) ? (1) : (0); 
    output_stream.push_back(bit);

    sample_index += sample_per_symbol;
  }

  std::ofstream atg("atan.txt");
  if (atg.is_open()) {
    for (const auto &s : phase_samples) {
      atg << s << "\n\r";
    }
  }
  return output_stream;
}

/*Простая имплементация ФНЧ, которая необходима перед фильтрацией*/
std::vector<std::complex<int16_t>> &
GMSK::LPF(std::vector<std::complex<int16_t>> &in) {
  /*TODO: Need to  write LPF filter! */
  return in; 
}

bool GMSK::CheckQuadrantSequence(uint8_t last, uint8_t next) const {
  if (next > last && last != 4)
    return true;
  if (last == 4u && next == 0u)
    return true;
  return false;
}

} // namespace Demodulators
