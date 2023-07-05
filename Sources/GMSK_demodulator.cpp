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

/**
 * @brief В данной функции по битовой маске определяем, какие из битов никогда
 * не используются в отсчетах, и двигаем влево на это количество бит.
 * Своеобразная нормировка отсчетов.
 */
size_t GMSK::Normalization(std::vector<std::complex<int16_t>> &samples) {
  int16_t bitmask_shift = 0;

  for (const auto &s : samples) {
    bitmask_shift |= abs(s.real());
    bitmask_shift |= abs(s.imag());
  }

  uint32_t ones = (1u << 31);
  uint8_t shift_val = 0;

  while ((bitmask_shift & ones) == ones) {
    shift_val++;
  }

  for (auto &s : samples) {
    std::complex<int16_t> shifted =
        (s.real() << shift_val, s.imag() << shift_val);
    s = shifted;
  }
  return static_cast<size_t>(shift_val);
}

const std::vector<bool> &GMSK::Get_Bitstream() const { return output_stream; }

/**
 * @brief Демодуляция происходит здесь. Описание приложено в начале файла.
 */
std::vector<bool> GMSK::Demodulate(void) {
  std::vector<int> demodulatedBits;
  double phase = 0.0;
  double curr_phase = 0;
  size_t sample_per_symbol =
      sample_rate / baudrate; /*Количество отсчетов на один символ*/
  std::vector<double> phase_samples;
  for (const auto &s : input_stream) {
    phase_samples.push_back(
        (180.0 / M_PI) *
        atanf(static_cast<float>(s.imag()) / static_cast<float>(s.real())));
  }
  std::cout << "Phase samples size " << phase_samples.size();
  std::ofstream atg("atan.txt");
  if (atg.is_open()) {
    for (const auto &s : phase_samples) {
      atg << s << "\n\r";
    }
  }
  return output_stream;
}

/*Простая имплементация фильтра гаусса, которая необходима перед фильтрацией*/
std::vector<std::complex<int16_t>> &
GMSK::LPF(std::vector<std::complex<int16_t>> &in) {}

bool GMSK::CheckQuadrantSequence(uint8_t last, uint8_t next) const {
  if (next > last && last != 4)
    return true;
  if (last == 4u && next == 0u)
    return true;
  return false;
}

} // namespace Demodulators
