/**
 * @file GMSK_demodulator
 * @author FernandezKA
 * @date 2023-07-05
 * @brief Реализация GMSK демодулятора
 * @details
 *  1. Пропускаем I/Q через ФНЧ, таким образом убирая ВЧ;
 *  2. Находим значение фазы для каждого отсчета;
 *  3. Оцениваем разность фаз в точках семплирования, по этому значению смотрим,
 * был ли скачок фазы.
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
#include <vector>

namespace Demodulators {

GMSK::GMSK() {
  std::cout << "GMSK demodulator is created \n\r";
  input_stream.reserve(frame_size);
}

GMSK::~GMSK() { std::cout << "GMSK demodulator is removed \n\r"; }

/**
 * @brief Добавление битового потока для обработки. Все дальнейшие операции производятся над ним. 
 * 
 * @param input обрабатываемый битовый поток
 * @return int 
 */
int GMSK::Add_Samples(std::vector<std::complex<double>> &input) {
  input_stream.clear();
  for (const auto &s : input)
    input_stream.push_back(s);
  return (input_stream.size() == input.size()) ? (0) : (-1);
}

const std::vector<bool> &GMSK::Get_Bitstream() const { return output_stream; }

/**
 * @brief Имплементация Binary GMSK/FM демодулятора. 
 * 1. Детектировали фазу у каждого отсчета
 * 2. Проинтегрировали значение в каждом символьном интервале 
 * 3. Сравнили крайнее значение куммулятивной суммы с прошлым временым интервалом. По результату, если 
 * значение фазы увеличилось - 1 бит, иначе 0 бит был передан. 
 * @return std::vector<bool> битовый поток после демодуляции
 */
std::vector<bool> GMSK::Demodulate(void) {
  std::vector<double> phase_samples;
  for (auto &s : input_stream) {
    phase_samples.push_back(abs(std::atan2(s.imag(), s.real())));
  }

  auto cum_sum = CumSum(phase_samples, samples_per_symbol);
  std::ofstream phase("phase.txt");
  for (auto s : cum_sum) {
    phase << s << std::endl;
  }

  for (size_t i = samples_per_symbol - 1; i < phase_samples.size();
       i += samples_per_symbol) {
    static double last_phase = 0;
    bool bit = (cum_sum.at(i) - last_phase > 0);
    last_phase = cum_sum.at(i);
    output_stream.push_back(bit);
  }
  std::cout << "Phase samples size " << phase_samples.size() << std::endl;

  return output_stream;
}
/**
 * @brief Реализация ФНЧ, который используется перед демодуляцией сигнала. 
 * 
 * @param in I/Q отсчеты
 * @return std::vector<std::complex<int16_t>>& Обработанный поток I/Q отсчетов. 
 */
std::vector<std::complex<int16_t>> &
GMSK::LPF(std::vector<std::complex<int16_t>> &in) {
  /*TODO: Need to  write LPF filter! */
  return in;
}

/**
 * @brief Используется вместо интегратора, для усреднения значений на одном символьном интервале. 
 * 
 * @param in входные данные
 * @param reset_index размер окна суммирования
 * @return std::vector<double> 
 */
std::vector<double> GMSK::CumSum(std::vector<double> &in,
                                 size_t reset_index) const {
  std::vector<double> cum_sum;
  cum_sum.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    static double sum = 0;
    if (!(i % reset_index))
      sum = 0;
    sum += in.at(i);
    cum_sum.push_back(sum);
  }
  return cum_sum;
}

/**
 * @brief В АИС используется вставка нулевого бита, если мы находим блок последовательных единичек, размером более 5. 
 * В этой функции мы убираем эти вставки, таким образом восстанавливая исходный битовый поток. 
 * TODO: Это не физический уровень! Семантически оно не должно быть тут. 
 * @param in битовый поток, после демодуляции
 * @return std::vector<bool> 
 */
std::vector<bool> GMSK::RemoveInsertions(const std::vector<bool> &in) const {
  std::vector<bool> parsed;
  parsed.reserve(in.size());
  size_t counter = 0;
  for (const auto &s : in) {
    if (s) {
      ++counter;
    }
    if (counter == 5 && !s) {
      counter = 0;
      continue;
    } else {
      parsed.push_back(s);
    }
  }
  return parsed;
}

} // namespace Demodulators
