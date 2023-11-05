/**
 * @file GMSK_demodulator.hpp
 * @author FernandezKA (fernandes.kir@yandex.ru)
 * @brief GMSK demodulator
 * @version 0.1
 * @date 2023-07-22
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>
#include <iterator>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

template <typename sample, typename ext_sample> class GMSK {
private:
  std::vector<sample> m_iq;
  std::vector<ext_sample> m_ext_iq;
  std::vector<bool> m_bitstream;
  std::vector<ext_sample> preamble_template;
  size_t preamble_zeros_length = 0;

protected:
  size_t m_f_sampling;
  size_t m_baudrate;
  size_t m_sample_per_symbol;

  const size_t MINIMAL_AIS_PACKET_SIZE = 70;
  const size_t PREAMBLE_LEN = 24;

  const std::vector<bool> START_STOP_FLAG = {0, 1, 1, 1, 1, 1, 1, 0};

public:
  /**
   * @param f_samling frequency of sampling [Hz]
   * @param BAUDRATE baudrate for current demodulator. For example, AIS used
   * BAUDRATE = 9600U
   */
  GMSK(size_t f_sampling, size_t baudrate)
      : m_f_sampling(f_sampling), m_baudrate(baudrate) {
    m_sample_per_symbol = m_f_sampling / m_baudrate;
    auto preamble_info = get_preamble(m_sample_per_symbol);
    preamble_template = preamble_info.first;
    preamble_zeros_length = preamble_info.second;
    assert(preamble_template.size() != 0);
  }
  ~GMSK() {}

  /**
   * Block of user functions
   */
  /**
   * @brief загружает I/Q отсчеты для демодуляции
   *
   * @param rhs Raw I/Q
   * @return true
   * @return false
   */
  bool add_samples(const std::vector<sample> &rhs) {
    m_iq = rhs;
    m_ext_iq.clear();
    m_ext_iq.resize(m_iq.size());
    for (size_t i = 0; i < m_iq.size(); ++i) {
      m_ext_iq.at(i) = get_convert(m_iq.at(i));
    }
    get_log(m_iq, "iq.txt");
    return 0u != m_iq.size();
  }

  inline std::vector<bool> get_bitstream() const { return m_bitstream; }

  inline size_t get_baudrate() const { return m_baudrate; }

  inline size_t get_sample_per_symbol() const { return m_sample_per_symbol; }

  inline double get_gain() const {
    return (2 * static_cast<double>(get_baudrate())) / M_PI;
  }

  /**
   * Block of internal conversion functions
   */

  /**
   * @brief Используется для преобразования отсчетов в расширенный формат,
   *  чтобы избежать потерь при умножении из - за ограничений в разрядности.
   *
   * @param rhs исходный отсчет
   * @return ext_sample преобразованный отсчет
   */
  inline ext_sample get_convert(sample rhs) const {
    double real = static_cast<double>(rhs.real());
    double imag = static_cast<double>(rhs.imag());
    ext_sample res(real, imag);
    return res;
  }

  /**
   * @brief умножение двух чисел между собой, второе из которых является
   * комплексно сопряженным.
   *
   * @param iq_0
   * @param iq_1
   * @return ext_sample
   */
  inline ext_sample multiply_conjugate(ext_sample iq_0, ext_sample iq_1) const {
    return iq_0 * conj(iq_1);
  }

  /**
   * @brief Определение среднего значения на заданном интервале
   *
   * @return {state, average}
   */
  inline sample get_average(const std::vector<sample> &sample_,
                            size_t initial_pos_, size_t length_) const {
    int sum_re = 0;
    int sum_im = 0;
    for (size_t i = initial_pos_; i < initial_pos_ + length_; ++i) {
      sum_re += sample_.at(i).real();
      sum_im += sample_.at(i).imag();
    }
    sum_re /= static_cast<int>(length_);
    sum_im /= static_cast<int>(length_);
    sample result(sum_re, sum_im);
    return result;
  }

  /**
   * Block of packet detector part of code
   */

  /**
   * @brief В данной функции мы демодулируем сообщения, детектируя начало пакета
   */
  std::pair<size_t, size_t> try_to_find_packet() {
    auto packet_index = detect_packet();
    std::cout << "Packet detected from " << packet_index.first << " to "
              << packet_index.second << std::endl;
    double norm_coeff = abs(get_average(
        m_iq, packet_index.first, packet_index.second - packet_index.first));
    size_t packet_precision_shift = static_cast<size_t>(
        get_shift_of_packet(m_ext_iq, packet_index.first, norm_coeff));
    std::cout << "Maximal correlation with template preamble on position (from "
                 "amplitude detected position) "
              << packet_precision_shift << std::endl;
    return {packet_index.first + packet_precision_shift,
            packet_index.second + packet_precision_shift};
  }

  /**
   * @brief Первичное определение начала пакета по нарастанию амплитуды. Если
   * значение амплитуды для конкретного участка выше среднего значения по
   * реализации, то мы считаем, что в этом месте гипотетически может быть пакет.
   * @return позиция начала пакета (включая преамбулу) в реализации, которую мы
   * обрабатываем.
   */
  std::pair<size_t, size_t> detect_packet() const {
    sample average = get_average(m_iq, 0, m_iq.size());
    std::cout << "Average of samples: " << average << std::endl;
    size_t packet_begin = 0;
    size_t packet_end = 0;
    for (size_t i = 0; i < m_iq.size(); ++i) {
      if (abs(m_iq.at(i)) > 2 * abs(average)) /**3 dB threshold value*/
      {
        if (!packet_begin)
          packet_begin = i;
        packet_end = i;
      }
    }
    return {packet_begin, packet_end};
  }

  /**
   * @brief Уточнение позиции начала пакета при помощи корелляции с
   * детерминированной нормированной преамбуле.
   * @return Сдвиг позиции положения пакета в правую сторону, необходимый для
   * достижения максимальной корелляции преамбул.
   */
  int get_shift_of_packet(std::vector<ext_sample> &realization,
                          size_t packet_position, double normalization_coeff) {
    auto packet_iterator = realization.begin() + packet_position;
    std::vector<ext_sample> preamble_normed(preamble_template.size());
    for (size_t i = 0; i < preamble_template.size(); ++i) {
      preamble_normed.at(i) = preamble_template.at(i) * normalization_coeff;
    }
    auto preamble_iterator =
        preamble_normed
            .begin(); /*Домножаем преамбулу на размерность принятой реализации*/
    auto correlation = get_correlation(preamble_iterator, packet_iterator,
                                       preamble_normed.size());
    get_log(correlation, "correlation.txt");
    double max_corr = 0;
    size_t max_corr_pos = 0;
    for (size_t i = 0; i < correlation.size(); ++i) {
      if (correlation.at(i) > max_corr) {
        max_corr = correlation.at(i);
        max_corr_pos = i;
      }
    }
    return max_corr_pos;
  }

  /**
   * Block of functions for correlation finding
   */

  void fft(std::vector<std::complex<double>> &x, size_t n, int sign) {
    // Проверка на то, что n является степенью двойки
    if ((n & (n - 1)) != 0) {
      throw std::invalid_argument("n must be a power of 2");
    }

    // Выполнение FFT
    fft_rec(x, n, sign);
  }

  void fft_rec(std::vector<std::complex<double>> &x, int n, int sign) {
    if (n == 1) {
      return;
    }

    // Разделение входного вектора на четные и нечетные элементы
    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);

    for (size_t i = 0; i < n / 2; i++) {
      even[i] = x[2 * i];
      odd[i] = x[2 * i + 1];
    }

    // Рекурсивное выполнение FFT для четных и нечетных подвекторов
    fft_rec(even, n / 2, sign);
    fft_rec(odd, n / 2, sign);

    // Вычисление значения W_n^k
    std::complex<double> W_n_k =
        std::polar(1.0, sign * 2 * M_PI / static_cast<double>(n));

    // Выполнение бабочечной операции
    std::complex<double> w = 1.0;
    for (size_t i = 0; i < n / 2; i++) {
      std::complex<double> t = w * odd[i];
      x[i] = even[i] + t;
      x[i + n / 2] = even[i] - t;
      w *= W_n_k;
    }
  }

  void ifft(std::vector<std::complex<double>> &x, int n) {
    fft(x, n, -1);

    // Нормализация результата обратного FFT
    for (size_t i = 0; i < n; i++) {
      x[i] /= n;
    }
  }

  /**
   * @brief Реализация получения корелляции с использованием БПФ.
   *Корреляция двух сигналов может быть вычислена при помощи БПФ для комплексных
   чисел следующим образом:
    1. Дополнить каждый сигнал нулями до длины, равной сумме их длин минус один.
    2. Применить БПФ к каждому из дополненных сигналов.
    3. Умножить полученные спектры двух сигналов поэлементно.
    4. Применить обратное БПФ к полученному произведению спектров.
    5. Взять модуль результата обратного БПФ как значение корреляции для каждого
   сдвига.

    Этот алгоритм основан на свойстве корреляции, что она равна обратному
   преобразованию Фурье произведения преобразований Фурье двух сигналов. Для
   комплексных чисел необходимо использовать сопряжение вместо транспонирования
   при умножении спектров.

    Применение БПФ и ОБПФ позволяет вычислить корреляцию за O(n log n) времени
   вместо O(n^2) при использовании прямого вычисления корреляции.
   * @return значения корелляции для каждого сдвига
   */
  std::vector<double>
  get_correlation(std::vector<std::complex<double>>::const_iterator x,
                  std::vector<std::complex<double>>::const_iterator y, int n) {
    std::vector<double> result(static_cast<std::vector<double>::size_type>(n));

    // Дополнение векторов нулями
    std::vector<std::complex<double>> x_padded(2 * n);
    std::vector<std::complex<double>> y_padded(2 * n);
    std::copy(x, x + n, x_padded.begin());
    std::copy(y, y + n, y_padded.begin());

    std::cout << "Size of preamble is " << n << std::endl;
    std::cout << "Size of vectors before FFT " << x_padded.size() << " "
              << y_padded.size() << std::endl;

    get_log(x_padded, "x_padded.txt");
    get_log(y_padded, "y_padded.txt");

    // Выполнение FFT для дополненных векторов
    fft(x_padded, 2 * n, 1);
    fft(y_padded, 2 * n, 1);

    get_log(x_padded, "x_spectrum.txt");
    get_log(y_padded, "y_spectrum.txt");

    // Вычисление произведения спектров двух векторов поэлементно
    for (size_t i = 0; i < 2 * n; i++) {
      x_padded[i] *= std::conj(y_padded[i]);
    }

    get_log(x_padded, "spectrum_multiply.txt");

    // Выполнение обратного FFT для произведения спектров
    ifft(x_padded, 2 * n);

    // Вывод действительной части результата обратного FFT как значение
    // корреляции для каждого сдвига
    for (size_t i = 0; i < n; i++) {
      result[i] = abs(x_padded[i]);
    }

    get_log(result, "corr.txt");

    return result;
  }

  /**
   * @brief Здесь мы формируем массив с преамбулой, для которого мы будем
   * применять функцию корелляции, чтобы детектировать начало пакета.
   */
  inline std::pair<std::vector<ext_sample>, size_t>
  get_preamble(size_t samples_per_symbol) const {
    /*Период повторения синусоиды в преамбуле составляет 4 * кол - вот отсчетов
     * на символьный интервал Количество отсчетов на символьный интервал мы
     * можем посчитать, как Fs * BAUDRATE, для нашей прошивки это ~10.4 отсчетов
     * на символьный интервал Последнюю позицию дополняем нулями до размера
     * ближайшей степени двойки, чтобы использовать БПФ в дальшейшем.
     */
    const double preamble_length = static_cast<double>(PREAMBLE_LEN);
    std::vector<ext_sample> template_preamble(
        static_cast<size_t>(preamble_length) * samples_per_symbol);
    for (size_t i = 0; i < template_preamble.size(); ++i) {
      double t = static_cast<double>(i) / (samples_per_symbol * PREAMBLE_LEN);
      double phase = 2 * M_PI * PREAMBLE_LEN / 4 * t + M_PI_2;
      double real = cos(phase);
      double imag = sin(phase);
      ;
      ext_sample preamble_sample = {real, imag};
      template_preamble.at(i) = preamble_sample;
    }
    size_t minimal_power = 0;
    while (template_preamble.size() > pow(2, minimal_power))
      minimal_power++;
    size_t count_zeros_added = pow(2, minimal_power) - template_preamble.size();
    while (template_preamble.size() < pow(2, minimal_power))
      template_preamble.push_back(0);
    std::cout << "Template of preamble size is " << template_preamble.size()
              << std::endl;
    get_log(template_preamble, "template_preamble.txt");
    return {template_preamble, count_zeros_added};
  }

  /**
   * Main logic of GMSK demodulator
   */

  /**
   * @brief Используется для принятия жестких решений на основе предоставляемой
   * метрики.
   *
   * @param metrics вектор с метрическими значениями, на основе которых
   * выносятся решения
   * @return std::vector<bool> демодулированный битовый поток
   */
  std::vector<bool> hard_decision(const std::vector<double> &metrics) const {
    std::vector<bool> decisions;
    decisions.reserve(metrics.size() / get_sample_per_symbol());
    for (size_t i = get_sample_per_symbol();
         i < metrics.size() - get_sample_per_symbol();
         i += get_sample_per_symbol()) {
      double sum = 0;
      for (size_t j = i; j < i + get_sample_per_symbol(); ++j) {
        sum += metrics.at(j);
      }
      decisions.push_back(sum < 0);
    }
    return decisions;
  }
  /**
   * @brief Основное тело демодулятора GMSK/FM binary
   *
   * @return std::vector<bool> демодулированный битовый поток
   */
  std::vector<bool> demodulate(std::pair<size_t, size_t> packet_detected) {
    size_t realization_size =
        packet_detected.second - packet_detected.first + 1;
    std::cout << "Samples size is " << realization_size << std::endl;
    if (0u == realization_size)
      return {0};
    get_log(m_iq, "iq.txt");

    std::vector<ext_sample> diff(realization_size);
    std::vector<double> metrics(realization_size);

    if (packet_detected.second > m_iq.size())
      packet_detected.second = m_iq.size();
    for (size_t i = packet_detected.first + 1; i < packet_detected.second - 1;
         ++i) {
      ext_sample shift_sample = get_convert(m_iq.at(i));
      ext_sample curr_sample = get_convert(m_iq.at(i - 1));
      ext_sample mul_conj = multiply_conjugate(shift_sample, curr_sample);
      metrics.at(i - 1 - packet_detected.first) =
          get_gain() * std::atan2(mul_conj.imag(), mul_conj.real());
      diff.at(i - 1 - packet_detected.first) = mul_conj;
    }

    get_log(diff, "diff.txt");
    get_log(metrics, "metrics.txt");

    m_bitstream = hard_decision(metrics);

    return m_bitstream;
  }

  /**
   * Block of code for working with channel layout
   */

  /**
   * @brief Кодирует исходный битовый поток с использованием NZRI
   *
   * @return битовый поток после кодирования
   */
  inline std::vector<bool> NRZI_encode(const std::vector<bool> &rhs) const {
    std::vector<bool> coded_;
    coded_.reserve(rhs.size());
    bool curr_bit = true;
    for (size_t i = 1; i < rhs.size(); ++i) {
      if (rhs.at(i) == 1) {
        curr_bit = !curr_bit;
      }
      coded_.push_back(curr_bit);
    }
    return coded_;
  }

  /**
   * @brief Декодирует битовый поток при помощи NRZI в исходный битовый поток
   *
   * @return Исходный битовый поток до декодирования
   */
  inline std::vector<bool> NRZI_decode(const std::vector<bool> &encoded) const {
    std::vector<bool> decoded_;
    decoded_.reserve(encoded.size());
    bool currentBit = true; // start with high bit
    for (int i = 0; i < encoded.size(); i++) {
      if (encoded[i]) {
        // if signal changes, add 1 to decoded string
        decoded_.push_back(currentBit);
        currentBit = !currentBit;
      } else {
        // if signal stays the same, add 0 to decoded string
        decoded_.push_back(currentBit);
      }
    }
    return decoded_;
  }

  /**
   * @brief Декодирует битовый поток при помощи NRZI в исходный битовый поток
   *
   * @return Исходный битовый поток до декодирования
   */
  inline std::vector<bool>
  NRZI_decode(std::vector<bool>::const_iterator encoded, size_t len) const {
    bool last_bit = encoded[0];
    std::vector<bool> decoded_;
    decoded_.reserve(len);
    for (size_t i = 1; i < len; ++i) {
      if (last_bit ==
          encoded[static_cast<std::iterator<std::random_access_iterator_tag,
                                            bool>::difference_type>(i)]) {
        decoded_.push_back(1);
      } else {
        decoded_.push_back(0);
      }
      last_bit = encoded[i];
    }
    return decoded_;
  }

  /**
   * @brief Remove insertions each of five bits into bitstream
   * @return std::vector<bool> with removed bits
   */

  std::vector<bool>
  get_remove_insertions(const std::vector<bool> &stream) const {
    std::vector<bool> bitstream;
    bitstream.reserve(stream.size());
    size_t counter = 0;
    bool last_bit = stream.at(0);
    bitstream.push_back(stream.at(0));
    for (size_t i = 1; i < stream.size(); ++i) {
      if (stream.at(i) == last_bit) {
        ++counter;
      } else {
        counter = 0;
      }
      if (counter < 5)
        bitstream.push_back(stream.at(i));
    }
    return bitstream;
  }

  /**
   * @brief Check preamble into demodulated bitstream (AFTER NRZI decoding!!!)
   * @return true preamble has been detected
   * @return false preamble hasn't been detected
   */
  inline bool is_preamble_detect(const std::vector<bool> &rhs) const {
    size_t errors = 0;
    if (rhs.size() < PREAMBLE_LEN - 1)
      return false;
    for (size_t i = 1; i < rhs.size(); ++i) {
      if (rhs.at(i - 1) == rhs.at(i))
        errors++;
      if (2u < errors)
        return false;
    }
    return true;
  }

  /**
   * @brief Chech Start/Stop flags
   * @return true start/stop flag exist
   * @return false mistake into flag detecting
   */

  inline bool is_start_stop_flag(std::vector<bool>::const_iterator rhs) const {
    size_t error = 0u;
    bool status = true;
    for (size_t i = 0; i < START_STOP_FLAG.size(); ++i) {
      status &= rhs[i] == START_STOP_FLAG.at(i);
    }
    return status;
  }

  /**
   * Block of function for debug
   */

  /**
   * @brief Логгирование данных в файл
   * @param rhs логгируемые данные
   * @param path путь к файлу логирования
   * @return true
   * @return false
   */
  template <typename T>
  bool get_log(const std::vector<sample> &rhs, std::string path) const {
    std::ofstream dump(path, std::ios::app);
    if (dump.is_open()) {
      dump.write(reinterpret_cast<char *>(&rhs[0]), sizeof(T) * rhs.size());
    } else {
      return false;
    } 
  }
};

/*End of GMSK_demodulator.hpp file*/
