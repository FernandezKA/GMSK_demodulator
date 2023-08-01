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

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fftw3.h>
#include <bitset>
#include <algorithm>
#include <iterator>

using namespace std;

template <typename sample, typename ext_sample>
class GMSK
{
private:
  std::vector<sample> m_iq;
  std::vector<ext_sample> m_ext_iq;
  std::vector<bool> m_bitstream;
  std::vector<ext_sample> preamble_template;

protected:
  size_t m_f_sampling;
  size_t m_baudrate;
  size_t m_sample_per_symbol;

  const size_t MINIMAL_AIS_PACKET_SIZE = 70;
  const size_t PREAMBLE_LEN = 24;

public:
  GMSK(size_t f_sampling, size_t baudrate)
      : m_f_sampling(f_sampling), m_baudrate(baudrate)
  {
    m_sample_per_symbol = m_f_sampling / m_baudrate;
    preamble_template = get_preamble(m_sample_per_symbol);
    assert(preamble_template.size() != 0);
  }
  ~GMSK() {}

  /**
   * @brief загружает I/Q отсчеты для демодуляции
   *
   * @param rhs Raw I/Q
   * @return true
   * @return false
   */
  bool add_samples(const std::vector<sample> &rhs)
  {
    m_iq = rhs;
    m_ext_iq.clear();
    m_ext_iq.resize(m_iq.size());
    for (size_t i = 0; i < m_iq.size(); ++i)
    {
      m_ext_iq.at(i) = get_convert(m_iq.at(i));
    }
    return m_iq.size() == rhs.size();
  }

  inline std::vector<bool> get_bitstream() const { return m_bitstream; }

  /**
   * @brief Используется для принятия жестких решений на основе предоставляемой
   * метрики.
   *
   * @param metrics вектор с метрическими значениями, на основе которых
   * выносятся решения
   * @return std::vector<bool> демодулированный битовый поток
   */
  std::vector<bool> hard_decision(const std::vector<double> &metrics) const
  {
    std::vector<bool> decisions;
    decisions.reserve(metrics.size() / get_sample_per_symbol());
    /*Значения берем из середины символьного интервала*/
    for (size_t i = get_sample_per_symbol() / 2; i < metrics.size();
         i += get_sample_per_symbol())
    {
      decisions.push_back(metrics.at(i) > 0);
    }
    return decisions;
  }

  inline size_t get_baudrate() const { return m_baudrate; }

  inline size_t get_sample_per_symbol() const { return m_sample_per_symbol; }

  inline double get_gain() const
  {
    return (2 * static_cast<double>(get_baudrate())) / M_PI;
  }

  /**
   * @brief Используется для преобразования отсчетов в расширенный формат,
   *  чтобы избежать потерь при умножении из - за ограничений в разрядности.
   *
   * @param rhs исходный отсчет
   * @return ext_sample преобразованный отсчет
   */
  inline ext_sample get_convert(sample rhs) const
  {
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
  inline ext_sample multiply_conjugate(ext_sample iq_0, ext_sample iq_1) const
  {
    return iq_0 * conj(iq_1);
  }

  /**
   * @brief Основное тело демодулятора GMSK/FM binary
   *
   * @return std::vector<bool> демодулированный битовый поток
   */
  std::vector<bool> demodulate()
  {
    std::cout << "Samples size is " << m_iq.size() << std::endl;
    assert(0u != m_iq.size());

    get_log(m_iq, "iq.txt");

    std::vector<ext_sample> diff(m_iq.size());
    std::vector<double> metrics(m_iq.size());

    for (size_t i = 1; i < m_iq.size(); ++i)
    {
      ext_sample shift_sample = get_convert(m_iq.at(i));
      ext_sample curr_sample = get_convert(m_iq.at(i - 1));
      ext_sample mul_conj = multiply_conjugate(shift_sample, curr_sample);
      metrics.at(i - 1) =
          get_gain() * std::atan2(mul_conj.imag(), mul_conj.real());
      diff.at(i - 1) = mul_conj;
    }

    get_log(diff, "diff.txt");
    get_log(metrics, "metrics.txt");

    m_bitstream = hard_decision(metrics);

    return m_bitstream;
  }

  /**
   * @brief Определение среднего значения на заданном интервале
   *
   * @return {state, average}
   */
  inline sample get_average(const std::vector<sample> &sample_, size_t initial_pos_, size_t length_) const
  {
    int sum_re = 0;
    int sum_im = 0;
    for (size_t i = initial_pos_; i < initial_pos_ + length_; ++i)
    {
      sum_re += sample_.at(i).real();
      sum_im += sample_.at(i).imag();
    }
    sum_re /= length_;
    sum_im /= length_;
    return static_cast<sample>(sum_re, sum_im);
  }

  /**
   * @brief В данной функции мы демодулируем сообщения, детектируя начало пакета
   */
  bool demodulate_without_noise()
  {
    size_t avg_packet_position = detect_packet();

    auto iq_begin = m_ext_iq.begin();

    auto corr = get_correlation(iq_begin + avg_packet_position, preamble_template.begin(), preamble_template.size());

    // std::cout << "Maximal correlation with template preamble on position " << max(corr) << std::endl;

    return true;
  }

  /**
   * @brief Первичное определение начала пакета по нарастанию амплитуды. Если значение амплитуды для конкретного участка выше среднего значения
   * по реализации, то мы считаем, что в этом месте гипотетически может быть пакет. После чего проверяем, соответствует ли длительность импульса
   * размеру пакета АИС, или это был шум. После этого, если два критерия выполнены, мы пытаемся про помощи коррелляции с преамбулой детектировать
   * точное начало пакета, с которого будем демодулировать сообщение.
   * @return позиция начала пакета (включая преамбулу) в реализации, которую мы обрабатываем.
   */
  size_t detect_packet() const
  {
    size_t packet_begin_index = 0;
    sample average = get_average(m_iq, 0, m_iq.size());
    for (size_t i = 0; i < m_iq.size(); ++i)
    {
      if (abs(m_iq.at(i)) > abs(average))
      {
        auto local_average = get_average(m_iq, i, m_sample_per_symbol * MINIMAL_AIS_PACKET_SIZE);
        if (abs(local_average) / abs(average) >= 2)
        { /*Пока что просто считаем начало пакета по уровню 3 dB, без порога*/
          std::cout << "Packet begin on position : " << i << std::endl;
          packet_begin_index = i;
        }
      }
    }
    return packet_begin_index;
  }

  /**
   * @brief Здесь мы формируем массив с преамбулой, для которого мы будем применять функцию корелляции, чтобы детектировать начало пакета.
   *
   */
  inline std::vector<ext_sample> get_preamble(size_t samples_per_symbol) const
  {
    /*Период повторения синусоиды в преамбуле составляет 4 * кол - вот отсчетов на символьный интервал
     * Количество отсчетов на символьный интервал мы можем посчитать, как Fs * BAUDRATE, для нашей прошивки это ~10.4 отсчетов на символьный интервал
     */
    const double preamble_length = static_cast<double>(PREAMBLE_LEN);
    std::vector<ext_sample> template_preamble(static_cast<size_t>(preamble_length) * samples_per_symbol);
    for (size_t i = 0; i < template_preamble.size(); ++i)
    {
      double real = cos((2 * M_PI * preamble_length * static_cast<double>(i) / static_cast<double>(template_preamble.size())) / pow(
                                                                                                                                    static_cast<double>(samples_per_symbol), 2));
      double imag = sin((2 * M_PI * preamble_length * static_cast<double>(i) / static_cast<double>(template_preamble.size())) / pow(
                                                                                                                                    static_cast<double>(samples_per_symbol), 2));
      ext_sample preamble_sample = {real, imag};
      template_preamble.at(i) = preamble_sample;
    }
    return template_preamble;
  }

  /**
   * @brief Кодирует исходный битовый поток с использованием NZRI
   *
   * @return битовый поток после кодирования
   */
  inline std::vector<bool> NRZI_encode(const std::vector<bool> &rhs) const
  {
    std::vector<bool> coded_;
    coded_.reserve(rhs.size());
    bool curr_bit = true;
    for (size_t i = 1; i < rhs.size(); ++i)
    {
      if (rhs.at(i) == 1)
      {
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
  inline std::vector<bool> NRZI_decode(const std::vector<bool> &encoded) const
  {
    std::vector<bool> decoded_;
    decoded_.reserve(encoded.size());
    bool currentBit = true; // start with high bit
    for (int i = 0; i < encoded.size(); i++)
    {
      if (encoded[i])
      {
        // if signal changes, add 1 to decoded string
        decoded_.push_back(currentBit);
        currentBit = !currentBit;
      }
      else
      {
        // if signal stays the same, add 0 to decoded string
        decoded_.push_back(currentBit);
      }
    }
    return decoded_;
  }

  void fft(std::vector<std::complex<double>> &x, size_t n, int sign)
  {
    // Проверка на то, что n является степенью двойки
    if ((n & (n - 1)) != 0)
    {
      throw std::invalid_argument("n must be a power of 2");
    }

    // Выполнение FFT
    fft_rec(x, n, sign);
  }

  void fft_rec(std::vector<std::complex<double>> &x, size_t n, int sign)
  {
    if (n == 1)
    {
      return;
    }

    // Разделение входного вектора на четные и нечетные элементы
    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);

    for (size_t i = 0; i < n / 2; i++)
    {
      even[i] = x[2 * i];
      odd[i] = x[2 * i + 1];
    }

    // Рекурсивное выполнение FFT для четных и нечетных подвекторов
    fft_rec(even, n / 2, sign);
    fft_rec(odd, n / 2, sign);

    // Вычисление значения W_n^k
    std::complex<double> W_n_k = std::polar(1.0, sign * 2 * M_PI / n);

    // Выполнение бабочечной операции
    std::complex<double> w = 1.0;
    for (size_t i = 0; i < n / 2; i++)
    {
      std::complex<double> t = w * odd[i];
      x[i] = even[i] + t;
      x[i + n / 2] = even[i] - t;
      w *= W_n_k;
    }
  }

  void ifft(std::vector<std::complex<double>> &x, int n)
  {
    fft(x, n, -1);

    // Нормализация результата обратного FFT
    for (size_t i = 0; i < n; i++)
    {
      x[i] /= n;
    }
  }

  /**
   * @brief Реализация получения корелляции с использованием БПФ.
   *Корреляция двух сигналов может быть вычислена при помощи БПФ для комплексных чисел следующим образом:
    1. Дополнить каждый сигнал нулями до длины, равной сумме их длин минус один.
    2. Применить БПФ к каждому из дополненных сигналов.
    3. Умножить полученные спектры двух сигналов поэлементно.
    4. Применить обратное БПФ к полученному произведению спектров.
    5. Взять действительную часть результата обратного БПФ как значение корреляции для каждого сдвига.

    Этот алгоритм основан на свойстве корреляции, что она равна обратному преобразованию Фурье произведения преобразований Фурье двух сигналов.
    Для комплексных чисел необходимо использовать сопряжение вместо транспонирования при умножении спектров.

    Применение БПФ и IFFT позволяет вычислить корреляцию за O(n log n) времени вместо O(n^2) при использовании прямого вычисления корреляции.
   * @return значения корелляции для каждого сдвига
   */
  std::vector<double> get_correlation(std::vector<std::complex<double>>::const_iterator x, std::vector<std::complex<double>>::const_iterator y, size_t n)
  {
    std::vector<double> result(n);

    // Дополнение векторов нулями
    std::vector<std::complex<double>> x_padded(2 * n);
    std::vector<std::complex<double>> y_padded(2 * n);
    std::copy(x, x + n, x_padded.begin());
    std::copy(y, y + n, y_padded.begin());

    // Выполнение FFT для дополненных векторов
    fft(x_padded, 2 * n, 1);
    fft(y_padded, 2 * n, 1);

    // Вычисление произведения спектров двух векторов поэлементно
    for (size_t i = 0; i < 2 * n; i++)
    {
      x_padded[i] *= std::conj(y_padded[i]);
    }

    // Выполнение обратного FFT для произведения спектров
    ifft(x_padded, 2 * n);

    // Вывод действительной части результата обратного FFT как значение корреляции для каждого сдвига
    for (size_t i = 0; i < n; i++)
    {
      result[i] = x_padded[i].real();
    }

    return result;
  }

  /**
   * @brief Логгирование данных в файл
   * @param rhs логгируемые данные
   * @param path путь к файлу логирования
   * @return true
   * @return false
   */
  bool get_log(const std::vector<sample> &rhs, std::string path) const
  {
    std::ofstream dump(path);
    if (dump.is_open())
    {
      for (const auto &s : rhs)
      {
        dump << s.real() << " " << s.imag() << std::endl;
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  bool get_log(const std::vector<ext_sample> &rhs, std::string path) const
  {
    std::ofstream dump(path);
    if (dump.is_open())
    {
      for (const auto &s : rhs)
      {
        dump << s.real() << " " << s.imag() << std::endl;
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  bool get_log(const std::vector<double> &rhs, std::string path) const
  {
    std::ofstream dump(path);
    if (dump.is_open())
    {
      for (const auto &s : rhs)
      {
        dump << s << std::endl;
      }
      return true;
    }
    else
    {
      return false;
    }
  }
};