Руководство пользователя Bell Optimizer
1. Введение
Bell Optimizer — это программа для автоматического проектирования колоколов с заданными акустическими характеристиками. Программа самостоятельно перебирает различные варианты геометрии, запускает расчёты в решателе CalculiX и выбирает оптимальную форму.
Программа состоит из двух частей:
• 
bell_opt.exe — оптимизатор, ищет оптимальную геометрию колокола
• 
bell_synth.exe — синтезатор, создаёт звуковой файл из результатов расчёта
 
2. Системные требования
• 
Операционная система: Windows (UCRT64 окружение)
• 
Установленный решатель CalculiX (для оптимизатора)
• 
Компилятор GCC и система сборки CMake (для сборки из исходников)
 
3. Установка и сборка
3.1. Подготовка окружения
Откройте терминал UCRT64 (MSYS2) и установите необходимые библиотеки:
bash
pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake
pacman -S mingw-w64-ucrt-x86_64-gsl mingw-w64-ucrt-x86_64-nlopt
pacman -S mingw-w64-ucrt-x86_64-fftw mingw-w64-ucrt-x86_64-eigen3
pacman -S mingw-w64-ucrt-x86_64-calculix
3.2. Сборка проекта
bash
cd /path/to/bell_optimizer
mkdir build
cd build
cmake ..
make -j4
После сборки в папке build/ появятся исполняемые файлы:
• 
bell_opt.exe — оптимизатор
• 
bell_synth.exe — синтезатор
 
4. Структура проекта
text
bell_optimizer/
├── build/                      # Сборка
│   ├── bell_opt.exe            # Оптимизатор
│   └── bell_synth.exe          # Синтезатор
├── bells/                      # Рабочая папка для расчётов (рекомендуется)
│   ├── bell.inp                # Входной файл для CalculiX
│   ├── bell.dat                # Результаты расчёта
│   ├── best_geometry.json      # Лучшая найденная геометрия
│   ├── best_bell.inp           # Копия bell.inp
│   └── best_bell.dat           # Копия bell.dat
├── targets.json                # Целевые параметры (частота, соотношения)
├── evolution.json              # Параметры эволюционного алгоритма
├── synthesis.json              # Параметры синтезатора
└── optimize.cmd                # Командный файл для запуска
 
5. Файлы конфигурации
5.1. targets.json — Целевые параметры
Определяет, какой звук должен получиться у колокола.
json
{
    "prime_freq_target": 432.0,
    "freq_ratios": [0.5, 1.25, 1.5, 2.0],
    "amp_ratios": [0.8, 0.6, 0.4, 0.2],
    "freq_weights": [2.0, 1.0, 1.0, 1.0],
    "amp_weights": [1.0, 1.0, 1.0, 1.0],
    "prime_weight": 0.0,
    "k_prime": 20.0,
    "k_freq": 20.0,
    "k_amp": 5.0
}
Параметры:
Поле
Описание
prime_freq_target 
Целевая частота Prime (основного тона) в герцах
freq_ratios 
Соотношения частот: [Hum/Prime, Tertia/Prime, Quinta/Prime, Nominal/Prime]
amp_ratios 
Соотношения амплитуд: [Prime/Hum, Tertia/Hum, Quinta/Hum, Nominal/Hum]
freq_weights 
Веса ошибок для частот: [Hum, Tertia, Quinta, Nominal]
amp_weights 
Веса ошибок для амплитуд: [Prime, Tertia, Quinta, Nominal]
prime_weight 
Вес ошибки Prime (обычно 0, т.к. масштабируется)
k_prime
Крутизна сигмоиды для Prime
k_freq
Крутизна сигмоиды для частот
k_amp 
Крутизна сигмоиды для амплитуд
Пример для православного колокола 432 Гц:
json
{
    "prime_freq_target": 432.0,
    "freq_ratios": [0.5, 1.25, 1.5, 2.0],
    "amp_ratios": [0.8, 0.6, 0.4, 0.2],
    "freq_weights": [2.0, 1.0, 1.0, 1.0],
    "amp_weights": [1.0, 1.0, 1.0, 1.0],
    "prime_weight": 0.0,
    "k_prime": 20.0,
    "k_freq": 20.0,
    "k_amp": 5.0
}
 
5.2. evolution.json — Параметры оптимизации
Управляет процессом поиска оптимальной формы.
json
{
    "population_size": 30,
    "generations": 30,
    "num_elite": 6,
    "mutation_probability": 0.1,
    "mutation_strength": 0.05,
    "num_knots": 5
}
Параметры:
Поле
Описание
population_size
Количество кандидатов в популяции (больше = шире поиск)
generations 
Количество поколений эволюции (больше = точнее)
num_elite 
Количество лучших, переходящих в следующее поколение
mutation_probability
Вероятность мутации (0.0-1.0)
mutation_strength
Сила мутации
num_knots 
Количество узлов сплайна (включая макушку и устье)
 
5.3. synthesis.json — Параметры синтезатора
Управляет синтезом звука из результатов расчёта.
json
{
    "sample_rate": 96000,
    "oversample_factor": 8,
    "duration": 0.0,
    "fade_in": 0.01,
    "fade_out": 0.2,
    "loss_factor": 0.0002,
    "nonlinear_damping": 0.5,
    "tail_factor": 3.0,
    "envelope_threshold": 0.001,
    "boost_amplitude": 0.0,
    "boost_freq_factor": 0.5,
    "slowdown_factor": 0.15,
    "impact_position": 0.70,
    "impact_duration": 0.0005,
    "impact_force": 1.0,
    "max_frequency": 20000.0,
    "use_acoustic_model": false,
    "use_acoustic_coupling": false,
    "use_helmholtz_resonance": false,
    "helmholtz_amplitude": 0.04,
    "helmholtz_damping": 0.3,
    "radiation_efficiency": 0.4,
    "radiation_damping_factor": 0.005,
    "use_mount_model": false,
    "mount_damping": 0.2,
    "reverb_mix": 0.0,
    "reverb_delay": 0.1,
    "reverb_feedback": 0.5,
    "use_rms_normalization": false,
    "target_rms": 0.35,
    "target_peak": 0.95
}
Параметры:
Поле
Описание
sample_rate 
Частота дискретизации аудио (Гц)
oversample_factor
Коэффициент oversampling
duration 
Длительность синтеза (сек). 0.0 = автоматическая
fade_in
Длительность атаки (сек)
fade_out 
Длительность затухания (сек)
loss_factor
Коэффициент потерь для расчёта T60
nonlinear_damping
Нелинейное затухание
tail_factor 
Множитель естественной длительности
envelope_threshold
Порог для определения длительности
boost_amplitude 
Ускорение атаки (0 = выключено)
slowdown_factor 
Замедление затухания в конце
impact_position
Позиция удара (0-1)
impact_force 
Сила удара
max_frequency
Максимальная частота для синтеза
use_acoustic_model
Включить акустическую модель
use_mount_model
Включить модель крепления
reverb_mix
Уровень реверберации (0-1)
use_rms_normalization 
Включить RMS нормализацию
target_peak
Целевой пиковый уровень
 
6. Запуск оптимизатора
6.1. Подготовка
1. 
Создайте рабочую папку (например, bells/)
2. 
Поместите в неё файлы конфигурации:
◦ 
targets.json
◦ 
evolution.json
◦ 
(опционально) synthesis.json — для синтеза после оптимизации
6.2. Запуск
bash
cd bells
..\build\bell_opt.exe
6.3. Что происходит
1. 
Программа загружает targets.json и evolution.json
2. 
Если есть best_geometry.json — продолжает оптимизацию с сохранённого результата
3. 
Запускает эволюционный поиск оптимальной формы
4. 
Каждое поколение показывает:
◦ 
Номер поколения
◦ 
Лучшую ошибку
◦ 
Параметры геометрии
◦ 
Частоты всех тонов
5. 
Сохраняет лучший результат в best_geometry.json
6.4. Пример вывода
text
============================================================
=========  BELL OPTIMIZER - Orthodox Bell Tuning  ==========
============================================================

  Material: E=11000.0 MPa, nu=0.34, rho=8.60e-09 t/mm^3
  Target Prime: 432.0 Hz
  Target ratios (relative to Prime):
    Hum     = 0.50  (216.0 Hz)
    Tertia  = 1.25  (540.0 Hz)
    Quinta  = 1.50  (648.0 Hz)
    Nominal = 2.00  (864.0 Hz)
  Population size: 30
  Generations: 30
============================================================

 Gen        Error          H         Rm      t_top      t_bot          Hum        Prime       Tertia       Quinta      Nominal
------------------------------------------------------------------------------------------------------------------------
   0       0.1889     1322.1      710.6       44.9       92.9        231.0        432.0        523.1        673.9        875.9
   1       0.1768     1298.0      697.7       44.0       91.2        235.3        440.0        532.8        686.4        892.1
 
7. Запуск синтезатора
7.1. Подготовка
После работы оптимизатора в папке bells/ появляется файл bell.dat с частотами мод.
7.2. Запуск
bash
cd bells
..\build\bell_synth.exe bell.dat
Или с указанием пути к файлу:
bash
..\build\bell_synth.exe best_bell.dat
7.3. Результаты
Программа создаёт:
Файл
Описание
bell_XXXhz.wav 
Звуковой файл колокола
envelopes_all.csv
Огибающие всех мод
spectrum_XXXhz.csv 
Спектр в CSV формате
Где XXX — частота Prime в герцах.
7.4. Пример вывода
text
======================================================================
🔔 СИНТЕЗ КОЛОКОЛА ИЗ .dat (С FFT)
======================================================================

✅ Загружена конфигурация из synthesis.json

📋 ПАРАМЕТРЫ:
   Основная частота (Prime): 432.0 Гц
   Количество мод: 9
   Сила удара: 1.0
   Коэффициент потерь: 2.00e-04
   Длительность: 4.0 с
   Fade in: 0.010 с, Fade out: 0.200 с

🎵 СПЕКТР МОД (ТОЛЬКО ИЗ .dat):
   ----------------------------------------------------------------------
   № | Частота(Гц) | Амплитуда | T60(с) | Название
   ----------------------------------------------------------------------
    1 |      216.00 |  0.800000 |   5.20  | Hum
    2 |      432.00 |  1.000000 |   2.60  | Prime
    3 |      540.00 |  0.600000 |   2.08  | Tertia
    4 |      648.00 |  0.400000 |   1.73  | Quinta
    5 |      864.00 |  0.200000 |   1.30  | Nominal
   ----------------------------------------------------------------------

   Сохранён: bell_432hz.wav
   Размер: 1536 KB

======================================================================
✅ ГОТОВО: bell_432hz.wav
   Источник: bell.dat
   Мод: 9
   Длительность: 4.0 с
   CSV спектр: spectrum_432hz.csv
   CSV огибающих: envelopes_all.csv
======================================================================
 
8. Частоты православных колоколов
 
Тон
Коэффициент
Частота, Гц (при Prime = 432)
Описание
Hum (гул)
0.50 
216.0
Самый низкий тон, основа звучания
Prime (основной)
1.00 
432.0
Главный тон колокола
Tertia (терция)
1.25 
540.0
Терция, придаёт мягкость
Quinta (квинта)
1.50 
648.0
Квинта, добавляет яркость
Nominal (номинальный)
2.00
864.0
Верхний тон, завершает спектр
Важно: Prime всегда = 1.0 и не включается в массив freq_ratios. Массив содержит только соотношения остальных тонов: [Hum/Prime, Tertia/Prime, Quinta/Prime, Nominal/Prime].
 
9. Типовые ошибки и их решение
Ошибка
Решение
Файл targets.json не найден 
Создать targets.json в рабочей папке
Файл evolution.json не найден 
Создать evolution.json в рабочей папке
CalculiX failed 
Установить CalculiX или проверить путь в PATH
No modes found 
Проверить правильность bell.inp или увеличить количество мод
No valid candidates
Увеличить population_size или ослабить ограничения
SAVE FAILED
Проверить права на запись или закрыть best_geometry.json в редакторе
 
10. Продвинутые настройки
10.1. Изменение целевой частоты
В targets.json измените prime_freq_target:
json
"prime_freq_target": 440.0   // Для 440 Гц
10.2. Изменение соотношений тонов
Для музыкального колокола:
json
"freq_ratios": [0.5, 1.0, 1.5, 2.0],
"amp_ratios": [0.5, 0.8, 0.4, 0.2]
10.3. Более точная оптимизация
В evolution.json:
json
"population_size": 50,
"generations": 50,
"num_elite": 8
10.4. Более быстрая оптимизация
В evolution.json:
json
"population_size": 15,
"generations": 10,
"num_elite": 3
10.5. Включение эффектов в синтезаторе
В synthesis.json:
json
"use_acoustic_model": true,
"use_acoustic_coupling": true,
"use_helmholtz_resonance": true,
"use_mount_model": true,
"reverb_mix": 0.12
 
11. Командный файл для запуска
Создайте optimize.cmd в корневой папке:
batch
@echo off
cd bells
..\build\bell_opt.exe
pause
Для синтеза звука:
batch
@echo off
cd bells
..\build\bell_synth.exe bell.dat
pause
 
12. Заключение
Bell Optimizer позволяет автоматически находить оптимальную геометрию колокола с заданными акустическими характеристиками. Программа полностью управляется через JSON-файлы, что делает её гибкой и настраиваемой под любые задачи.
