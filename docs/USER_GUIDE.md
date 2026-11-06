# User guide

## 1. Build

On Windows with MSYS2 UCRT64 and CMake available:

```bat
build.cmd
```

The executable is `build/bell.exe`. The script configures, builds and runs the registered CTest suite. `test.cmd` reruns the suite without rebuilding; `test_convergence.cmd [all|mesh|thickness]` runs the two mesh-refinement diagnostics; `test_ccx_current_geometry.cmd` runs the optional external CalculiX comparison. For an application-only build, configure CMake with `-DBUILD_SMOKE_TESTS=OFF` and build target `bell`.

## 2. Commands

```text
bell optimize [--solver=ccx|internal] [--workdir <directory>]
bell synth [file.dat] [--workdir <directory>]
bell synth-geometry [geometry.json] [--ellipticity q] [--output file.wav] [--workdir <directory>]
bell --help
bell --version
```

The optional `--workdir` changes the process working directory before the command runs. Paths to configuration and input files are therefore resolved relative to that directory.

## 3. Optimization

Typical internal-solver workflow from the project root:

```bat
bell.cmd optimize --solver=internal --workdir bells
```

The optional CalculiX reference workflow is:

```bat
bell.cmd optimize --solver=ccx --workdir bells
```

The `bells/` directory contains geometry seeds, evolution settings and modal target files. The optimizer writes its best geometry and modal data there. The Western `western_five_partial` objective is the production objective. Chinese target files are experimental specifications and are not a production two-tone optimizer path.

## 4. Synthesis from modal data

```bat
bell.cmd synth best_bell.dat --workdir bells
```

The command writes a mono WAV and associated analysis output. For the bundled listening example, run:

```bat
make_audio.cmd
```

The helper places generated listening files under `audio/`; it does not change the source geometry or tuning configurations.

## 5. Synthesis from geometry

```bat
bell.cmd synth-geometry best_geometry.json --ellipticity 0.85 --output bell_elliptic.wav --workdir bells
```

This workflow evaluates the supplied geometry and synthesizes a listening result. The current ellipticity-aware geometry workflow still uses the documented reduced-order coordinate-morphing path for its legacy modal-frequency estimate unless the direct FEM evaluation option is explicitly selected in the relevant configuration/API. Do not interpret the morphing estimate as a fresh FEM solve on the morphed geometry.

## 6. Configuration and output level

- `evolution*.json` controls evolutionary search settings.
- `targets*.json` defines modal tuning targets.
- `synthesis.json` controls sample rate, oversampling, duration, losses, normalization and optional reverb.
- `best_geometry.json` stores the optimized geometry; `best_bell.dat` stores modal data for synthesis.

`target_peak` controls the final digital peak level. Peak normalization is a presentation-level operation; it does not represent calibrated sound pressure level. Reverb is a post-effect and is disabled by default. The generated WAV is mono.

## 7. Physical interpretation

The production optimizer's meridional generalized influence is a simplified excitation proxy. The Chinese localized-strike model is a separate diagnostic and is not yet used by the production optimizer. Likewise, current acoustic transfer and synthesis are reduced-order models, not calibrated predictions of sound pressure at a specified distance.

For mathematical assumptions and known limits, see `MATHEMATICAL_AND_ALGORITHMIC_FOUNDATIONS.md` and `VALIDATION.md`.

## Материал

Материал задаётся в `bells/evolution*.json` только именем предустановки:

```json
"material": "bell_bronze_reference"
```

Solver получает только механические свойства предустановки: `E`, `nu` и `rho`. Химический состав сплава не является параметром BellOptimizer и не вводится в конфигурацию.

Доступные предустановки в текущей версии:

- `bell_bronze_reference` — текущая производственная эталонная бронза проекта: `E = 11000 MPa`, `nu = 0.34`, `rho = 8.60e-9 tonne/mm^3`. Это инженерная эталонная модель, а не историческая сертификация состава древней бронзы.
- `bronze_aluminum_reference` — механическая эталонная модель алюминиевой бронзы, перенесённая из предоставленного CalculiX-файла `Bronze_aluminum_alloy.inp`: `E = 110000 MPa`, `nu = 0.30`, `rho = 7.40e-9 tonne/mm^3`.

Если `material` не указан, используется `bell_bronze_reference`. Неизвестное имя предустановки является ошибкой конфигурации.

CCX/CalculiX-файлы используются как внешний справочный источник механических параметров; их тепловые, электрические и прочие поля не импортируются, если они не используются физической моделью BellOptimizer.

Стереофонический прототип для `synth-clapper` включается флагом `--stereo`. Параметры `--listener-distance 3.0` и `--ear-spacing 0.18` задают расстояние до центра головы и расстояние между виртуальными ушами; `--head-yaw-deg` задаёт поворот линии ушей относительно положения «лицом к колоколу» (положительное значение смещает правое ухо дальше от колокола). Расчёт FEM, контакта языка и модального движения общий для двух каналов; комплексное азимутальное поле каждой моды вычисляется один раз, затем оценивается для двух ушей. Это пока геометрическое стерео без модели головы, ушных раковин, HRTF и помещения. При нулевом повороте головы и источнике строго перед слушателем геометрические пути до ушей равны, поэтому без модели головы каналы могут быть очень похожи. Ненулевой `--head-yaw-deg` создаёт геометрическую разницу путей.

Пример: `bell synth-clapper bells/best_geometry.json --output audio/bell_stereo.wav --stereo --listener-polar-deg 30 --listener-azimuth-deg 20 --head-yaw-deg 10`.

Для отдельной команды прослушивания внутреннего языка (`synth-clapper`) у `bell_bronze_reference` дополнительно задан справочный коэффициент структурных потерь `eta_loss = 1.8e-4`. Он получен пересчётом среднего коэффициента модального демпфирования колокольной бронзы из литературы (`zeta = 0.9e-4`, `eta_loss ≈ 2*zeta`), а не измерением конкретной плавки. У `bronze_aluminum_reference` коэффициент не задан: отсутствие данных не заменяется выдуманным значением. `--loss-factor` в `synth-clapper` служит только для явных исследований чувствительности. Это свойство пока не используется оптимизатором.

