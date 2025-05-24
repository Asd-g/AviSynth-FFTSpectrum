## Description

A filter that displays the FFT frequency spectrum of a given clip. Supposedly useful for determining original resolution of upscaled anime content.

This is [a port of the VapourSynth plugin FFTSpectrum](https://github.com/Beatrice-Raws/FFTSpectrum).

[fftw3](https://github.com/FFTW/fftw3) is used.

### Requirements:

- AviSynth 2.60 / AviSynth+ 3.4 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

### Usage:

```
FFTSpectrum (clip, bool "grid", int "opt")
```

### Parameters:

- clip<br>
    A clip to process. It must be in YUV 8-bit planar format.

- grid<br>
    Whether a grid with origin at the center of the image and spacing of 100 pixels should be drawn over the resulting spectrum.<br>
    Default: False.

- opt<br>
    Sets which cpu optimizations to use.<br>
    -1: Auto-detect.<br>
    0: Use C++ code.<br>
    1: Use SSE2 code.<br>
    2: Use AVX2 code.<br>
    3: Use AVX-512 code.<br>
    Default: 1.

### Building:

```
Requirements:
- Git
- C++20 compiler
- CMake >= 3.16
- Ninja
- FFTW
```

|    Option   |        Description       | Default value |
|:-----------:|:------------------------:|:-------------:|
| STATIC_FFTW | Link against static FFTW |      OFF      |


```
git clone https://github.com/Asd-g/AviSynth-FFTSpectrum
cd AviSynth-FFTSpectrum
cmake -B build -G Ninja -DCMAKE_PREFIX_PATH=<path_to_the_fftwf_installation>
ninja -C build
```
