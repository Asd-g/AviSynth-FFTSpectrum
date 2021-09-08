## Description

A filter that displays the FFT frequency spectrum of a given clip. Supposedly useful for determining original resolution of upscaled anime content.

This is [a port of the VapourSynth plugin FFTSpectrum](https://github.com/Beatrice-Raws/FFTSpectrum).

[fftw3](https://github.com/FFTW/fftw3) is used.

### Requirements:

- AviSynth 2.60 / AviSynth+ 3.4 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

- libfftw3f-3.dll or fftw3.dll

### Usage:

```
FFTSpectrum (clip, bool "grid")
```

### Parameters:

- clip\
    A clip to process. It must be in YUV 8-bit planar format.
    
- grid\
    Whether a grid with origin at the center of the image and spacing of 100 pixels should be drawn over the resulting spectrum.\
    Default: False.
