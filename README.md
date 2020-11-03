# Description

A filter that displays the FFT frequency spectrum of a given clip. Supposedly useful for determining original resolution of upscaled anime content.

This is [a port of the VapourSynth plugin FFTSpectrum](https://github.com/Beatrice-Raws/FFTSpectrum).

[fftw3](https://github.com/FFTW/fftw3) is used.

libfftw3f-3.dll or fftw3.dll required.

# Usage

```
FFTSpectrum (clip, bool "grid")
```

## Parameters:

- clip\
    A clip to process. It must be in YUV 8-bit planar format.
    
- grid\
    Whether a grid with origin at the center of the image and spacing of 100 pixels should be drawn over the resulting spectrum.\
    Default: False.
