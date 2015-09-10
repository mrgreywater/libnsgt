# libnsgt 0.01
This is an implementation of the ['Nonstationary Gabor transform'](http://www.univie.ac.at/nonstatgab/) (Inversible Constant Q Transform) in C.
It supports forward and backward transform of the whole or streamed signal. It currently only does complex to complex transforms with double floating point precision. 

Usage
=====
The library consists of the files nsgt.h and nsgt.c. 
To make deployment as easy as possible, the project contains a sln file for [Visual Studio 2015](https://www.visualstudio.com/de-de/downloads/download-visual-studio-vs.aspx) and cmakelists.txt for [CLion](https://www.jetbrains.com/clion/) to compile and execute the added examples.

Usage of the library should be relatively straight forwards, and basic examples are found in the examples folder.

Examples
========
Executing the examples will yield something similar to this:
```
Test full nsgt
Initialization in 393.626 ms
Forward in 28.885ms, Backwards in 26.427ms, Error: -152.032db

Test nsgt stream (sliCQ)
Initialization in 132.662 ms
Forward in 5.999ms, Backwards in 4.990ms,
Forward in 5.374ms, Backwards in 4.906ms, Error: -152.613db
Forward in 5.254ms, Backwards in 4.995ms, Error: -152.779db
Forward in 5.255ms, Backwards in 4.871ms, Error: -153.877db
Forward in 5.717ms, Backwards in 4.921ms, Error: -152.243db
Forward in 5.305ms, Backwards in 5.449ms, Error: -152.238db
Forward in 5.392ms, Backwards in 5.009ms, Error: -152.481db
Forward in 5.422ms, Backwards in 4.861ms, Error: -152.391db
Forward in 5.220ms, Backwards in 4.856ms, Error: -150.750db
```
Also the following spectrogram output:
![Glockenspiel spectrogram, linear scaled](http://i.imgur.com/XqNUzdq.png)
![Glockenspiel spectrogram, decibel scaled](http://i.imgur.com/YWVcFLx.png)
![Glockenspiel spectrogram, phon scaled](http://i.imgur.com/GuxAzmg.png)
![Glockenspiel spectrogram, sone scaled](http://i.imgur.com/yjmF0Cj.png)
Information
===========
The papers and original Matlab implementation can be found here: http://www.univie.ac.at/nonstatgab/

A Python implementation can be found here: https://github.com/grrrr/nsgt

Dependencies
============

- FFTW3 for the fft transform ([Windows](http://www.fftw.org/install/windows.html) / [Linux](https://packages.debian.org/sid/libfftw3-dev))
- libsndfile for the examples ([Windows](http://www.mega-nerd.com/libsndfile/#Download) / [Linux](https://packages.debian.org/de/sid/libsndfile-dev))

Roadmap and contributing
=======
The following items need to be addressed:
- apart from double floating point precision, the library needs to support float and long double precision.
- real-to-complex forward and complex-to-real backwards tranforms needs to be implemented to save memory and time
- other fft libraries should be supported, one of them possibly be included to get rid of the dependencies
- other frequency scalings should be supported (currently only octave)
- support fftw wisdom and serialization

If you want to contribute, please fork this repository and make a pull request with your changes. It is greatly appreciated.

License
=======
[MIT License](LICENSE)

note that FFTW is licensed with [GPLv2](http://www.fftw.org/doc/License-and-Copyright.html)
