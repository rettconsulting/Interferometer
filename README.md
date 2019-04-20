# Interferometer

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**This is a demonstration of fourier transform frequency scanning interferomtery.**

You can use a swept frequency laser source to measure suface figure and homogeneity of transmisive optics like windows and lenses. This is a different technique than phase shifting interferometry which requires you to physically perturb the cavity length to creaete swept fringe images.

This is based on Leslie Deck of the Zygo Corporation's published work:
**"Fourier-transform phase-shifting interferometry"**
https://www.osapublishing.org/ao/abstract.cfm?uri=ao-42-13-2354

This technique is used in some of their optical surface measurement tools, such as this:
**"Zygo Verifire™ MST"**
https://www.zygo.com/?/met/interferometers/verifire/mst/

This was written in 2008 or so. This is not an example of excellent programming (I write much better software every day...) but is an demo implementation of the tricky physical and mathematical technique desribed in Deck's paper.

This has been updated to a Visual Studio 2017 C++ solution, you can clone and compile it easily.

This uses CImage++ (for visualization) and FFTW (for fourier transform) libraries:
**http://cimg.eu/**
**http://www.fftw.org/**

Copyright © Alex Martin - alex@rettc.com - http://www.rettc.com

![image](https://raw.githubusercontent.com/mrlucretius/Interferometer/master/interferometer-with-unwrapping.jpg "Interferometer in Action")
