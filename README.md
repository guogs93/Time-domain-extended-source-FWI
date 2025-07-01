This repository contains the source code for time-domain extended-source full-waveform inversion. 

The code was developed by Gaoshan Guo during his PhD and postdoc (2021.01~2024.12) at the CNRS–Géoazur laboratory, under the supervision of Stéphane Operto, and supported by the WIND consortium (https://www.geoazur.fr/WIND/bin/view).

We gratefully acknowledge the sponsorship of AkerBP, ExxonMobil, Petrobras, Shell, and SINOPEC.

A detailed README file and demonstration examples will be released upon acceptance of the associated manuscript.

Usage
This code is intended for researchers working in seismic imaging of OBS data. 

Compilation and running instructions depend on the computing environment—please refer to your corresponding Makefile.inc_* file.

To clean compiled objects and binaries, run: make clean && make all

If you use this code in your research or publication, please cite the following paper:

```bibtex
@article{Guo_2025_RIC,
  author = {Guo, Gaoshan and Operto, Stéphane},
  title = {Robust imaging of the crust from long-offset ocean-bottom seismometer data using time-domain extended-source full-waveform inversion: Application to the eastern Nankai Trough, Japan},
  journal = {submitted to GRL},
}
