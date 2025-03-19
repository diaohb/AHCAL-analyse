# Analyse

Writen by Hongbin

## Compile

```
cd [Dir]
mkdir build && cd build
cmake ..
make
```

## Usage

* ***digitize.cxx*** is to digitize MC to beam test format file
* `digi [input file] [pedestal file] [dac file] [mip file] [spe file] [lowgain adc file] [sipm model file] [output file]`
* ***Calib.cxx*** is to reconstruct energy from beam test format file 
* `Calib [input file] [pedestal file] [dac file] [mip file] [output file]`

## Calibration

in calibration, there are some root file that calibrated from beam test file, may be used in analyse code