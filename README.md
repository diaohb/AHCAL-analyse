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
* ***mip.cxx*** is to calibrate mip spectrum form MC format file
* `mip [input file]`
* ***mip_data.cxx*** is to calibrate mip spectrum form beam test format file
* `mipdata [input file] [pedestal file]`
* ***pedestal.cxx*** is to calibrate pedestal. 
* `pedestal [input file] [is force trigger mode]`
* ***digitize.cxx*** is to digitize MC to beam test format file
* `digi [input file] [pedestal file] [dac file] [mip file] [spe file] [lowgain adc file] [sipm model file] [output file]`
* ***Calib.cxx*** is to reconstruct energy from beam test format file 
* `Calib [input file] [pedestal file] [dac file] [mip file] [output file]`

* ***analyse.cxx*** is to analyse data after reconstruction and get some infomation of it, like hit number, hit layer, Fractal dimension(FD)
* `analyse [input file]`

* ***disp.cxx*** is to display an event in 3D
* `disp [input file] [event number]`
* ***disp2d.cxx*** is to display an event in 2D
* `disp [input file] [event number] [layer number]`

## Calibration

in calibration, there are some root file that calibrated from beam test file, may be used in analyse code