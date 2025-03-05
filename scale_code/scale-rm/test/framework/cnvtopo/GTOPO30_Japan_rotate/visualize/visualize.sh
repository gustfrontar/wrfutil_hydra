#! /bin/bash -x

### Visalization ###
echo "+visualize by gpview"
rm -f dcl.pdf

# boundary conditions
gpview topo.pe\*.nc@topo --nocont --aspect=4 --wsn 2 --srange=10: || exit
convert -density 150 -rotate 90 +antialias dcl.pdf topo.png
rm -f dcl.pdf
