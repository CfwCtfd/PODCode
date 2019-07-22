#!/bin/bash
#
g++ -g -c -I$HOME/libcpp/include pod.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling file_read.cpp"
  exit
fi
#
g++ pod.o $HOME/libcpp/lib/linpack_d.o $HOME/libcpp/lib/blas1_d.o $HOME/libcpp/lib/blas0.o  -lm
#g++ svd_demo.o $HOME/libcpp/lib/linpack_d.o $HOME/libcpp/lib/blas1_d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_demo.o + linpack_d.o + blas1_d.o."
  exit
fi
#
rm pod.o
#
chmod ugo+x a.out
mv a.out ~/bin/podunac_c++

#
echo "Executable installed "
