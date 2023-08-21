#!/bin/bash

DIR="obj";

if [ ! -d $DIR ]
then
     mkdir $DIR;
else
	rm -r $DIR; mkdir $DIR;
fi

cd obj;

h5pcc -O3 -fPIC -c ../src/h5_wrapper.c;
h5pcc -O3 -shared -o libh5_wrapper.so h5_wrapper.o;

h5pcc -O3 -fPIC -c ../src/fourier_transform.c;
h5pcc -O3 -shared -o libfourier_transform.so fourier_transform.o;

h5pcc -O3 -fopenmp -fPIC -c ../src/ft_module.c;
h5pcc -O3 -fopenmp -shared -o libft_module.so ft_module.o;

h5pcc -O3 -fopenmp -L. -lh5_wrapper -lfourier_transform -lft_module -lm -o main;

mv main	../;
