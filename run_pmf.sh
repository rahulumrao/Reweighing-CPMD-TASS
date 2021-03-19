#!/bin/bash

#gfortran  -fcheck=all Probability_analysis.F90 -o Probability_analysis.x
#for i in -0.8 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4
#for i in 2.4
#do
echo -e "\e[0;31m=============================================================================================\e[0m"
echo -e "\e[0;34m                      ENTERING THE DIRECTORY umb_"${i}   "\e[0m"
echo -e "\e[0;31m=============================================================================================\e[0m"
#cd cv_file/umb_${i}

  bin/Probability_analysis.x ,\
 -T0 300                 ,\
 -T 1000                 ,\
 -bias_fact 1500         ,\
 -tmin 5000              ,\
 -ncv 5                  ,\
 -UCV 1                  ,\
 -MTD n                  ,\
 -MCV 0                  ,\
 -tool pmf		 ,\
 -interpolate		 ,\
 -nr 14			 ,\
 -grid 1.5 4.5 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 1.0 6.0 0.05 ,\
#mv Pu_2D.dat PROB_${i}
#cp PROB_${i} ../../PROB
#cd ../
#done

