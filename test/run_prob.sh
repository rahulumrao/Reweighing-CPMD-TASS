#!/bin/bash
# This bash script will create unbias probability which can be used in WHAM to compute free energy
# For multiple umbrella window one has to run this script in each umbrella folder to generate 
# individual umbrella probabilities
cd cv_file/umb_2.2/ 	# specify the cv file location/folder
  ../../bin/Probability_analysis.x ,\
 -T0 300                 ,\
 -T 1000                 ,\
 -bias_fact 1500         ,\
 -tmin 5000              ,\
 -ncv 5                  ,\
 -UCV 1                  ,\
 -MTD n                  ,\
 -MCV 0                  ,\
 -tool probT		 ,\
 -Prob_nD 1              ,\
 -CV_num 1               ,\
 -grid 1.0 5.0 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 1.0 6.0 0.05 ,\
 -pfrqMD 1 ,\
 -dtMTD 200

# RAHUL VERMA
# DEPARTMENT OF CHEMISTRY
# IIT KANPUR, INDIA
# Email : vrahul@iitk.ac.in
