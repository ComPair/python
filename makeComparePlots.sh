#!/bin/bash

if [ $# -gt 0 ]; then

	python compareFigureOfMerit.py -f $1/*Aeff* -yl "Effective Area (cm$^2$)"
	python compareFigureOfMerit.py -f $1/*AngRes* -yl "Angular Resolution"
	python compareFigureOfMerit.py -f $1/*EnRes* -yl "Energy Resolution"

else
    echo "Usage: $0 <path to files for making comparision>"
fi

