#!/bin/bash

ifort mod1.f90 m_config.f90 saida.f90 langevin.f90 -O3 -o langevin 
#ifort mod1.f90 m_config.f90 saida.f90 langevin.f90  -O0 -g -traceback -check all -o langevin