#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 13:03:06 2019

Este programa lê arquivos CSV onde cada arquivo é um passo de tempo e cada 
linha é uma partícula e exporta em vtk o reultado. 

@author: gabriel
"""

from evtk.hl import pointsToVTK
import numpy as np
import csv, os, shutil
import configparser
from zipfile import ZipFile
import time
try:
    import progressbar
    pbar = True 
except:
    print('progressbar not available. Try one:')
    print('conda install -c conda-forge progressbar2')
    print('conda install -c conda-forge/label/gcc7 progressbar2')
    pbar = False

via = os.getcwd()

aux = True
folder = 'result'
i = 1

#Ler a quantidade de arquivos de settings.

config = configparser.ConfigParser()
config.read('settings.ini')
print('Reading settings.ini')
N = int(config['global']['N'].split()[0])
nimpre =  int(config['global']['nimpre'].split()[0])

#Checar se existe arquivo
while aux:
    try:
        os.mkdir(folder)
        aux = False
    except FileExistsError:
        print('Folder {} exists'.format(folder))
        opt = input('Overwrite? [y/n] ')
        if opt == 'y' or opt == 'Y':
            aux =False 
            shutil.rmtree(folder)
            os.mkdir(folder)
        else:
            folder = input('Enter new folder name:\n')
            if len(folder) == 0:  
                folder = 'result' + '_' + str(i) 
                aux = True
                i = i+1       

print('\nThe results will be saved at {}\n'.format(folder))


    
a = os.listdir('temp')
no_out_files = nimpre
if len(a)/2 < nimpre:
    no_out_files = int(len(a)/2)
    print("Propable incomplete execution. Processing {} files".format(len(a)/2))
    nimpre = int(len(a)/2)-1 
    input("nimpre = {}. OK?\n".format(nimpre))

tipo = np.zeros(N)
rsol = np.zeros(N)
j,k = 0,0

        
# Ler os arquivos e colocá-los em vetores

x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)


vx = np.zeros(N)
vy = np.zeros(N)
# vz = np.zeros(N)
omega = np.full(N,np.nan)

# cx = np.zeros(N)
# cy = np.zeros(N)
# cz = np.zeros(N)

nID = np.linspace(1,N,N)
zip_positions = ZipFile(via+'/'+folder+'/positions.zip','a')
zip_velocities = ZipFile(via+'/'+folder+'/velocities.zip','a')    
    
    
print('Converting...')
position_list = zip_positions.namelist()
aux1 = 0
flag = True
if pbar:
    bar = progressbar.ProgressBar(max_value=nimpre)
for fnum in range(0,nimpre+1):
    with open('temp/position.csv.'+str(fnum),encoding='utf-8') as file_locus:
        csv_lector = csv.reader(file_locus,delimiter = ',')
        i = 0
        for linea in csv_lector:
            x[i] = linea[0]
            y[i] = linea[1]
            i += 1

    while flag:
        if any('position.csv.'+str(fnum+aux1) in s for s in position_list):
            aux1 += 1
        else:
            flag = False
    zip_positions.write(via+'/temp/position.csv.'+str(fnum), 'position.csv.'+str(fnum+aux1))
    # shutil.move(via+'/temp/position.csv.'+str(fnum),via+'/'+folder+'/position.csv.'+str(fnum)) 
    with open('temp/velocity.csv.'+str(fnum),encoding='utf-8') as file_velocitas:
        csv_lector = csv.reader(file_velocitas,delimiter = ',')
        i = 0
        for linea in csv_lector:
            vx[i] = linea[0]
            vy[i] = linea[1]
            i = i+1
            
    zip_velocities.write(via+'/temp/velocity.csv.'+str(fnum),'velocity.csv.'+str(fnum+aux1))
    # shutil.move(via+'/temp/velocity.csv.'+str(fnum),via+'/'+folder+'/velocity.csv.'+str(fnum))

    fin = 0
    grupo = 0
    
    ini = 0      
    fin = N
    xs = x[ini:fin]
    ys = y[ini:fin]
    zs = z[ini:fin]
    # vxs = vx[ini:fin]
    # vys = vy[ini:fin]
    # "Vx" : vxs, "Vy" : vys
    vel = (vx[ini:fin],vy[ini:fin],np.zeros(fin-ini))
    nIDs = nID[ini:fin]
    pointsToVTK(via+'/'+folder +'/grupo'+ str(grupo) + '_' +str(fnum+aux1), xs, ys, zs, data = {"Vel" : vel , "nID" : nIDs })               

    grupo += 1

    if pbar:
        bar.update(fnum)
    if fnum == nimpre+1:
        time.sleep(0.1)


zip_positions.close()
zip_velocities.close()

shutil.rmtree('temp')

with open("settings.txt","a") as settingstxt:
    settingstxt.write("[out_files]\n")
    settingstxt.write("out_files = {}".format(no_out_files))

try:
    shutil.move(via+'/settings.txt',via+'/'+folder+'/settings.txt') 
except:
    print('No settings file found!\n')

include_folder = True 
if os.path.isfile('.gitignore'):
    with open('.gitignore','a+') as file:
        for line in file:
            if line == folder:
                include_folder = False
        if include_folder:
            file.write(folder+'\n')
else:
    with open('.gitignore','a+') as file:
        file.write(folder+'\n')

 
    #with open('cell.csv.'+str(fnum),encoding='utf-8') as file_velocitas:
        #csv_lector = csv.reader(file_velocitas,delimiter = ',')
        #i = 0
        #for linea in csv_lector:
            #cx[i] = int(float(linea[0]))
            #cy[i] = int(float(linea[1]))
            
            #i = i+1
            
    #shutil.move(via+'/cell.csv.'+str(fnum),via+'/'+folder+'/cell.csv.'+str(fnum))  
