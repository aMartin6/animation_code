import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.path import Path
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial.distance import pdist, squareform
import numpy as np
import colloid_plot_library as cpl
import data_importerDM as di
import sys
from math import *

def calc_energy(x_p, y_p, x_d, y_d, spring, disorder_mag,
               particle_rad, disorder_rad):

    '''
    Use the scipy library on the given numpy arrays for 
    the positions of the particles
    
    x_p:      numpy array containing all x-positions of particles
    y_p:      numpy array containing all y-positions of particles
    x_d:      numpy array containing all x-positions of disorder
    y_d:      numpy array containing all y-positions of disorder
    spring:   spring constant of particle-particle reactions
    disorder_mag: magnitude of potential for disorder
    '''
    energy = []
    #Calculate energy from neighboring particles
    for i in range(len(x_p)):
        disorder_energy = 0
        particle_energy = 0
        for p in range(len(x_p)):
            if (i != p):
                distance = np.sqrt((x_p[i] - x_p[p])**2 + ((y_p[i] - y_p[p])**2))
            else:
                distance = 100
            if (distance < 2*particle_rad):
                particle_energy += 0.5 * spring * (particle_rad-distance)**2
                
    #Calculate energy from overlapping disorder
        for p in range(len(x_d)):
                distance = np.sqrt((x_p[i] - x_d[p])**2 + ((y_p[i] - y_d[p])**2))
                if (distance < disorder_rad):
                    disorder_energy += abs(disorder_mag * (disorder_rad-distance)**2)

        energy.append(disorder_energy + particle_energy)
    return energy

if __name__ == "__main__":
    #all possibilities
    data_types = [0,1,2] #["smtest", "ascii", "binary"]

    #the one we will use
    data_type = data_types[0]

    if data_type == 0:
        print("Reading directly from smtest (binary)")
        print("Writing velocity_data/XV*npy files")

    elif data_type == 1:
        print("Reading from velocity_data/XV* files (ascii)")
        print("Writing velocity_data/XV*npy files")

    elif data_type == 2:
        print("Reading velocity_data/XV*npy files (binary)")


    #--------------------------------------------------------------
    #get data for initial frame, 
    #---------------------------------------------------------------
    inputfile = "Pa0"
    
    (Sx, Sy, radius, maxtime, writemovietime ) = cpl.get_input_data(inputfile)

    out_file1="Energies.txt"
    out_file2="../../AverageEnergies.txt"

    #---------------------------
    #get and parse data
    #---------------------------

    datafile_prefix = "velocity_data/XV_data_t="
    pin_file = "pin_array.dat"
    try: 
        pin_data = di.get_data(pin_file,5,sep=" ")
    except:
        print("No pinning data in expected format")
        sys.exit()       
    xd = pin_data[1]
    yd = pin_data[2]
    pin_rad = pin_data[3][0]
    pin_mag = pin_data[4][0]
    total_energy = 0
    mintime = 500
    EnergyArray = []

    for time in range(mintime,maxtime,5000):
        id,type,xp,yp = cpl.get_and_parse_data(data_type, time,
                                               movie_type="cmovie")
        #Calculate particle energies
        kspring = 50
        energy_data=calc_energy(xp,yp,xd,yd,kspring,pin_mag,radius,pin_rad)
        EnergyArray.append(sum(energy_data))
        total_energy += sum(energy_data)
        with open(out_file1, "a") as myfile:
            myfile.write("%f" %time)
            myfile.write(" ")
            myfile.write("%f" %total_energy)
            myfile.write("\n")
    
    Average = (total_energy * 5000) / (maxtime - 500)
    with open(out_file2, "a") as myfile:
        myfile.write("%f" %Average)
    
    SumValues = 0
    Count = 0
    for time in range(mintime,maxtime,5000):
        SumValues += abs(EnergyArray[Count] - Average)
        Count += 1
        
    Stddev1 = sqrt((SumValues*5000)/(maxtime-500))
    with open(out_file2, "a") as myfile:
        myfile.write(" ")
        myfile.write("%f" %Stddev1)
        myfile.write("\n")

    sys.exit()
