import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.path import Path
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial.distance import pdist, squareform
import numpy as np
import colloid_plot_library as cpl
import data_importerDM as di
import sys

plt.rc('font',size=22)

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


def plot_energy(ax,energy,x_p,y_p,size):
    
    #vmin=1e-7
    #vmax=1e-3
    #norm1=matplotlib.colors.LogNorm(vmin,vmax)
    ax.scatter(x_p,y_p,s=size,zorder=10,c=energy)

    return

if __name__ == "__main__":

    #all possibilities
    data_types = [0,1,2] #["smtest", "ascii", "binary"]

    #the one we will use
    data_type = data_types[2]

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

    #plot at the value starttime
    starttime=45000 #

    SX = (Sx[1]-Sx[0])
    SY = (Sy[1]-Sy[0])

    size=30  #hard coded by what "looks good"
    #---------------------------
    #set up a 1x1 plot in a subroutine
    #---------------------------

    fig,ax1 = cpl.format_plot(Sx=Sx,Sy=Sy)

    #---------------------------
    #plot the pinning array
    #---------------------------

    cpl.plot_pins(ax1,size=size)
    out_file="particle_energies.png"

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

    #if data_type == 0 , this will process all of smtest
    id,type,xp,yp = cpl.get_and_parse_data(data_type,
                                           starttime,
                                           movie_type="cmovie")
    #Calculate particle energies
    kspring = 50
    energy_data=calc_energy(xp,yp,xd,yd,kspring,pin_mag,radius,pin_rad)

    #plot it
    plot_energy(ax1,energy_data,xp,yp,size)
    ########################################################

    #add a label such as (a)
    label = False

    if label:
        xt=0.05
        yt=0.9
        ax1.text(xt,yt,label,transform = ax1.transAxes)
    #plt.axes().set_aspect('equal')
    #ax.set_aspect_ratio

    fig.tight_layout()
    plt.savefig(out_file)
    
    sys.exit()
