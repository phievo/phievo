from phievo import __silent__,__verbose__
if __verbose__:
    print('Execute plotdata.py')

from phievo.initialization_code import *
import numpy
from matplotlib import pyplot as plt
import subprocess
import linecache, copy

from phievo.Networks.palette import *

""" Various utilities to plot data"""

"""Data should be organized like this :

A row corresponds to a time step. The first number of the row is the index of the time step. Then, concentrations of proteins are displayed in order. A typical line of the File to plot is organized like this :

Time_step P_0_Cell_0 ... P_N_Cell_0 P_0_Cell_1 ... P_N_Cell_1 P_0_Cell_2 ... P_N_Cell_2 ...

where P_i_Cell_j is the concentration of protein i in cell j at the time step Time_Step. You therefore need to give the number of proteins per cell as an argument of the fonction to now which  float corresponds to which cells and proteins
"""
ylimmax = 1.2

def Plot_Data(Name, ncell, size, nstep, list_species=[],list_time=[], list_output=[], npoints=500):
    """ Plot the concentrations of all proteins in one cell vs time
	ncell = which cell [0.,,)
	size = number of proteins
	nstep = number of time steps (beginning from 0) to plot
	npoints= number of points for each curve
	"""
    colors = color_generate(size)
    data = open(Name, "r")
    legendkey = []
    if list_species==[]:
        list_toplot=range(size)
    else:
        list_toplot=list_species
    result = numpy.zeros((size, npoints),
                         dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i at time j in the cell ncell given as an argument
    linecache.clearcache()
    fig = plt.figure()
    ax = fig.gca()
    for indexstep in range(npoints):
        line = linecache.getline(Name, int(indexstep * nstep / npoints))  # reads a line
        s = line.split()  # splits each line and puts it in an array
        for index in s:
            for ngene in range(size):
                try:
                    result[ngene, indexstep] = float(s[1 + ncell * size + ngene])
                    # with previous convention s[1+ncell*size+ngene] contains P_ngene_Cell_ncell
                except Exception:
                    msg = "Error in Plot_data, on string {}, unable to reach index {}"
                    msg = msg.format(s,1 + ncell * size + ngene)
                    display_error(msg)
    for i in range(size):
        style = '-'
        if i in list_output:
            style = '-'
        if i in list_toplot:
            print(len(result[i]))
            plt.plot(result[i], style, color=colors[i], lw=4.0)
        # plt.ylim(ymax=1.2)
        if (len(list_time) > 0) and (list_time[i] >= 0):
            position = float(list_time[i])
            ax.text(position - 15, result[i][list_time[i]] + 0.1, str(i), fontsize=20, fontweight='bold',
                       color=colors[i])
        legendkey.append("Species %i" % i)
    fontsize = 20
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label1.set_fontsize(fontsize)
    #     tick.label1.set_fontweight('bold')
    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label1.set_fontsize(fontsize)
    #     tick.label1.set_fontweight('bold')
    ax.set_xlabel('Time, Cell=' + str(ncell), fontsize=20, fontweight='bold')
    ax.set_ylabel('Concentration', fontsize=20, fontweight='bold')
    ax.legend(legendkey, loc='best')
    plt.show()
    return fig

def Plot_Profile(Name, ncelltot, size, nline, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins at time nline, size is the number of variables (proteins) in the cell, ncelltot the total number of cells in the embryo"""

    result = numpy.zeros((size, ncelltot),dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i in cell j at time defined by nline
    colors = color_generate(size)
    legendkey = []
    linestyle = ['-', '--', ':', '-.']
    cursor = 0  # cursor will be the index of the column
    linecache.clearcache()
    line = linecache.getline(Name, nline - 1)  # get the line we want to plot
    plt.rc('axes', linewidth=2)
    fig = plt.figure()
    ax = fig.gca()
    s = line.split()
    if list_species==[]:
        list_toplot=range(size)
    else:
        list_toplot=list_species
    list_toplot=range(size)
    for index in s:  # takes each variable in the line
        if (cursor > 0):
            ncell = int((cursor - 1) / size)  #computes the corresponding cell : column 1 to size corresponds to cell 0, size+1 to 2*size corresponds to cell 1, ...
            ngene = cursor - 1 - size * ncell  #computes the corresponding gene
            try:
                result[ngene, ncell] = float(index)  #assigns it to the table result

            except Exception:
                display_error(str(ngene)+str(ncell)+str(index)+"Error in Plot_Profile")
        cursor += 1
    n_position = numpy.zeros((ncelltot), dtype=float)
    # print "Tab filled"
    for i in range(size):
        style = '--'
        #print i
        if i in list_output:
            style = '-'
        if i in list_toplot:
            ax.plot(result[i], style, color=colors[i], lw=4.0)
        legendkey.append("Species %i" % i)
    fontsize = 20
    # for tick in ax.xaxis.get_major_ticks():
    #     tick.label1.set_fontsize(fontsize)
    #     tick.label1.set_fontweight('bold')
    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label1.set_fontsize(fontsize)
    #     tick.label1.set_fontweight('bold')
    ax.set_xlabel('Position', fontsize=20, fontweight='bold')
    ax.set_ylabel('Concentration', fontsize=20, fontweight='bold')
    #plt.legend(loc=0)
    ax.legend(legendkey)
    #ax.set_yscale('log')
    plt.show()
    return fig



def Plot_pMHC(Name, ncelltot, size, ntau, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins at time nline, size is the number of variables (proteins) in the cell, ncelltot the total number of cells in the embryo"""

    result = numpy.zeros((size, ncelltot,ntau),dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i in cell j at time defined by nline
    colors = color_generate(size)
    plt.clf()
    plt.hold(True)
    legendkey = []
    linestyle = ['-', '--', ':', '-.']
    cursor = 0  # cursor will be the index of the column
    linecache.clearcache()
    line = linecache.getline(Name, 1)  # get the line we want to plot
    s = line.split()
    if list_species==[]:
        list_toplot=range(size)
    else:
        list_toplot=list_species
    for index in s:  # takes each variable in the line
        if (cursor > 0):
            ncell = int(( cursor - 1) / (ntau*size) ) #computes the index of corresponding initial condition : column 1 to size*tau corresponds to CI 0, size*tau+1 to 2*size*tau corresponds to CI 2
            idtau=int((cursor - 1 - size * ntau*ncell)/size)#computes the corresponding tau
            ngene=cursor-1-size * ntau*ncell-size*idtau
            #ngene = int((cursor - 1 - size * ntau*ncell)/ntau)  #computes the corresponding species
            #idtau=cursor-1-size * ntau*ncell-ngene*ntau #computes the corresponding tau
            #print ngene,ncell,idtau
            try:
                    result[ngene, ncell, idtau] = float(index)  #assigns it to the table result
            except Exception:
                display_error(str(ngene)+' '+str(ncell)+' '+str(idtau)+' '+str(index)+"\nError in Plot_Profile")

        cursor += 1
    n_position = numpy.zeros((ncelltot), dtype=float)
    # print "Tab filled"
    for k in range(ntau):
        plt.subplot(1,ntau,k+1)
        plt.rc('axes', linewidth=2)
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim([0.01,100000])
        for i in range(size):
            style = '--'
            #print i
            if i in list_output:
                style = '-'
            if i in list_toplot:
                plt.plot(result[i,:,k], style, color=colors[i], lw=4.0)
                legendkey.append("Species %i" % i)
        fontsize = 20
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')

        plt.xlabel('Position', fontsize=20, fontweight='bold')
        plt.ylabel('Concentration', fontsize=20, fontweight='bold')
        plt.xlim(0,20)
    plt.legend(legendkey,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


def net_test_plot(compile_and_integrate, net, prmt, namefolder, generation, my_ntries=1, index_cell=-1,list_species=[],):
    """Given network, prmt dictionary, run CCode and send the network graph and time history tonamefolder, with suitable name built from generation number
    """
    print("Plotting")
    net.draw()

    plt.figure()
    # comment from here not to display results
    my_prmt = copy.copy(prmt)
    my_prmt['ntries'] = my_ntries
    nnetwork = 1000  # dummy network number, so not to overwrite others in /Workplace
    c_datafile = namefolder + "/Buffer"  # name used by CCode for output
    compile_and_integrate(net, my_prmt, nnetwork, 1, net.Cseed)
    for i in range(my_ntries):
        subprocess.run(["mv","Buffer%d"%i,namefolder])
    size = len(net.list_types['Species'])
    ncelltot = prmt['ncelltot']
    nstep = prmt['nstep']
    list_output = []
    for output in net.list_types['Output']:
        list_output.append(output.int_id())

    if 'discrete' in prmt:
        nstep = size * nstep

    if (ncelltot == 1) or (index_cell > 0):
        for i in range(my_ntries):
            name_c = c_datafile + str(i)
            plt.figure()
            Plot_Data(name_c, max(0, index_cell), size, nstep, list_species=list_species,list_output=list_output)
            return;

    for i in range(my_ntries):
            name_c = c_datafile + str(i)
            plt.figure()
            Plot_Profile(name_c, ncelltot, size, nstep, list_species=list_species,list_output=list_output)

def plot_linear(compile_and_integrate, net, k, prmt, nnetwork, namefolder):
    compile_and_integrate(net, prmt, nnetwork, 1)
    size = len(net.list_types['Species'])
    ncell = prmt['ncelltot']
    nstep = prmt['nstep']

    for i in range(prmt['ntries']):
        Plot_Profile("Buffer%i" % i, ncell, size, nstep)
        plt.savefig(namefolder + "/Step%iProfile%i" % (k, i))


# ############### print instructions ########################################

if __name__ == '__main__':
    print("to plot a temporal  profile from data file 'Foo', type  Plot_Data('Foo',ncell,size,nstep) where ncell is the index of the cell you want to plot the data from, size is the number of species (i.e. variables per cell) and nstep the last time point")
    print("to plot a spatial  profile from data file 'Foo', type  Plot_Profile('Foo',ncell,size,nline) where ncell is the number of cells in the 'embryo', size is the number of species (i.e. variables per cell) and nline the index of the line you want to plot ")
