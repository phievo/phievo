import numpy
import pylab
import config
from popen2 import Popen3
import linecache, copy
from classes_eds2 import draw_Network



exec('from '+config.name_pretty_graph+' import *')
exec('from '+config.name_deriv2+' import compile_and_integrate')



from palette import *

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
    pylab.clf()
    linecache.clearcache()
    ax = pylab.gca()
    for indexstep in range(npoints):
        line = linecache.getline(Name, indexstep * nstep / npoints)  # reads a line
        s = line.split()  # splits each line and puts it in an array
        for index in s:
            for ngene in range(size):
                try:
                    result[ngene, indexstep] = eval(s[
                        1 + ncell * size + ngene])  # with previous convention s[1+ncell*size+ngene] contains P_ngene_Cell_ncell
                except:
                    print(s)
                    print(size, ncell, ngene, indexstep)
    for i in range(size):
        style = '-'
        if i in list_output:
            style = '-'
        if i in list_toplot:
            pylab.plot(result[i], style, color=colors[i], lw=4.0)
        # pylab.ylim(ymax=1.2)
        if (len(list_time) > 0) and (list_time[i] >= 0):
            position = float(list_time[i])
            pylab.text(position - 15, result[i][list_time[i]] + 0.1, str(i), fontsize=20, fontweight='bold',
                       color=colors[i])
        legendkey.append("Species %i" % i)
    fontsize = 20
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

    pylab.xlabel('Time, Cell=' + str(ncell), fontsize=20, fontweight='bold')
    pylab.ylabel('Concentration', fontsize=20, fontweight='bold')

    pylab.legend(legendkey, loc='best')
    pylab.show()



def Plot_Profile(Name, ncelltot, size, nline, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins at time nline, size is the number of variables (proteins) in the cell, ncelltot the total number of cells in the embryo"""

    result = numpy.zeros((size, ncelltot),
                         dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i in cell j at time defined by nline
    colors = color_generate(size)
    pylab.clf()
    legendkey = []
    linestyle = ['-', '--', ':', '-.']
    cursor = 0  # cursor will be the index of the column
    linecache.clearcache()
    line = linecache.getline(Name, nline - 1)  # get the line we want to plot
    pylab.rc('axes', linewidth=2)
    ax = pylab.gca()
    s = line.split()
    if list_species==[]:
        list_toplot=range(size)
    else:
        list_toplot=list_species
    for index in s:  # takes each variable in the line
        if (cursor > 0):
            ncell = int((
                            cursor - 1) / size)  #computes the corresponding cell : column 1 to size corresponds to cell 0, size+1 to 2*size corresponds to cell 1, ...
            ngene = cursor - 1 - size * ncell  #computes the corresponding gene
            try:
                result[ngene, ncell] = eval(index)  #assigns it to the table result
                #if (result[ngene,ncell]<0.1):
                #result[ngene,ncell]=10*eval(index)#Be careful of this when plotting a profile (discontinuity from 1->0.1

                #if (ngene==1) or (ngene==3):
                #result[ngene,ncell]=0.4*eval(index)

            except:
                print(ngene, ncell, index)
                print("Error in Plot_Profile")
        cursor += 1
    n_position = numpy.zeros((ncelltot), dtype=float)
    # print "Tab filled"
    for i in range(size):
        style = '--'
        #print i
        if i in list_output:
            style = '-'
        if i in list_toplot:
            pylab.plot(result[i], style, color=colors[i], lw=4.0)
        #pylab.ylim(ymax=ylimmax)
        legendkey.append("Species %i" % i)
        #if (len(list_AP)>0) and (list_AP[i]>=0):
        #    position=float(list_AP[i])+n_position[list_AP[i]]
        #    if (position==0):
        #        position=0.5+n_position[0]

        #    pylab.text(position,result[i][list_AP[i]]+0.1,str(i),fontsize=20,fontweight='bold',color=colors[i])
        #    n_position[list_AP[i]]=n_position[list_AP[i]]+0.7
    fontsize = 20
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

    pylab.xlabel('Position', fontsize=20, fontweight='bold')
    pylab.ylabel('Concentration', fontsize=20, fontweight='bold')
    pylab.legend(legendkey, loc=position)
    #ax.set_yscale('log')
    pylab.show()
    


def Plot_pMHC(Name, ncelltot, size, ntau, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins at time nline, size is the number of variables (proteins) in the cell, ncelltot the total number of cells in the embryo"""

    result = numpy.zeros((size, ncelltot,ntau),dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i in cell j at time defined by nline
    colors = color_generate(size)
    pylab.clf()
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
                    result[ngene, ncell, idtau] = eval(index)  #assigns it to the table result

            except:
                print(ngene, ncell, idtau)
                print(eval(index))
                print("Error in Plot_Profile")
        cursor += 1
    n_position = numpy.zeros((ncelltot), dtype=float)
    # print "Tab filled"
    #for k in xrange(ntau):
    for k in range(2):
        pylab.subplot(1,2,k+1)
        #pylab.subplot(1,ntau,k+1)
        pylab.rc('axes', linewidth=2)
        ax = pylab.gca()
        ax.set_yscale('log')
        pylab.ylim([0.01,100000])
        for i in range(size):
            style = '--'
            #print i
            if i in list_output:
                style = '-'
            if i in list_toplot:
                #pylab.plot(result[i,:,k], style, color=colors[i], lw=4.0)
                pylab.plot(result[i,k,:], style, color=colors[i], lw=4.0)
                legendkey.append("Species %i" % i)
        fontsize = 20
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')

        pylab.xlabel('Position', fontsize=20, fontweight='bold')
        pylab.ylabel('Concentration', fontsize=20, fontweight='bold')
        pylab.xlim(0,20)
    pylab.legend(legendkey,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    pylab.show()
        

def net_test_plot(net, prmt, namefolder, generation, my_ntries=1, index_cell=-1,list_species=[]):
    """Given network, prmt dictionary, run CCode and send the network graph and time history tonamefolder, with suitable name built from generation number
    """
    print("Plotting")
    P = pretty_graph(net)
    P.write_dot(namefolder + "/GraphGen%i.dot" % generation)
    pylab.figure()
    # comment from here not to display results
    my_prmt = copy.copy(prmt)
    my_prmt['ntries'] = my_ntries
    nnetwork = 1000  # dummy network number, so not to overwrite others in /Workplace
    c_datafile = namefolder + "/Buffer"  # name used by CCode for output
    compile_and_integrate(net, my_prmt, nnetwork, 1, net.Cseed)
    for i in range(my_ntries):
        child = Popen3("mv Buffer%i " % i + namefolder)
        child.wait()
    size = len(net.list_types['Species'])
    ncelltot = prmt['ncelltot']
    nstep = prmt['nstep']
    list_output = []
    for output in net.list_types['Output']:
        list_output.append(output.int_id())

    if 'discrete' in prmt:
        nstep = size * nstep
    if 'pMHC' in net.list_types: 
        #quick and dirty way to check if this is a problem with pMHC, in the future we might want specific plot_data for individual problems
        pMHC_size = len(net.list_types['pMHC'])
        ntau=len(prmt['tau_off'])
        ncelltot=(size+pMHC_size+1)*ntau
        name_c = c_datafile + str(i)
        Plot_pMHC(name_c, ncelltot, size+pMHC_size+1, ntau,list_species=list_species,list_output=list_output)
        return;

    if (ncelltot == 1) or (index_cell > 0):
        for i in range(my_ntries):
            name_c = c_datafile + str(i)
            pylab.figure()
            Plot_Data(name_c, max(0, index_cell), size, nstep, list_species=list_species,list_output=list_output)
            return;

    for i in range(my_ntries):
            name_c = c_datafile + str(i)
            pylab.figure()
            Plot_Profile(name_c, ncelltot, size, nstep, list_species=list_species,list_output=list_output)


def plot_linear(net, k, prmt, nnetwork, namefolder):
    compile_and_integrate(net, prmt, nnetwork, 1)
    size = len(net.list_types['Species'])
    ncell = prmt['ncelltot']
    nstep = prmt['nstep']
    
    for i in range(prmt['ntries']):
        Plot_Profile("Buffer%i" % i, ncell, size, nstep)
        pylab.savefig(namefolder + "/Step%iProfile%i" % (k, i))


# ############### print instructions ########################################

if __name__ == '__main__':
    print("to plot a temporal  profile from data file 'Foo', type  Plot_Data('Foo',ncell,size,nstep) where ncell is the index of the cell you want to plot the data from, size is the number of species (i.e. variables per cell) and nstep the last time point")
    print("to plot a spatial  profile from data file 'Foo', type  Plot_Profile('Foo',ncell,size,nline) where ncell is the number of cells in the 'embryo', size is the number of species (i.e. variables per cell) and nline the index of the line you want to plot ")
