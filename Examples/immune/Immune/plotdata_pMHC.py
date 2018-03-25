import numpy
import pylab
#from popen2 import Popen3
import linecache, copy
from phievo.initialization_code import *
#from phievo.Networks.classes_eds2 import draw_Network
import matplotlib.gridspec as gridspec


#exec('from '+name_pretty_graph+' import *')
import phievo.Networks.deriv2 as deriv2

from phievo.AnalysisTools import palette



""" Functions to  plot results for the immune recognition problem"""

ylimmax = 1.2





def Plot_pMHC(Name, ncelltot, total_size,ntau,list_ligands, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins as function of ligand concentration (list_ligands) contain the list, for the differeent taus (number ntau)"""
    result = numpy.zeros((total_size, ncelltot,ntau),dtype=float)  # creates a table result to store data, results[i,j] will contain the concentration of protein i in cell j at time defined by nline
    colors = palette.color_generate(total_size)
    pylab.clf()
    legendkey = []
    linestyle = ['-', '--', ':', '-.']
    cursor = 0  # cursor will be the index of the column
    linecache.clearcache()
    print(Name)
    line = linecache.getline(Name, 1)  # get the line we want to plot
    s = line.split()
    if list_species==[]:
        list_toplot=range(total_size)
    else:
        list_toplot=list_species
    for index in s:  # takes each variable in the line
        if (cursor > 0):
            ncell = int(( cursor - 1) / (ntau*total_size) ) #computes the index of corresponding initial condition : column 1 to size*tau corresponds to CI 0, size*tau+1 to 2*size*tau corresponds to CI 2
            idtau=int((cursor - 1 - total_size * ntau*ncell)/total_size)#computes the corresponding tau
            ngene=cursor-1-total_size * ntau*ncell-total_size*idtau
            try:
                    result[ngene, ncell, idtau] = eval(index)  #assigns it to the table result

            except:
                print(ngene, ncell, idtau)
                print(eval(index))
                print("Error in Plot_Profile")
        cursor += 1
    n_position = numpy.zeros((ncelltot), dtype=float)
    # print "Tab filled"
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.2, right=0.8, wspace=0.0)
    for k in range(ntau):
        #pylab.subplot(1,ntau,k+1)
        #pylab.rc('axes', linewidth=2)
        ax = pylab.subplot(gs[0, k])
        #ax = pylab.gca()
        ax.set_yscale('log')
        pylab.ylim([0.01,1000])
        for i in range(total_size):
            style = '--'
            #print i
            if i in list_output:
                style = '-'
            if i in list_toplot:
                pylab.plot(result[i,:,k], style, color=colors[i], lw=4.0)
                legendkey.append("%i" % i)
        fontsize = 20
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            tick.label1.set_fontweight('bold')
        if (k==0):
            pylab.xlabel('Ligand #, 3s', fontsize=20, fontweight='bold')
        if (k==1):
            pylab.xlabel('Ligand #,  10s', fontsize=20, fontweight='bold')
        if (k==0):
            pylab.ylabel('Output', fontsize=20, fontweight='bold')
        else:
            ax.axes.get_yaxis().set_visible(False)
        pylab.xlim(0,19)
        sublist=[0,4,8,12,16]
        sublist_ligands=[list_ligands[index] for index in sublist]
        pylab.xticks(sublist,sublist_ligands)
    pylab.legend(legendkey,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    pylab.show()
        

