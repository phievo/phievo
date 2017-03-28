import numpy
import pylab
from popen2 import Popen3
import linecache, copy
from phievo.initialization_code import *
from phievo.Networks.classes_eds2 import draw_Network



exec('from '+name_pretty_graph+' import *')
import phievo.Networks.deriv2 as deriv2



from palette import *

""" Functions to  plot results for the immune recognition problem"""

ylimmax = 1.2





def Plot_pMHC(Name, ncelltot, size, ntau,list_ligands, position='best', list_species=[],list_AP=[], list_output=[]):
    """ Plot profile of all proteins as function of ligand concentration (list_ligands) contain the list, for the differeent taus (number ntau)"""

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
    for k in range(ntau):
        pylab.subplot(1,ntau,k+1)
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
                pylab.plot(result[i,:,k], style, color=colors[i], lw=4.0)
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
        sublist=[0,4,8,12,16,19]
        sublist_ligands=[list_ligands[index] for index in sublist]
        pylab.xticks(sublist,sublist_ligands)
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
    
    pMHC_size = len(net.list_types['pMHC'])
    ntau=len(prmt['tau_off'])
    ncelltot=(size+pMHC_size+1)*ntau
    name_c = c_datafile + str(i)
    Plot_pMHC(name_c, ncelltot, size+pMHC_size+1, ntau,my_prmt['Ligands'],list_species=list_species,list_output=list_output)
    return;

  

