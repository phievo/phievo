from phievo.Networks.deriv2 import *

def all_params2C(net, prmt, print_buf,Cseed=0):
    """ collect all the numerical constants neelocalneig,diff,index_ligand,ded in C code, for a network, net, 
        prmt dictionary from initialziation, and boolean print_buf to control printing
    of time history by C codes 
    """
    hdr = []   # collect lines of output as list then join, speed issue

    # various sizes/lengths mostly from prmt dict
    size=len(net.dict_types['Species'])
    
    #for s in net.dict_types['Species']:
    #    s.print_node()
    pMHC_size = len(net.dict_types['pMHC'])
    #print len(net.dict_types['pMHC'])
    #print net.dict_types['pMHC']
    hdr.append("#define SIZE %i" %(size+pMHC_size+1) )
    #Index are organized in the following way
    #from 0 to size-1-> agonist, kinases,phosphatases and complex of the agonist cascades. Order is imposed by evolution
    # from size: self ligand
    # size+1->pmHC_size_1: cascade of the self ligands
    
    
    #hdr.append("#define NSTEP %i" %prmt['nstep'] )
    hdr.append("#define NNEIGHBOR %i" %prmt['nneighbor'] )
    hdr.append("#define NTRIES %i" %prmt['ntries'] )
    
    ##############################################################################
    
    
    # Added or modified for the pMHC-TCR simulation.
    
    # Related to concentrations and mechanics of the simulation.
    hdr.append("#define NTAU %i" %len(prmt['tau_off']) ) # Number of different dissociation times.
    hdr.append("#define NLIGANDS %i" %len(prmt['Ligands']) ) # Number of different ligand concentrations.
    hdr.append("#define NLIGANDS_SELF %i" %len(prmt['Ligand_self']))
    hdr.append("#define TAU_SELF %f" %prmt['tau_self'])
    
    hdr.append("#define NCELLTOT NTAU*NLIGANDS*NLIGANDS_SELF") # Total number of cells given by the number of dissociation times multiplied by the number of ligand concentrations.
    hdr.append("#define RECEPTOR %i" %prmt['Receptor_Conc'])
    
    hdr.append("#define NPARTITION_OUTPUT %i" %prmt['n_partition_output'])
    
    # specifically related to the adaptive integrator.
    hdr.append("#define MAX_TIME %f" %prmt['max_time'])
    hdr.append("#define TOL %.15f" %prmt['tolerance'])
    hdr.append("#define BELOW %f" %prmt['scalebelow'])
    hdr.append("#define ABOVE %f" %prmt['scaleabove'])
    hdr.append("#define TINY %.15f" %prmt['tiny'])
    
    # List of ID of species in the KPR cascade for init_history file.
    if 'pMHC' in net.dict_types:
        list_pMHC=net.dict_types['pMHC']
        hdr.append("#define pMHC_LENGTH %i" %len(list_pMHC) )
        hdr.append("static int pMHC_LIST[pMHC_LENGTH];")
    
    # List of all phosphatase or kinase (unphosphorylated).
    nP = 0
    if 'Phosphatase' in net.dict_types:
        P_list = net.dict_types['Phosphatase']
        for i in range(len(P_list)):
            if(P_list[i].n_phospho == 0):
                nP += 1    # all the phosphatase (not phosphorylated).
            
    nK = 0
    if 'Kinase' in net.dict_types:
        K_list = net.dict_types['Kinase']
        for i in range(len(K_list)):
            if(K_list[i].n_phospho == 0):
                if not K_list[i].isinstance('pMHC'):
                    nK += 1  
         
    S_list = net.dict_types['Species']
    counter = 0
    for i in range(len(S_list)):
        if (S_list[i].isinstance('Kinase') or S_list[i].isinstance('Phosphatase')):
            if not S_list[i].isinstance('pMHC'):
                out = net.graph.successors(S_list[i])
                for j in range(len(out)):
                    if out[j].isinstance('Initial_Concentration'):
                        counter += 1
    
    if counter != (nK + nP):
        print(nK+nP)
        print(counter)
        print("Error in counting the number of types of kinase and phosphatase.")
        
    hdr.append("#define KP_LENGTH %i" %(nK + nP) )
    hdr.append("static int KP_LIST[KP_LENGTH];")
    hdr.append("static double KP_CONC[SIZE];")
    
    # Lists containing the ligand concentrations and dissociation times of interest. Needed for the init_history file and compte_deriv_inC.
    # See initialization file for the entry in dictionary parameter corresponding to Ligands and tau_off.
    hdr.append("static int LIGAND_LIST[NLIGANDS];")
    hdr.append("static int LIGAND_SELF_LIST[NLIGANDS_SELF];")
    hdr.append("static double TAU_LIST[NTAU];")

    ##############################################################################
    if 'langevin_noise' in prmt:
        hdr.append("#define  CONCENTRATION_SCALE %f" %prmt['langevin_noise'] )
    else:
        hdr.append("#define  CONCENTRATION_SCALE 1.0" )
    
    # optional generic parameters for specific C subroutines as dict or list.
    # define NFREE_PRMT is flag in Ccode that free_prmt as list is being used
    if 'free_prmt' in prmt:            
        if isinstance(prmt['free_prmt'], dict):
            hdr.append("#define NFREE_PRMT 0" )
            for key,value in prmt['free_prmt'].items():
                hdr.append("#define %s %s" %( key.upper(), value) )                   
        elif isinstance(prmt['free_prmt'], list):
            hdr.append("#define NFREE_PRMT %i" %len(prmt['free_prmt']) )
            str_pp = ', '.join( [str(nn) for nn in prmt['free_prmt'] ] )
            hdr.append("static double free_prmt[] = {%s};" %str_pp )
        else:
            print('unexpected type of prmt[free_prmt] in deriv2.py.  Skipping', prmt['free_prmt'])
    else:
        hdr.append("#define NFREE_PRMT 0" )  # need define dummy length variable incase test in C
        hdr.append("static double free_prmt[] = {}; " )  # dummy defintion for C
    
    # IO stuff,
    trackin  = track_variable(net, 'Input')
    trackout = track_changing_variable(net, 'Output')  # this routine does not use n_put attribute of Species.
    #Ligand + Diffusible species
    tracklig = track_changing_variable(net, 'Ligand' )
    trackdiff = track_changing_variable(net, 'Diffusible' )
    
    # poor way to count the outputs (to specialized...)
    if 'Output' in net.dict_types:
        Outs = net.dict_types['Output']
        Out = Outs[0]  # at most one output...
        if Out.isinstance('pMHC'):
            hdr.append("#define    NOUTPUT %i" %(len(trackout)+1) )
        else:
            hdr.append("#define    NOUTPUT %i" %len(trackout) )

    hdr.append("#define    NINPUT %i" %len(trackin) )
    hdr.append("#define    NLIGAND %i" %len(tracklig) )
    hdr.append("#define    NDIFFUSIBLE %i" %len(trackdiff) )
    # misc other numbers
    if (Cseed==0):
        hdr.append("#define SEED %i" %(int(random.random()*1000000)) ) #seed for the C rand
    else:
        hdr.append("#define SEED %i" %Cseed )  
    
    hdr.append("#define PRINT_BUF %i" %print_buf )
    hdr.append("#define DT %f" %prmt['dt'] )

    # the mapping of input/output indices
    str_in = ', '.join( [str(nn) for nn in trackin] )
    hdr.append("static int trackin[] = {%s};" %str_in )
    str_out = ', '.join( [str(nn) for nn in trackout] )
    
    # very unmodular way of adding the ouput arising from the self ligands...
    if 'Output' in net.dict_types:
        Outs = net.dict_types['Output']
        Out = Outs[0]  # at most one output...
        if Out.isinstance('pMHC'):  # if the output is in the KPR cascade, need to add the other species for output.
            n = Out.n_phospho
            number_species = len(net.dict_types['Species'])
            str_out = str_out+", %i"%(n + number_species + 1)
            
    hdr.append("static int trackout[] = {%s};" %str_out )
    str_lig = ', '.join( [str(nn) for nn in tracklig] )
    hdr.append("static int tracklig[] = {%s};" %str_lig )

    str_diff = ', '.join( [str(nn) for nn in trackdiff] )
    hdr.append("static int trackdiff[] = {%s};" %str_diff )    
    str_diff_constant=', '.join([str(net.dict_types['Species'][nn].diffusion) for nn in trackdiff])
    hdr.append("static double diff_constant[] = {%s};" %str_diff_constant )#table containing diffusion constants of ligands
    list_ext=[]
    if 'Ligand' in net.dict_types:
        for nn in net.dict_types['Ligand']:
            if nn.isinstance('Diffusible'):
                list_ext.append(str(1))
            else:
                list_ext.append(str(0))

    str_ext=', '.join(list_ext)
    hdr.append("static int externallig[] = {%s};" %str_ext )#table containing external tags for ligands
    hdr.append("\n\n")
    
    return '\n'.join(hdr)    # note added the \n here between all elements of hdr

def write_deriv_inC(net,prmt,programm_file):
    """Write the integration equations in the C-file
    
    This function is a default and should be updated in Networks/interaction.py
    
    Args:
        net (Network): the network under study
        programm_file (TextIOWrapper): the built_integrator file
    Return:
        None: directly write the string in the C-file
    """
    add = programm_file.write #create a bound method for readibility
    add("#define ID_INTEGRATOR\n")
    add('void declare_list(){\n')
    for j in range(len(prmt['Ligands'])):
        add("LIGAND_LIST[%i] = %f;\n" %(j, prmt['Ligands'][j]) )
    for j in range(len(prmt['Ligand_self'])):
        add("LIGAND_SELF_LIST[%i] = %f;\n" %(j, prmt['Ligand_self'][j] ))
    for j in range(len(prmt['tau_off'])):
        add("TAU_LIST[%i] = %f;\n" %(j, prmt['tau_off'][j]) )
    list_pMHC=net.dict_types['pMHC']
    for j in range(len(list_pMHC)):
        add("pMHC_LIST[%i] = %i;\n" %(j, list_pMHC[j].int_id()) )
    
    S_list = net.dict_types['Species']
    counter = 0
    for i in range(len(S_list)):
        if (S_list[i].isinstance('Kinase') or S_list[i].isinstance('Phosphatase') ):
            if not S_list[i].isinstance('pMHC'):
                out = net.graph.successors(S_list[i])
                for j in range(len(out)):
                    if out[j].isinstance('Initial_Concentration'):
                        add("KP_LIST[%i] = %i;\n" %(counter, S_list[i].int_id()) )
                        add("KP_CONC[%i] = %f;\n" %(S_list[i].int_id(), out[j].conc) )
                        counter += 1
    add('}\n\n')
    add(compute_deriv_inC(net))

def write_program(programm_file,net, prmt, print_buf, Cseed=0):
    """Write the built_integrator of the network in the C file
    
    Collect python encoded C and the stored files selected via cfile
    dictionary and write them in the correct order.
    
    Args:
        programm_file (TextIOWrapper): the built_integrator file
        net (Network): -
        prmt (dict): passed to all_params2C
        print_buf (bool): passed to all_params2C
        Cseed (int): passed to all_params2C
    
    Return:
        str: the C programm as a python string
    """
    # these have to be loaded in this order due to implicit type def's
    required_files2 = ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main' ]
    programm_file.write(all_params2C(net, prmt, print_buf, Cseed))    
    programm_file.write(open(cfile['header']).read())
    programm_file.write(open(cfile['utilities']).read())
    programm_file.write('/***** end of header, begining of python computed functions ***/\n\n')
    write_deriv_inC(net,prmt,programm_file) #define in Networks/interaction.py
    programm_file.write('/***** end of python computed functions, beginning problem specific fns ***/\n\n')
    for file_name in required_files2:
        if file_name in cfile and cfile[file_name].endswith('.c'):  # omit files = ' ' etc
            programm_file.write(open(cfile[file_name]).read())

def modifier(deriv2):
    deriv2.all_params2C = all_params2C
    deriv2.write_deriv_inC = write_deriv_inC
    deriv2.write_program = write_program
    return deriv2
