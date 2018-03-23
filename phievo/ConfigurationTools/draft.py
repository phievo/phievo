


if PRINT:
    print("Kinetic parameters:")
    for key,val in mutation.dictionary_ranges.items():
        print("\t",key,":",val)

    print("\nMutation parameters:")
    for key,val in initialization.dictionary_mutation.items():
        print("\t",key,":",val)    

    print("\npfile")
    for key,val in pfile.items():
        print("\t",key,":",val)

    print("\ncfile")
    for key,val in cfile.items():
        print("\t",key,":",val)

    print("\nprmt")
    for key,val in prmt.items():
        print("\t",key,":",val)

    print("\nrestart")
    for key,val in restart.items():
        print("\t",key,":",val)
