import os

_default_text = "This is a STOP file.\nDeleting this file will stop the running evolution."

def test_STOP_file(file_name,current_state=None):
    """
        Test whether the stop file generated at the beginning of each evolution is still present.
        If not stop the running evolution.


        Args:
            file_name: path of the STOP file
    """
    if not os.path.isfile(file_name):
        __stop_evolution__()
    elif current_state:
        edit_state_STOP_file(file_name,current_state)

def edit_state_STOP_file(file_name,current_state):
    """
    Write the current state of the evolution in the Stop file.

    Args:
        file_name: path of the STOP file
        current_state: dictionnary of information regarding the state of the evolution.
                        The dictionnary needs at least the keys "seed", and "generation","fitness".
    """
    data = _default_text
    data+="\n\n --- Current State ---"
    seed = current_state.pop("seed")
    generation = current_state.pop("generation")
    fitness = current_state.pop("fitness")
    data+= "\nSeed:{0}\ngeneration = {1} (fitness: {2})".format(seed,generation,fitness)
    for key,val in current_state.items():
        data+=("\n{0} : {1}".format(key,val))
    data+="\n"
    with open(file_name,"w") as stop_file:
        stop_file.write(data)

def create_STOP_file(model_dir):
    """
        Write the stop file used to interrupt a simulation

        Args:
            model_dir: path of the project directory. The stop file is written there.
        Returns:
            the path of the stop file
    """
    path = os.path.abspath(os.path.join(model_dir,"STOP.txt"))
    with open(path,"w") as stop_file:
        stop_file.write(_default_text)
    return path

def __stop_evolution__():
    print("\nSTOP file deleted, this stops the evolution.")
    os._exit(1)
