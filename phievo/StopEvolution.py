import os

def test_STOP_file(file_name):
    """
        Test whether the stop file generated at the beginning of each evolution is still present.
        If not stop the running evolution.
        

        Args:
            file_name: path of the STOP file
    """
    if not os.path.isfile(file_name):
        __stop_evolution__()

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
        stop_file.write("This is a STOP file.\n")
        stop_file.write("Deleting this file will stop the running evolution.")
    return path

def __stop_evolution__():
    os._exit(1)
