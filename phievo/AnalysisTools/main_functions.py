import numpy as np


def smoothing(array,param):
    """Smoothen an array by averaging over the neighbourhood
    
    Args:
        array (list): the to be smoothed array
        param (int): the distance of the neighbourhood
    
    Returns:
        list: of same size as array
    """
    length = len(array)
    return [np.mean(array[max(0,i-param):i+param+1]) for i in range(length)]
