import numpy as np

def intersection(lst1, lst2):
    """Intersection of two lists

    Args:
        lst1 (_type_): _description_
        lst2 (_type_): _description_

    Returns:
        _type_: _description_
    """    
    return list(set(lst1) & set(lst2))

def confidence_interval(a):
    """_summary_

    Args:
        a (_type_): _description_

    Returns:
        _type_: _description_
    """    
    sorted = np.sort(a)
    boundary = int(np.around(0.026*len(a)))
    lower = sorted[boundary]
    higher = sorted[-boundary]
    return lower, higher