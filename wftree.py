""" estimate tree size of the whole population
when the population is haploid, of fixed size N,
and evolving according to the Wright-Fisher model
"""

import sys
import numpy as np
from scipy.stats import geom


def estimate_tree_size(pop__size):
    """ denote population size (pop__size) with N """
    ## """ initialise number of blocks
    __u = pop__size
    ## initialise tree size """
    tree_size= float(0)
    while __u > 1:
        tree_size += float(__u)
        ## sample __u parents from the population of size pop_sze
        ## and return  the number of unique parents
        if __u > 2:
            __u = np.unique( np.random.choice( pop__size,  __u ) ).size
        else :
        ## only two blocks left, so sample from a Geometric
        ## with probability  1/N of success
            __u = 1
            tree_size += float(geom.rvs( 1.0/float(pop__size), size=1))
    print( tree_size/(2.0* float(pop__size)*np.log(float( pop__size))))

if __name__ == '__main__':
   for r in np.arange(1000):
        estimate_tree_size(int(sys.argv[1]))
