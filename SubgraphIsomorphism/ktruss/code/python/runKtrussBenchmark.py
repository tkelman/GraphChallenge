import numpy as np
import scipy as sp
from ktruss import ktruss

#Use the pandas package if available
#import pandas as pd

inc_mtx_file = '../../../data/amazon0302_inc.tsv';

E_expected=np.array([ (1, 1, 0 ,0 ,0), (0, 1, 1 ,0 ,0), (1, 0,0,1,0),(0,0,1,1,0),(1,0,1,0,0),(0,0,0,0,0) ])

E=ktruss(inc_mtx_file,3)

#if sum(sum(E.toarray()-E_expected)):
#    print "unable to verify"
#else:
#    print E

###################################################
# Graph Challenge benchmark
# Developer: Dr. Vijay Gadepally (vijayg@mit.edu)
# MIT
###################################################
# (c) <2015> Vijay Gadepally
###################################################
