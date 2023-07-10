import GPy
from sklearn.metrics import DistanceMetric
from numpy import fliplr, pi
from GPy.kern import Kern
from GPy.kern import RBF
from GPy.kern import Exponential

class RBFhaversine(RBF):
    # NOT VALID
    #
    # do I need to do this for every type of covariance function
    # or is there a general way of just specifying a new distance metric
    #
    #def __init__(self, input_dim, variance=1., lengthscale=None, ARD=False, active_dims=None, name='Exponential', *args, **kwargs):
#        print("Warning - ARD must be False at the moment")
        #assert ARD=FALSE
#        super(Exponentialhaversine, self).__init__(input_dim, variance, lengthscale, ARD, active_dims, name, *args, **kwargs)

    def _unscaled_dist(self, X, X2=None):
        """
        Compute the haversine distance between each row of X and X2, or between
        each pair of rows of X if X2 is None. First column must be longitude and
        the second latitude, both in degrees.
        """
        haversine = DistanceMetric.get_metric('haversine')
        if X2 is None:
            return 6371.*haversine.pairwise(fliplr(X)*pi/180.)
            # note sklearn haversine distance requires (lat, long) whereas we are
            # working with (long, lat). numpy.fliplr switches the columns of X.
        else:
            return 6371.*haversine.pairwise(fliplr(X)*pi/180., fliplr(X2)*pi/180.)

class Exponentialhaversine(Exponential):
    def __init__(self, input_dim, variance=1., lengthscale=None, ARD=False, active_dims=None, name='Exponential', *args, **kwargs):
        print("Warning - ARD must be False at the moment")
        super(Exponentialhaversine, self).__init__(input_dim, variance, lengthscale, ARD, active_dims, name, *args, **kwargs)

    def _unscaled_dist(self, X, X2=None):
        """
        Compute the haversine distance between each row of X and X2, or between
        each pair of rows of X if X2 is None. First column must be longitude and
        the second latitude, both in degrees.
        """
        #print("WARNING - ARD=False only")
        haversine = DistanceMetric.get_metric('haversine')
        if X2 is None:
            return 6371.*haversine.pairwise(fliplr(X)*pi/180.)
            # note sklearn haversine distance requires (lat, long) whereas we are
            # working with (long, lat). numpy.fliplr switches the columns of X.
        else:
            return 6371.*haversine.pairwise(fliplr(X)*pi/180., fliplr(X2)*pi/180.)
