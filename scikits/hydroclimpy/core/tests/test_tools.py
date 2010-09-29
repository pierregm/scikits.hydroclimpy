import numpy as np
import numpy.ma as ma
import scikits.timeseries as ts

from numpy.ma.testutils import *

from scikits.hydroclimpy import Bag, Cluster, flatten_sequence



class TestFunctions(TestCase):
        #
    def test_flatten_sequence(self):
        "Test flattening a sequence"
        test = flatten_sequence([0,[1,2,[3,4]],[5,6]])
        control = [0,1,2,3,4,5,6,]
        self.failUnless([_ for _ in test] == control)


class TestCluster(TestCase):
    #
    def test_cluster(self):
        sequence = [0, 0, 1, 2, 2, 2, 3, 4, 3, 4, 4, 4]
        klust = Cluster(sequence, 0)
        assert_equal(klust.clustered,
                     [[0,0,], [1,], [2,2,2,], [3,], [4,], [3,], [4,4,4,]])
        assert_equal(klust.uniques,
                     [0, 1, 2, 3, 4, 3, 4])
        assert_equal(len(klust.uniques), len(klust.clustered))
        assert_equal(klust.slices,
                     [slice(0,2), slice(2,3), slice(3,6), slice(6,7),
                      slice(7,8), slice(8,9), slice(9, 12)])
        #
        x = [ 1.8, 1.3, 2.4, 1.2, 2.5, 3.9, 1. , 3.8, 4.2, 3.3,
              1.2, 0.2, 0.9, 2.7, 2.4, 2.8, 2.7, 4.7, 4.2, 0.4]
        klust = Cluster(x, 1)
        assert_equal(klust.starts,
                     [ 0, 2, 3, 4, 5, 6, 7, 10, 13, 17, 19])
        assert_equal(Cluster(x, 1.5).starts,
                     [ 0, 6, 7, 10, 13, 17, 19])
        assert_equal(Cluster(x, 2.5).starts,
                     [ 0, 6, 7, 19])
        assert_equal(Cluster(x,2.5, np.greater).starts,
                     [ 0, 1, 2, 3, 4, 5, 8, 9,10,11,12,13,14,15,16,17,18])
    #
    def test_mark_greaterthan(self):
        y = [ 0, -1, 0, 0, 0, 1, 1, -1, -1, -1, 
              1,  1, 0, 0, 0, 0, 1,  1,  0,  0]
        klust = Cluster(y,0)
        assert_equal(klust.starts, [ 0,  1,  2,  5,  7, 10, 12, 16, 18])
        assert_equal(klust.mark_greaterthan(3),
                     [1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,1,1,1])
    #
    def test_grouped_slices(self):
        "Make sure we're not massing w/ .uniques ..."
        y = [ 0, -1, 0, 0, 0, 1, 1, -1, -1, -1, 
              1,  1, 0, 0, 0, 0, 1,  1,  0,  0]
        klust = Cluster(y,0)
        uniques_before = klust.uniques
        klust.grouped_slices()
        uniques_after = klust.uniques
        assert_equal(uniques_before, uniques_after)



class TestBag(TestCase):
    #
    def setUp(self):
        self.bag = Bag(a=1, b=2, c=3)
    #
    def test_getitem(self):
        "Try getitem/getattr"
        bag = self.bag
        self.failUnless(bag['a'] == 1)
        self.failUnless(bag.a == 1)
    #
    def test_creation(self):
        "Test the creation of a Bag"
        # With a dictionary
        bag = Bag(self.bag)
        self.failUnless(bag==self.bag)
        # With optional parameters
        bag = Bag(**self.bag)
        self.failUnless(bag==self.bag)
        # With a list of tuples
        bag = Bag([(k,v) for (k,v) in self.bag.items()])
        self.failUnless(bag==self.bag)
    #
    def test_set_by_attribute(self):
        "Try setting by attribute"
        bag = self.bag
        bag.d = 4
        self.failUnless(bag['d'] == 4)
        self.failUnless(bag.d == 4)
        #
        bag.a = -1
        self.failUnless(bag['a'] == 1)
    #
    def test_set_specials(self):
        "Try setting special functions"
        bag = self.bag
        try:
            bag.__dict__ = 7
        except TypeError:
            pass
        #
        bag.pop = 7
        self.failUnless(hasattr(bag.pop, '__call__'))

    def test_pickling(self):
        bag = self.bag
        import cPickle
        pickled = cPickle.loads(cPickle.dumps(bag))
        self.failUnless(pickled==bag)




if __name__ == '__main__':
    run_module_suite()
