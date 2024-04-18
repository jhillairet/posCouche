# -*- coding: utf-8 -*-
"""
Tests for ic_utils

@author: J.Hillairet
"""
import unittest

from numpy.testing import assert_almost_equal
from posCouche.ic_utils import IC_resonance_radius

class TestIC_resonance_radius(unittest.TestCase):
    def test_IC_resonance_radius_default_values(self):
        assert_almost_equal(2.529, IC_resonance_radius(Itor=1250, f=55, n=1, A=1), decimal=3)
    
    def test_IC_resonance_radius_default_values(self):
        assert_almost_equal(2.529, IC_resonance_radius(Itor=1250, f=55, n=1, A=1), decimal=3)   

    