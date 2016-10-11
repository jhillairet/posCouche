# -*- coding: utf-8 -*-
"""
Testing routines for formeTS

@author: J.Hillairet
"""
import pytest

from posCouche.formeTS import formeTS

def test_negative_shot_number():
    with pytest.raises(ValueError):
        formeTS(-1)

def test_non_integer_shot_number():
    with pytest.raises(ValueError):
        formeTS(10.2)
        


    