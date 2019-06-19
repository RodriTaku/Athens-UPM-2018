# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 17:34:07 2018

@author: RodrigoTakuro
"""

class paramStruct:
    """paramStruct Class for cartPend
    
    Parameters
    ----------
    cartMass : int, float
        cart mass M (default=1./2)
    bobMass : int,float
        bob mass m (default=1./4)
    poleLength : int, float
        pole length L (default=1.)
    gravity : int, float
        gravitational constant (default=3./4)
    """
    def __init__(self, cartMass=1./2, bobMass=1./4, poleLength=1., gravity=3./4):
        self.M = float(cartMass)
        self.m = float(bobMass)
        self.L = float(poleLength)
        self.g = float(gravity)