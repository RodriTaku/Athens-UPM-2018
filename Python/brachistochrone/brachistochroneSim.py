# -*- coding: utf-8 -*-
"""
Created on March, 2018

@author: RodrigoTakuro
"""
#from math import pi
from numpy import array, ones, zeros, append, sqrt, sin, cos, linspace, nan
from scipy.optimize import root
from scipy.sparse import csr_matrix

def brachobj(Z, x, y, v, g, Jrow, Jcol):
    """Brachistochrone ojective or indicator function
    
    Parameters
    ----------
    Z : list or ndarray
        internal state
    x : list or ndarray
        horizontal positions
    y : list or ndarray
        vertical positions (positive is downwards)
    v : list or ndarray
        norms of the velocity vector
    g : float
        gravity
    Jrow : list or numpy.array
        row indices of nonzero entries of the jacobian
    Jcol : list or numpy.array
        column indices of nonzero entries of the jacobian
        
    Returns
    -------
    J : float
        Integrator equations to be solved
    """
    return 10*Z[0] - Z[15]*(Z[4] - Z[11] + Z[0]*g*cos(Z[5])) - Z[22]*(Z[11] - Z[18] + Z[0]*g*cos(Z[12])) - Z[29]*(Z[18] - Z[25] + Z[0]*g*cos(Z[19])) - Z[36]*(Z[25] - Z[32] + Z[0]*g*cos(Z[26])) - Z[43]*(Z[32] - Z[39] + Z[0]*g*cos(Z[33])) - Z[50]*(Z[39] - Z[46] + Z[0]*g*cos(Z[40])) - Z[57]*(Z[46] - Z[53] + Z[0]*g*cos(Z[47])) - Z[64]*(Z[53] - Z[60] + Z[0]*g*cos(Z[54])) - Z[68]*(Z[60] - Z[65] + Z[0]*g*cos(Z[61])) - Z[13]*(Z[2] - Z[9] + Z[0]*Z[4]*sin(Z[5])) - Z[20]*(Z[9] - Z[16] + Z[0]*Z[11]*sin(Z[12])) - Z[27]*(Z[16] - Z[23] + Z[0]*Z[18]*sin(Z[19])) - Z[34]*(Z[23] - Z[30] + Z[0]*Z[25]*sin(Z[26])) - Z[41]*(Z[30] - Z[37] + Z[0]*Z[32]*sin(Z[33])) - Z[48]*(Z[37] - Z[44] + Z[0]*Z[39]*sin(Z[40])) - Z[55]*(Z[44] - Z[51] + Z[0]*Z[46]*sin(Z[47])) - Z[62]*(Z[51] - Z[58] + Z[0]*Z[53]*sin(Z[54])) - Z[67]*(Z[59] - y[10] + Z[0]*Z[60]*cos(Z[61])) - Z[8]*(v[0] - Z[4] + Z[0]*g*cos(Z[1])) - Z[66]*(Z[58] - x[10] + Z[0]*Z[60]*sin(Z[61])) - Z[7]*(y[0] - Z[3] + Z[0]*v[0]*cos(Z[1])) - Z[6]*(x[0] - Z[2] + Z[0]*v[0]*sin(Z[1])) - Z[14]*(Z[3] - Z[10] + Z[0]*Z[4]*cos(Z[5])) - Z[21]*(Z[10] - Z[17] + Z[0]*Z[11]*cos(Z[12])) - Z[28]*(Z[17] - Z[24] + Z[0]*Z[18]*cos(Z[19])) - Z[35]*(Z[24] - Z[31] + Z[0]*Z[25]*cos(Z[26])) - Z[42]*(Z[31] - Z[38] + Z[0]*Z[32]*cos(Z[33])) - Z[49]*(Z[38] - Z[45] + Z[0]*Z[39]*cos(Z[40])) - Z[56]*(Z[45] - Z[52] + Z[0]*Z[46]*cos(Z[47])) - Z[63]*(Z[52] - Z[59] + Z[0]*Z[53]*cos(Z[54]))

def brachfun(Z, x, y, v, g, Jrow, Jcol):
    """Brachistochrone control equations
    
    Parameters
    ----------
    Z : list or ndarray
        internal state
    x : list or ndarray
        horizontal positions
    y : list or ndarray
        vertical positions (positive is downwards)
    v : list or ndarray
        norms of the velocity vector
    g : float
        gravity
    Jrow : list or numpy.array
        row indices of nonzero entries of the jacobian
    Jcol : list or numpy.array
        column indices of nonzero entries of the jacobian
        
    Returns
    -------
    Fun : numpy.array
        Integrator equations to be solved
    """
    print "Objective function value: {}".format(brachobj(Z, x, y, v, g, Jrow, Jcol))
    return array([10 - Z[11]*Z[21]*cos(Z[12]) - Z[18]*Z[28]*cos(Z[19]) - Z[25]*Z[35]*cos(Z[26]) - Z[32]*Z[42]*cos(Z[33]) - Z[39]*Z[49]*cos(Z[40]) - Z[46]*Z[56]*cos(Z[47]) - Z[53]*Z[63]*cos(Z[54]) - Z[60]*Z[67]*cos(Z[61]) - Z[8]*g*cos(Z[1]) - Z[15]*g*cos(Z[5]) - Z[22]*g*cos(Z[12]) - Z[29]*g*cos(Z[19]) - Z[36]*g*cos(Z[26]) - Z[43]*g*cos(Z[33]) - Z[50]*g*cos(Z[40]) - Z[57]*g*cos(Z[47]) - Z[64]*g*cos(Z[54]) - Z[68]*g*cos(Z[61]) - Z[4]*Z[13]*sin(Z[5]) - Z[11]*Z[20]*sin(Z[12]) - Z[18]*Z[27]*sin(Z[19]) - Z[25]*Z[34]*sin(Z[26]) - Z[32]*Z[41]*sin(Z[33]) - Z[39]*Z[48]*sin(Z[40]) - Z[46]*Z[55]*sin(Z[47]) - Z[53]*Z[62]*sin(Z[54]) - Z[60]*Z[66]*sin(Z[61]) - Z[7]*v[0]*cos(Z[1]) - Z[6]*v[0]*sin(Z[1]) - Z[4]*Z[14]*cos(Z[5]),
                Z[0]*Z[8]*g*sin(Z[1]) - Z[0]*Z[6]*v[0]*cos(Z[1]) + Z[0]*Z[7]*v[0]*sin(Z[1]),
                Z[2] - x[0] - Z[0]*v[0]*sin(Z[1]),
                Z[3] - y[0] - Z[0]*v[0]*cos(Z[1]),
                Z[4] - v[0] - Z[0]*g*cos(Z[1]),
                Z[0]*Z[4]*Z[14]*sin(Z[5]) - Z[0]*Z[4]*Z[13]*cos(Z[5]) + Z[0]*Z[15]*g*sin(Z[5]),
                Z[6] - Z[13],
                Z[7] - Z[14],
                Z[8] - Z[15] - Z[0]*Z[14]*cos(Z[5]) - Z[0]*Z[13]*sin(Z[5]),
                Z[9] - Z[2] - Z[0]*Z[4]*sin(Z[5]),
                Z[10] - Z[3] - Z[0]*Z[4]*cos(Z[5]),
                Z[11] - Z[4] - Z[0]*g*cos(Z[5]),
                Z[0]*Z[11]*Z[21]*sin(Z[12]) - Z[0]*Z[11]*Z[20]*cos(Z[12]) + Z[0]*Z[22]*g*sin(Z[12]),
                Z[13] - Z[20],
                Z[14] - Z[21],
                Z[15] - Z[22] - Z[0]*Z[21]*cos(Z[12]) - Z[0]*Z[20]*sin(Z[12]),
                Z[16] - Z[9] - Z[0]*Z[11]*sin(Z[12]),
                Z[17] - Z[10] - Z[0]*Z[11]*cos(Z[12]),
                Z[18] - Z[11] - Z[0]*g*cos(Z[12]),
                Z[0]*Z[18]*Z[28]*sin(Z[19]) - Z[0]*Z[18]*Z[27]*cos(Z[19]) + Z[0]*Z[29]*g*sin(Z[19]),
                Z[20] - Z[27],
                Z[21] - Z[28],
                Z[22] - Z[29] - Z[0]*Z[28]*cos(Z[19]) - Z[0]*Z[27]*sin(Z[19]),
                Z[23] - Z[16] - Z[0]*Z[18]*sin(Z[19]),
                Z[24] - Z[17] - Z[0]*Z[18]*cos(Z[19]),
                Z[25] - Z[18] - Z[0]*g*cos(Z[19]),
                Z[0]*Z[25]*Z[35]*sin(Z[26]) - Z[0]*Z[25]*Z[34]*cos(Z[26]) + Z[0]*Z[36]*g*sin(Z[26]),
                Z[27] - Z[34],
                Z[28] - Z[35],
                Z[29] - Z[36] - Z[0]*Z[35]*cos(Z[26]) - Z[0]*Z[34]*sin(Z[26]),
                Z[30] - Z[23] - Z[0]*Z[25]*sin(Z[26]),
                Z[31] - Z[24] - Z[0]*Z[25]*cos(Z[26]),
                Z[32] - Z[25] - Z[0]*g*cos(Z[26]),
                Z[0]*Z[32]*Z[42]*sin(Z[33]) - Z[0]*Z[32]*Z[41]*cos(Z[33]) + Z[0]*Z[43]*g*sin(Z[33]),
                Z[34] - Z[41],
                Z[35] - Z[42],
                Z[36] - Z[43] - Z[0]*Z[42]*cos(Z[33]) - Z[0]*Z[41]*sin(Z[33]),
                Z[37] - Z[30] - Z[0]*Z[32]*sin(Z[33]),
                Z[38] - Z[31] - Z[0]*Z[32]*cos(Z[33]),
                Z[39] - Z[32] - Z[0]*g*cos(Z[33]),
                Z[0]*Z[39]*Z[49]*sin(Z[40]) - Z[0]*Z[39]*Z[48]*cos(Z[40]) + Z[0]*Z[50]*g*sin(Z[40]),
                Z[41] - Z[48],
                Z[42] - Z[49],
                Z[43] - Z[50] - Z[0]*Z[49]*cos(Z[40]) - Z[0]*Z[48]*sin(Z[40]),
                Z[44] - Z[37] - Z[0]*Z[39]*sin(Z[40]),
                Z[45] - Z[38] - Z[0]*Z[39]*cos(Z[40]),
                Z[46] - Z[39] - Z[0]*g*cos(Z[40]),
                Z[0]*Z[46]*Z[56]*sin(Z[47]) - Z[0]*Z[46]*Z[55]*cos(Z[47]) + Z[0]*Z[57]*g*sin(Z[47]),
                Z[48] - Z[55],
                Z[49] - Z[56],
                Z[50] - Z[57] - Z[0]*Z[56]*cos(Z[47]) - Z[0]*Z[55]*sin(Z[47]),
                Z[51] - Z[44] - Z[0]*Z[46]*sin(Z[47]),
                Z[52] - Z[45] - Z[0]*Z[46]*cos(Z[47]),
                Z[53] - Z[46] - Z[0]*g*cos(Z[47]),
                Z[0]*Z[53]*Z[63]*sin(Z[54]) - Z[0]*Z[53]*Z[62]*cos(Z[54]) + Z[0]*Z[64]*g*sin(Z[54]),
                Z[55] - Z[62],
                Z[56] - Z[63],
                Z[57] - Z[64] - Z[0]*Z[63]*cos(Z[54]) - Z[0]*Z[62]*sin(Z[54]),
                Z[58] - Z[51] - Z[0]*Z[53]*sin(Z[54]),
                Z[59] - Z[52] - Z[0]*Z[53]*cos(Z[54]),
                Z[60] - Z[53] - Z[0]*g*cos(Z[54]),
                Z[0]*Z[60]*Z[67]*sin(Z[61]) - Z[0]*Z[60]*Z[66]*cos(Z[61]) + Z[0]*Z[68]*g*sin(Z[61]),
                Z[62] - Z[66],
                Z[63] - Z[67],
                Z[64] - Z[68] - Z[0]*Z[67]*cos(Z[61]) - Z[0]*Z[66]*sin(Z[61]),
                x[10] - Z[58] - Z[0]*Z[60]*sin(Z[61]),
                y[10] - Z[59] - Z[0]*Z[60]*cos(Z[61]),
                Z[65] - Z[60] - Z[0]*g*cos(Z[61]),
                Z[68]])

    
def brachdfun(Z, x, y, v, g, Jrow, Jcol):
    """Brachistochrone control equations jacobian
    
    Parameters
    ----------
    Z : list or ndarray
        internal state
    x : list or ndarray
        horizontal positions
    y : list or ndarray
        vertical positions (positive is downwards)
    v : list or ndarray
        norms of the velocity vector
    g : float
        gravity
    Jrow : list or numpy.array
        row indices of nonzero entries of the jacobian
    Jcol : list or numpy.array
        column indices of nonzero entries of the jacobian
        
    Returns
    -------
    Fun : numpy.array
        Jacobian of the control equations
    """
    Jvec = array([ Z[8]*g*sin(Z[1]) - Z[6]*v[0]*cos(Z[1]) + Z[7]*v[0]*sin(Z[1]),
                -v[0]*sin(Z[1]),
                -v[0]*cos(Z[1]),
                -g*cos(Z[1]),
                Z[4]*Z[14]*sin(Z[5]) - Z[4]*Z[13]*cos(Z[5]) + Z[15]*g*sin(Z[5]),
                - Z[14]*cos(Z[5]) - Z[13]*sin(Z[5]),
                -Z[4]*sin(Z[5]),
                -Z[4]*cos(Z[5]),
                -g*cos(Z[5]),
                Z[11]*Z[21]*sin(Z[12]) - Z[11]*Z[20]*cos(Z[12]) + Z[22]*g*sin(Z[12]),
                - Z[21]*cos(Z[12]) - Z[20]*sin(Z[12]),
                -Z[11]*sin(Z[12]),
                -Z[11]*cos(Z[12]),
                -g*cos(Z[12]),
                Z[18]*Z[28]*sin(Z[19]) - Z[18]*Z[27]*cos(Z[19]) + Z[29]*g*sin(Z[19]),
                - Z[28]*cos(Z[19]) - Z[27]*sin(Z[19]),
                -Z[18]*sin(Z[19]),
                -Z[18]*cos(Z[19]),
                -g*cos(Z[19]),
                Z[25]*Z[35]*sin(Z[26]) - Z[25]*Z[34]*cos(Z[26]) + Z[36]*g*sin(Z[26]),
                - Z[35]*cos(Z[26]) - Z[34]*sin(Z[26]),
                -Z[25]*sin(Z[26]),
                -Z[25]*cos(Z[26]),
                -g*cos(Z[26]),
                Z[32]*Z[42]*sin(Z[33]) - Z[32]*Z[41]*cos(Z[33]) + Z[43]*g*sin(Z[33]),
                - Z[42]*cos(Z[33]) - Z[41]*sin(Z[33]),
                -Z[32]*sin(Z[33]),
                -Z[32]*cos(Z[33]),
                -g*cos(Z[33]),
                Z[39]*Z[49]*sin(Z[40]) - Z[39]*Z[48]*cos(Z[40]) + Z[50]*g*sin(Z[40]),
                - Z[49]*cos(Z[40]) - Z[48]*sin(Z[40]),
                -Z[39]*sin(Z[40]),
                -Z[39]*cos(Z[40]),
                -g*cos(Z[40]),
                Z[46]*Z[56]*sin(Z[47]) - Z[46]*Z[55]*cos(Z[47]) + Z[57]*g*sin(Z[47]),
                - Z[56]*cos(Z[47]) - Z[55]*sin(Z[47]),
                -Z[46]*sin(Z[47]),
                -Z[46]*cos(Z[47]),
                -g*cos(Z[47]),
                Z[53]*Z[63]*sin(Z[54]) - Z[53]*Z[62]*cos(Z[54]) + Z[64]*g*sin(Z[54]),
                - Z[63]*cos(Z[54]) - Z[62]*sin(Z[54]),
                -Z[53]*sin(Z[54]),
                -Z[53]*cos(Z[54]),
                -g*cos(Z[54]),
                Z[60]*Z[67]*sin(Z[61]) - Z[60]*Z[66]*cos(Z[61]) + Z[68]*g*sin(Z[61]),
                - Z[67]*cos(Z[61]) - Z[66]*sin(Z[61]),
                -Z[60]*sin(Z[61]),
                -Z[60]*cos(Z[61]),
                -g*cos(Z[61]),
                Z[8]*g*sin(Z[1]) - Z[6]*v[0]*cos(Z[1]) + Z[7]*v[0]*sin(Z[1]),
                Z[0]*Z[8]*g*cos(Z[1]) + Z[0]*Z[7]*v[0]*cos(Z[1]) + Z[0]*Z[6]*v[0]*sin(Z[1]),
                -Z[0]*v[0]*cos(Z[1]),
                Z[0]*v[0]*sin(Z[1]),
                Z[0]*g*sin(Z[1]),
                1,
                -1,
                1,
                -1,
                - Z[14]*cos(Z[5]) - Z[13]*sin(Z[5]),
                1,
                Z[0]*Z[14]*sin(Z[5]) - Z[0]*Z[13]*cos(Z[5]),
                -Z[0]*sin(Z[5]),
                -Z[0]*cos(Z[5]),
                -1,
                Z[4]*Z[14]*sin(Z[5]) - Z[4]*Z[13]*cos(Z[5]) + Z[15]*g*sin(Z[5]),
                Z[0]*Z[4]*Z[14]*cos(Z[5]) + Z[0]*Z[15]*g*cos(Z[5]) + Z[0]*Z[4]*Z[13]*sin(Z[5]),
                Z[0]*Z[14]*sin(Z[5]) - Z[0]*Z[13]*cos(Z[5]),
                -Z[0]*Z[4]*cos(Z[5]),
                Z[0]*Z[4]*sin(Z[5]),
                Z[0]*g*sin(Z[5]),
                -v[0]*sin(Z[1]),
                -Z[0]*v[0]*cos(Z[1]),
                1,
                -v[0]*cos(Z[1]),
                Z[0]*v[0]*sin(Z[1]),
                1,
                -g*cos(Z[1]),
                Z[0]*g*sin(Z[1]),
                1,
                1,
                -1,
                1,
                -1,
                - Z[21]*cos(Z[12]) - Z[20]*sin(Z[12]),
                1,
                Z[0]*Z[21]*sin(Z[12]) - Z[0]*Z[20]*cos(Z[12]),
                -Z[0]*sin(Z[12]),
                -Z[0]*cos(Z[12]),
                -1,
                Z[11]*Z[21]*sin(Z[12]) - Z[11]*Z[20]*cos(Z[12]) + Z[22]*g*sin(Z[12]),
                Z[0]*Z[11]*Z[21]*cos(Z[12]) + Z[0]*Z[22]*g*cos(Z[12]) + Z[0]*Z[11]*Z[20]*sin(Z[12]),
                Z[0]*Z[21]*sin(Z[12]) - Z[0]*Z[20]*cos(Z[12]),
                -Z[0]*Z[11]*cos(Z[12]),
                Z[0]*Z[11]*sin(Z[12]),
                Z[0]*g*sin(Z[12]),
                -Z[4]*sin(Z[5]),
                -Z[0]*Z[4]*cos(Z[5]),
                -1,
                -Z[0]*sin(Z[5]),
                1,
                -Z[4]*cos(Z[5]),
                Z[0]*Z[4]*sin(Z[5]),
                -1,
                -Z[0]*cos(Z[5]),
                1,
                -g*cos(Z[5]),
                Z[0]*g*sin(Z[5]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[28]*cos(Z[19]) - Z[27]*sin(Z[19]),
                1,
                Z[0]*Z[28]*sin(Z[19]) - Z[0]*Z[27]*cos(Z[19]),
                -Z[0]*sin(Z[19]),
                -Z[0]*cos(Z[19]),
                -1,
                Z[18]*Z[28]*sin(Z[19]) - Z[18]*Z[27]*cos(Z[19]) + Z[29]*g*sin(Z[19]),
                Z[0]*Z[18]*Z[28]*cos(Z[19]) + Z[0]*Z[29]*g*cos(Z[19]) + Z[0]*Z[18]*Z[27]*sin(Z[19]),
                Z[0]*Z[28]*sin(Z[19]) - Z[0]*Z[27]*cos(Z[19]),
                -Z[0]*Z[18]*cos(Z[19]),
                Z[0]*Z[18]*sin(Z[19]),
                Z[0]*g*sin(Z[19]),
                -Z[11]*sin(Z[12]),
                -Z[0]*Z[11]*cos(Z[12]),
                -1,
                -Z[0]*sin(Z[12]),
                1,
                -Z[11]*cos(Z[12]),
                Z[0]*Z[11]*sin(Z[12]),
                -1,
                -Z[0]*cos(Z[12]),
                1,
                -g*cos(Z[12]),
                Z[0]*g*sin(Z[12]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[35]*cos(Z[26]) - Z[34]*sin(Z[26]),
                1,
                Z[0]*Z[35]*sin(Z[26]) - Z[0]*Z[34]*cos(Z[26]),
                -Z[0]*sin(Z[26]),
                -Z[0]*cos(Z[26]),
                -1,
                Z[25]*Z[35]*sin(Z[26]) - Z[25]*Z[34]*cos(Z[26]) + Z[36]*g*sin(Z[26]),
                Z[0]*Z[25]*Z[35]*cos(Z[26]) + Z[0]*Z[36]*g*cos(Z[26]) + Z[0]*Z[25]*Z[34]*sin(Z[26]),
                Z[0]*Z[35]*sin(Z[26]) - Z[0]*Z[34]*cos(Z[26]),
                -Z[0]*Z[25]*cos(Z[26]),
                Z[0]*Z[25]*sin(Z[26]),
                Z[0]*g*sin(Z[26]),
                -Z[18]*sin(Z[19]),
                -Z[0]*Z[18]*cos(Z[19]),
                -1,
                -Z[0]*sin(Z[19]),
                1,
                -Z[18]*cos(Z[19]),
                Z[0]*Z[18]*sin(Z[19]),
                -1,
                -Z[0]*cos(Z[19]),
                1,
                -g*cos(Z[19]),
                Z[0]*g*sin(Z[19]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[42]*cos(Z[33]) - Z[41]*sin(Z[33]),
                1,
                Z[0]*Z[42]*sin(Z[33]) - Z[0]*Z[41]*cos(Z[33]),
                -Z[0]*sin(Z[33]),
                -Z[0]*cos(Z[33]),
                -1,
                Z[32]*Z[42]*sin(Z[33]) - Z[32]*Z[41]*cos(Z[33]) + Z[43]*g*sin(Z[33]),
                Z[0]*Z[32]*Z[42]*cos(Z[33]) + Z[0]*Z[43]*g*cos(Z[33]) + Z[0]*Z[32]*Z[41]*sin(Z[33]),
                Z[0]*Z[42]*sin(Z[33]) - Z[0]*Z[41]*cos(Z[33]),
                -Z[0]*Z[32]*cos(Z[33]),
                Z[0]*Z[32]*sin(Z[33]),
                Z[0]*g*sin(Z[33]),
                -Z[25]*sin(Z[26]),
                -Z[0]*Z[25]*cos(Z[26]),
                -1,
                -Z[0]*sin(Z[26]),
                1,
                -Z[25]*cos(Z[26]),
                Z[0]*Z[25]*sin(Z[26]),
                -1,
                -Z[0]*cos(Z[26]),
                1,
                -g*cos(Z[26]),
                Z[0]*g*sin(Z[26]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[49]*cos(Z[40]) - Z[48]*sin(Z[40]),
                1,
                Z[0]*Z[49]*sin(Z[40]) - Z[0]*Z[48]*cos(Z[40]),
                -Z[0]*sin(Z[40]),
                -Z[0]*cos(Z[40]),
                -1,
                Z[39]*Z[49]*sin(Z[40]) - Z[39]*Z[48]*cos(Z[40]) + Z[50]*g*sin(Z[40]),
                Z[0]*Z[39]*Z[49]*cos(Z[40]) + Z[0]*Z[50]*g*cos(Z[40]) + Z[0]*Z[39]*Z[48]*sin(Z[40]),
                Z[0]*Z[49]*sin(Z[40]) - Z[0]*Z[48]*cos(Z[40]),
                -Z[0]*Z[39]*cos(Z[40]),
                Z[0]*Z[39]*sin(Z[40]),
                Z[0]*g*sin(Z[40]),
                -Z[32]*sin(Z[33]),
                -Z[0]*Z[32]*cos(Z[33]),
                -1,
                -Z[0]*sin(Z[33]),
                1,
                -Z[32]*cos(Z[33]),
                Z[0]*Z[32]*sin(Z[33]),
                -1,
                -Z[0]*cos(Z[33]),
                1,
                -g*cos(Z[33]),
                Z[0]*g*sin(Z[33]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[56]*cos(Z[47]) - Z[55]*sin(Z[47]),
                1,
                Z[0]*Z[56]*sin(Z[47]) - Z[0]*Z[55]*cos(Z[47]),
                -Z[0]*sin(Z[47]),
                -Z[0]*cos(Z[47]),
                -1,
                Z[46]*Z[56]*sin(Z[47]) - Z[46]*Z[55]*cos(Z[47]) + Z[57]*g*sin(Z[47]),
                Z[0]*Z[46]*Z[56]*cos(Z[47]) + Z[0]*Z[57]*g*cos(Z[47]) + Z[0]*Z[46]*Z[55]*sin(Z[47]),
                Z[0]*Z[56]*sin(Z[47]) - Z[0]*Z[55]*cos(Z[47]),
                -Z[0]*Z[46]*cos(Z[47]),
                Z[0]*Z[46]*sin(Z[47]),
                Z[0]*g*sin(Z[47]),
                -Z[39]*sin(Z[40]),
                -Z[0]*Z[39]*cos(Z[40]),
                -1,
                -Z[0]*sin(Z[40]),
                1,
                -Z[39]*cos(Z[40]),
                Z[0]*Z[39]*sin(Z[40]),
                -1,
                -Z[0]*cos(Z[40]),
                1,
                -g*cos(Z[40]),
                Z[0]*g*sin(Z[40]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[63]*cos(Z[54]) - Z[62]*sin(Z[54]),
                1,
                Z[0]*Z[63]*sin(Z[54]) - Z[0]*Z[62]*cos(Z[54]),
                -Z[0]*sin(Z[54]),
                -Z[0]*cos(Z[54]),
                -1,
                Z[53]*Z[63]*sin(Z[54]) - Z[53]*Z[62]*cos(Z[54]) + Z[64]*g*sin(Z[54]),
                Z[0]*Z[53]*Z[63]*cos(Z[54]) + Z[0]*Z[64]*g*cos(Z[54]) + Z[0]*Z[53]*Z[62]*sin(Z[54]),
                Z[0]*Z[63]*sin(Z[54]) - Z[0]*Z[62]*cos(Z[54]),
                -Z[0]*Z[53]*cos(Z[54]),
                Z[0]*Z[53]*sin(Z[54]),
                Z[0]*g*sin(Z[54]),
                -Z[46]*sin(Z[47]),
                -Z[0]*Z[46]*cos(Z[47]),
                -1,
                -Z[0]*sin(Z[47]),
                1,
                -Z[46]*cos(Z[47]),
                Z[0]*Z[46]*sin(Z[47]),
                -1,
                -Z[0]*cos(Z[47]),
                1,
                -g*cos(Z[47]),
                Z[0]*g*sin(Z[47]),
                -1,
                1,
                1,
                -1,
                1,
                -1,
                - Z[67]*cos(Z[61]) - Z[66]*sin(Z[61]),
                1,
                Z[0]*Z[67]*sin(Z[61]) - Z[0]*Z[66]*cos(Z[61]),
                -Z[0]*sin(Z[61]),
                -Z[0]*cos(Z[61]),
                -1,
                Z[60]*Z[67]*sin(Z[61]) - Z[60]*Z[66]*cos(Z[61]) + Z[68]*g*sin(Z[61]),
                Z[0]*Z[60]*Z[67]*cos(Z[61]) + Z[0]*Z[68]*g*cos(Z[61]) + Z[0]*Z[60]*Z[66]*sin(Z[61]),
                Z[0]*Z[67]*sin(Z[61]) - Z[0]*Z[66]*cos(Z[61]),
                -Z[0]*Z[60]*cos(Z[61]),
                Z[0]*Z[60]*sin(Z[61]),
                Z[0]*g*sin(Z[61]),
                -Z[53]*sin(Z[54]),
                -Z[0]*Z[53]*cos(Z[54]),
                -1,
                -Z[0]*sin(Z[54]),
                1,
                -Z[53]*cos(Z[54]),
                Z[0]*Z[53]*sin(Z[54]),
                -1,
                -Z[0]*cos(Z[54]),
                1,
                -g*cos(Z[54]),
                Z[0]*g*sin(Z[54]),
                -1,
                1,
                1,
                -Z[60]*sin(Z[61]),
                -Z[0]*Z[60]*cos(Z[61]),
                -1,
                -Z[0]*sin(Z[61]),
                -Z[60]*cos(Z[61]),
                Z[0]*Z[60]*sin(Z[61]),
                -1,
                -Z[0]*cos(Z[61]),
                -g*cos(Z[61]),
                Z[0]*g*sin(Z[61]),
                -1,
                1])
    return csr_matrix((Jvec, (Jrow, Jcol))).todense()


def brachsim( x0, y0, x1, y1, h0, g ):
    """Brachistochrone solver using a procedure based on the
    Hamilton-Pontryagin action principle.
    
    Parameters
    ----------
    x0 : int, float
        initial horizontal position
    y0 : int, float
        initial vertical position
    x1 : int, float
        final horizontal position
    y1 : int, float
        final vertical position
    h0 : float
        initial guess for the time step
    g : float
        gravity
    
    Returns
    -------
    Sol : ndarray
        array containing the complete evolution of the system.
        Each column corresponds to a time slice.
        Rows correspond to x, y, v, u, px, py, pv in that order
        
    Example
    -------    
    >>> my_sol = brachsim( 0., 0., 3., 2., 0.5, 1 )
    """
    #Load row and column indices to reconstruct sparse jacobian
    idx = [  0,  4,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 73, 75, 76, 77]
    Jrow = [  1,  2,  3,  4,  5,  8,  9, 10, 11, 12, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 29, 30, 31, 32, 33, 36, 37, 38, 39, 40, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 57, 58, 59, 60, 61, 64, 65, 66, 67,  0,  1,  2,  3,  4,  2,  9,  3, 10,  0,  4,  5,  9, 10, 11,  0,  5,  8,  9, 10, 11,  0,  1,  6,  0,  1,  7,  0,  1,  8,  9, 16, 10, 17,  0, 11, 12, 16, 17, 18,  0, 12, 15, 16, 17, 18,  0,  5,  6,  8, 13,  0,  5,  7,  8, 14,  0,  5,  8, 15, 16, 23, 17, 24,  0, 18, 19, 23, 24, 25,  0, 19, 22, 23, 24, 25,  0, 12, 13, 15, 20,  0, 12, 14, 15, 21,  0, 12, 15, 22, 23, 30, 24, 31,  0, 25, 26, 30, 31, 32,  0, 26, 29, 30, 31, 32,  0, 19, 20, 22, 27,  0, 19, 21, 22, 28,  0, 19, 22, 29, 30, 37, 31, 38,  0, 32, 33, 37, 38, 39,  0, 33, 36, 37, 38, 39,  0, 26, 27, 29, 34,  0, 26, 28, 29, 35,  0, 26, 29, 36, 37, 44, 38, 45,  0, 39, 40, 44, 45, 46,  0, 40, 43, 44, 45, 46,  0, 33, 34, 36, 41,  0, 33, 35, 36, 42,  0, 33, 36, 43, 44, 51, 45, 52,  0, 46, 47, 51, 52, 53,  0, 47, 50, 51, 52, 53,  0, 40, 41, 43, 48,  0, 40, 42, 43, 49,  0, 40, 43, 50, 51, 58, 52, 59,  0, 53, 54, 58, 59, 60,  0, 54, 57, 58, 59, 60,  0, 47, 48, 50, 55,  0, 47, 49, 50, 56,  0, 47, 50, 57, 58, 65, 59, 66,  0, 60, 61, 65, 66, 67,  0, 61, 64, 65, 66, 67,  0, 54, 55, 57, 62,  0, 54, 56, 57, 63,  0, 54, 57, 64, 67,  0, 61, 62, 64,  0, 61, 63, 64,  0, 61, 64, 68]
    Jcol = [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, 30, 30, 31, 31, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, 36, 37, 37, 38, 38, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 45, 45, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 50, 50, 50, 50, 51, 51, 52, 52, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 57, 57, 57, 57, 58, 58, 59, 59, 60, 60, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 64, 64, 64, 64, 65, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, 68, 68]
    
    N = 11
    x = linspace(x0,x1,N)
    y = linspace(y0,y1,N)
    v = linspace(0,sqrt((x1-x0)/(h0*(N-1))**2 + (y1-y0)/(h0*(N-1))**2),N)
    px = ones(N)
    py = ones(N)
    pv = ones(N)
    u = zeros(N)
    px[0] = nan
    py[0] = nan
    pv[0] = nan
    u[-1] = nan
    
    Sol0 = array([x,y,v,u,px,py,pv])
    Sol0 = append(h0,Sol0.flatten('F'))
    Z0 = Sol0[idx]
    #pdb.set_trace()
    #Z = fsolve(brachfun, Z0, args=(x, y, v, g, Jrow, Jcol), fprime=brachdfun, full_output=0, xtol=1e-10)
    Z = root(brachfun, Z0, args=(x, y, v, g, Jrow, Jcol), jac=brachdfun, method='lm')
    
    Sol = Sol0
    Sol[idx] = Z.x

    return Sol[1:].reshape((7, N), order='F')