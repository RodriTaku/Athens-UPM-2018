# -*- coding: utf-8 -*-
"""
Created on March, 2018

@author: RodrigoTakuro
"""
#from math import pi
from numpy import array, matrix, zeros, append, sin, cos, linspace, pi
from scipy.optimize import fsolve
from scipy.sparse import csr_matrix
from parameterStructure import paramStruct

def cartPendfun(Y, X, T, Ffun, Ufun, dUfun, param, Jrow, Jcol, h):
    """Cart & Pendulum symplectic integrator (trapezoidal rule)
    
    Parameters
    ----------
    Y : numpy.array
        new internal state (x_{k+1}, theta_{k+1}, vx_{k+1}, vtheta_{k+1},
        VX_{k+1}_1, VX_{k+1}_2, VTheta_{k+1}_1, VTheta_{k+1}_2)
    X : numpy.array
        current external state
        (x_k, theta_k, vx_k, vtheta_k)
    T : numpy.array
        current and new time (t_k, t_{k+1})
    Ffun : function
        forcing function F(t)
    Ufun : function
        control function U(t,x,theta,vx,vtheta)
    dUfun : function
        control function gradient dU(t,x,theta,vx,vtheta)
    param : paramStruct class object
        parameters of the model
    Jrow : list or numpy.array
        row indices of nonzero entries of the jacobian
    Jcol : list or numpy.array
        column indices of nonzero entries of the jacobian
    h : float
        step size
        
    Returns
    -------
    Fun : numpy.array
        Integrator equations to be solved
    """
    F = array([Ffun(T[0]), Ffun(T[1])])
    U = array([Ufun(T[0], X[0], X[1], Y[4], Y[6]), Ufun(T[1], Y[0], Y[1], Y[5], Y[7])])
    return array([Y[0] - X[0] - (Y[4]*h)/2 - (Y[5]*h)/2,
                Y[1] - X[1] - (Y[6]*h)/2 - (Y[7]*h)/2,
                Y[4]*(param.M + param.m) - X[2]*(param.M + param.m) - (h*(F[0] + U[0]))/2 - X[3]*param.L*param.m*cos(X[1]) + Y[6]*param.L*param.m*cos(X[1]),
                Y[5]*(param.M + param.m) - X[2]*(param.M + param.m) - (h*(F[0] + U[0]))/2 - X[3]*param.L*param.m*cos(X[1]) + Y[7]*param.L*param.m*cos(Y[1]),
                2*Y[2]*(param.M/2 + param.m/2) - 2*X[2]*(param.M/2 + param.m/2) - (h*(F[0] + U[0]))/2 - (h*(F[1] + U[1]))/2 - X[3]*param.L*param.m*cos(X[1]) + Y[3]*param.L*param.m*cos(Y[1]),
                (h*param.L*(param.g*param.m*sin(X[1]) - cos(X[1])*F[0] + Y[4]*Y[6]*param.m*sin(X[1])))/2 - X[3]*param.L**2*param.m + Y[6]*param.L**2*param.m - X[2]*param.L*param.m*cos(X[1]) + Y[4]*param.L*param.m*cos(X[1]),
                (h*param.L*(param.g*param.m*sin(X[1]) - cos(X[1])*F[0] + Y[4]*Y[6]*param.m*sin(X[1])))/2 - X[3]*param.L**2*param.m + Y[7]*param.L**2*param.m - X[2]*param.L*param.m*cos(X[1]) + Y[5]*param.L*param.m*cos(Y[1]),
                (param.L*(2*Y[3]*param.L*param.m - h*cos(Y[1])*F[1] - 2*X[3]*param.L*param.m - h*cos(X[1])*F[0] - 2*X[2]*param.m*cos(X[1]) + 2*Y[2]*param.m*cos(Y[1]) + param.g*h*param.m*sin(X[1]) + param.g*h*param.m*sin(Y[1]) + Y[4]*Y[6]*h*param.m*sin(X[1]) + Y[5]*Y[7]*h*param.m*sin(Y[1])))/2])

    
def cartPenddfun(Y, X, T, Ffun, Ufun, dUfun, param, Jrow, Jcol, h):
    """Cart & Pendulum symplectic integrator jacobian (trapezoidal rule)
    Only the nonzero elements of the jacobian are computed and afterwards
    a sparse array is assembled. SciPy solvers do not yet allow sparse
    arrays as jacobians, so the array is turned to dense.
    
    Parameters
    ----------
    Y : numpy.array
        new internal state (x_{k+1}, theta_{k+1}, vx_{k+1}, vtheta_{k+1},
        VX_{k+1}_1, VX_{k+1}_2, VTheta_{k+1}_1, VTheta_{k+1}_2)
    X : numpy.array
        current external state
        (x_k, theta_k, vx_k, vtheta_k)
    T : numpy.array
        current and new time (t_k, t_{k+1})
    Ffun : function
        forcing function F(t)
    Ufun : function
        control function U(t,x,theta,vx,vtheta)
    dUfun : function
        control function gradient dU(t,x,theta,vx,vtheta)
    param : paramStruct class object
        parameters of the model
    Jrow : list or numpy.array
        row indices of nonzero entries of the jacobian
    Jcol : list or numpy.array
        column indices of nonzero entries of the jacobian
    h : float
        step size
        
    Returns
    -------
    dFun : numpy.array
        jacobian of cartPendfun with respect to Y
    """
    F = array([Ffun(T[0]), Ffun(T[1])])
    dU = array([dUfun(T[0], X[0], X[1], Y[4], Y[6]), dUfun(T[1], Y[0], Y[1], Y[5], Y[7])])
    Jvec = array([1,
                -(h*dU[1,1])/2,
                1,
                -Y[7]*param.L*param.m*sin(Y[1]),
                - (h*dU[1,2])/2 - Y[3]*param.L*param.m*sin(Y[1]),
                -Y[5]*param.L*param.m*sin(Y[1]),
                (param.L*(h*F[1]*sin(Y[1]) - 2*Y[2]*param.m*sin(Y[1]) + param.g*h*param.m*cos(Y[1]) + Y[5]*Y[7]*h*param.m*cos(Y[1])))/2,
                param.M + param.m,
                param.L*param.m*cos(Y[1]),
                param.L*param.m*cos(Y[1]),
                param.L**2*param.m,
                -h/2,
                param.M + param.m - (h*dU[0,3])/2,
                -(h*dU[0,3])/2,
                -(h*dU[0,3])/2,
                param.L*param.m*cos(X[1]) + (Y[6]*h*param.L*param.m*sin(X[1]))/2,
                (Y[6]*h*param.L*param.m*sin(X[1]))/2,
                (Y[6]*h*param.L*param.m*sin(X[1]))/2,
                -h/2,
                param.M + param.m,
                -(h*dU[1,3])/2,
                param.L*param.m*cos(Y[1]),
                (Y[7]*h*param.L*param.m*sin(Y[1]))/2,
                -h/2,
                param.L*param.m*cos(X[1]) - (h*dU[0,4])/2,
                -(h*dU[0,4])/2,
                -(h*dU[0,4])/2,
                param.m*param.L**2 + (Y[4]*h*param.m*sin(X[1])*param.L)/2,
                (Y[4]*h*param.L*param.m*sin(X[1]))/2,
                (Y[4]*h*param.L*param.m*sin(X[1]))/2,
                -h/2,
                param.L*param.m*cos(Y[1]),
                -(h*dU[1,4])/2,
                param.L**2*param.m,
                (Y[5]*h*param.L*param.m*sin(Y[1]))/2])
    return csr_matrix((Jvec, (Jrow, Jcol))).todense()


def cartPendsim( t0, q0, v0, ffun, ufun, dufun, param, h, nsteps ):
    """Cart & Pendulum simulation run using a symplectic integrator
    (trapezoidal rule)
    
    Parameters
    ----------
    t0 : int, float
        initial time (t_0)
    q0 : numpy.array
        initial configuration (x_0, theta_0)
    v0 : numpy.array
        initial velocities (vx_0, vtheta_0)
    ffun : function
        forcing function F(t). Must return single value
    ufun : function
        control function U(t,x,theta,vx,vtheta). Must return single value
    dufun : function
        gradient of control function dU(t,x,theta,vx,vtheta).
        Must return numpy.array of dimension 5
    param : paramStruct class object
        parameters of the model
    h : float
        step size
    nsteps : int
        number of steps of the iteration
    
    Returns
    -------
    sol : dictionary
        Container with time vector T, X array with states, F and U arrays with
        internal values of forcing and control
    """
    #Load row and column indices to reconstruct sparse jacobian
    Jrow = [ 0, 4, 1, 3, 4, 6, 7, 4, 7, 4, 7, 0, 2, 3, 4, 5, 6, 7, 0, 3, 4, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 3, 4, 6, 7]
    Jcol = [ 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7]
    
    #Initialize time vector container T
    T = t0 + linspace(0,nsteps,nsteps+1)*h    
    #Initialize external state container X:
    #The k-th column represents the k-th time-slice
    #(x_k, theta_k, vx_k, vtheta_k)
    X = zeros([4,nsteps+1])
    X[:,0] = append(q0,v0)
    #Initialize force container F:
    #The k-th column represents the values of the internal stages
    #(F_k^1, F_k^2)
    F = zeros([2,nsteps+1])
    #Initialize control container Y:
    #The k-th column represents the values of the internal stages
    #(U_k^1, U_k^2)
    U = zeros([2,nsteps+1])
    indsXY = range(0,4)
    #Initialize internal state container Y:
    #The k-th column represents the k-th time-slice
    #(x_k, theta_k, vx_k, vtheta_k, VX_k_1, VX_k_2, VTheta_k_1, VTheta_k_2)
    Y = zeros(8)
    Y[indsXY] = X[:,0]
    
    #Main loop
    for j in range(0, nsteps):
        #Call fsolve with 'cartPendfun' and 'cartPenddfun' as main inputs
        Y = fsolve(cartPendfun, Y, args=(X[:,j], T[j:j+2], ffun, ufun, dufun, param, Jrow, Jcol, h), fprime=cartPenddfun, full_output=0, xtol=1e-10)
        #Store results in containers
        F[:,j+1] = array([ffun(T[j]), ffun(T[j+1])])
        U[:,j+1] = array([ufun(T[j], X[0,j], X[1,j], Y[4], Y[6]), ufun(T[j+1], Y[0], Y[1], Y[5], Y[7])])
        X[:,j+1] = Y[indsXY]
    
    #return values in a dictionary
    return {'T':T, 'X':X, 'F':F ,'U':U }


def cartPendlinearization( X, param ):
    """Cart & Pendulum linearization
    
    Parameters
    ----------
    X : numpy.array
        state around which we will linearize
        (x, theta, vx, vtheta)
    param : paramStruct class object
        parameters of the model
        
    Returns
    -------
    A : numpy.matrix
        state linearization matrix
    B : numpy.matrix
        control linearization matrix
    """
    #x = X[0]
    theta = X[1]
    #vx = X[2]
    vtheta = X[3]
    
    M = param.M
    m = param.m
    g = param.g
    L = param.L

    A = matrix([[ 0,                                                                                                                                                                                                                                                            0, 1,                                                      0],
                [ 0,                                                                                                                                                                                                                                                            0, 0,                                                      1],
                [ 0,                                                          (m*(2*g*cos(theta)**2 - g + L*vtheta**2*cos(theta)))/(M + m - m*cos(theta)**2) - (m**2*cos(theta)*sin(theta)*(2*L*sin(theta)*vtheta**2 + 2*g*cos(theta)*sin(theta)))/(- m*cos(theta)**2 + M + m)**2, 0,        (2*L*m*vtheta*sin(theta))/(m*sin(theta)**2 + M)],
                [ 0, (2*m*cos(theta)*sin(theta)*(L*m*cos(theta)*sin(theta)*vtheta**2 + g*m*sin(theta) + M*g*sin(theta)))/(L*(- m*cos(theta)**2 + M + m)**2) - (2*(L*m*(2*cos(theta)**2 - 1)*vtheta**2 + g*m*cos(theta) + M*g*cos(theta)))/(L*(2*M + m - m*(2*cos(theta)**2 - 1))), 0, -(2*m*vtheta*sin(2*theta))/(2*M + m - m*cos(2*theta))]])
     
    B = matrix([                                        0,
                                                        0,
                              1/(M + m - m*cos(theta)**2),
                -cos(theta)/(L*(M + m - m*cos(theta)**2))])
    return (A, B)