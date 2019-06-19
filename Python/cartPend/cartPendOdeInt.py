# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 19:05:11 2018

@author: RodrigoTakuro
"""
from numpy import array, cos, sin, pi
from parameterStructure import paramStruct

def cartPendvecField( Y, t, f, u, du, param ):
    """Cart & Pendulum vector field (ODE) Python's
    non-symplectic integrator in odeint.
    
    Parameters
    ----------
    Y : numpy.array
        current state (x, theta, vx, vtheta)
    Y : float
        current time
    f : function
        forcing function f(t). Must return single value
    u : function
        control function u(t,x,theta,vx,vtheta). Must return single value
    du : function
        gradient of control function du(t,x,theta,vx,vtheta).
        Must return numpy.array of dimension 5
    param : paramStruct class object
        parameters of the model
        
    Returns
    -------
    vecField : numpy.array
        vector field
    """
    return array([                                                                                                                                                                                                                                   Y[2],
																																																										                                                                       Y[3],
																					                              -(-param.m*param.L*sin(Y[1])*Y[3]**2-cos(Y[1])*sin(Y[1])*param.g*param.m+f(t)*cos(Y[1])**2-f(t)-u(t,Y[0],Y[1],Y[2],Y[3]))/((1-cos(Y[1])**2)*param.m+param.M),
                  (-cos(Y[1])*sin(Y[1])*param.L*param.m**2*Y[3]**2-param.M*sin(Y[1])*param.g*param.m-sin(Y[1])*param.g*param.m**2+f(t)*param.M*cos(Y[1])-cos(Y[1])*param.m*u(t,Y[0],Y[1],Y[2],Y[3]))/(param.m*param.L*((1-cos(Y[1])**2)*param.m+param.M))])

def cartPendjacobian( Y, t, f, u, du, param ):
    """Cart & Pendulum jacobian (ODE) for Python's
    non-symplectic integrator in odeint.
    
    Parameters
    ----------
    Y : numpy.array
        current state (x, theta, vx, vtheta)
    Y : float
        current time
    f : function
        forcing function f(t). Must return single value
    u : function
        control function u(t,x,theta,vx,vtheta). Must return single value
    du : function
        gradient of control function du(t,x,theta,vx,vtheta).
        Must return numpy.array of dimension 5
    param : paramStruct class object
        parameters of the model
        
    Returns
    -------
    jac : numpy.array
        jacobian of the vector field
    """
    dU = du(t,Y[0],Y[1],Y[2],Y[3])
    return array([[                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0,                                                                   1,                                                                                                                                   0],
                  [                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0,                                                                   0,                                                                                                                                   1],
                  [                        dU[1]/(param.M - param.m*(cos(Y[1])**2 - 1)),                                                                                                                                                                                                   (dU[2] + param.g*param.m*cos(Y[1])**2 - param.g*param.m*sin(Y[1])**2 + 2*f(t)*cos(Y[1])*sin(Y[1]) + Y[3]**2*param.L*param.m*cos(Y[1]))/(param.M - param.m*(cos(Y[1])**2 - 1)) - (2*param.m*cos(Y[1])*sin(Y[1])*(f(t) - f(t)*cos(Y[1])**2 + Y[3]**2*param.L*param.m*sin(Y[1]) + param.g*param.m*cos(Y[1])*sin(Y[1]) + u(t, Y[0], Y[1], Y[2], Y[3])))/(param.M - param.m*(cos(Y[1])**2 - 1))**2,                        dU[3]/(param.M - param.m*(cos(Y[1])**2 - 1)),                                                   (dU[4] + 2*Y[3]*param.L*param.m*sin(Y[1]))/(param.M - param.m*(cos(Y[1])**2 - 1))],
                  [ -(cos(Y[1])*dU[1])/(param.L*(param.M - param.m*(cos(Y[1])**2 - 1))), (2*cos(Y[1])*sin(Y[1])*(param.g*param.m**2*sin(Y[1]) - param.M*f(t)*cos(Y[1]) + param.m*cos(Y[1])*u(t, Y[0], Y[1], Y[2], Y[3]) + param.M*param.g*param.m*sin(Y[1]) + Y[3]**2*param.L*param.m**2*cos(Y[1])*sin(Y[1])))/(param.L*(param.M - param.m*(cos(Y[1])**2 - 1))**2) - (param.M*f(t)*sin(Y[1]) + param.g*param.m**2*cos(Y[1]) + param.m*cos(Y[1])*dU[2] - param.m*sin(Y[1])*u(t, Y[0], Y[1], Y[2], Y[3]) + param.M*param.g*param.m*cos(Y[1]) + Y[3]**2*param.L*param.m**2*cos(Y[1])**2 - Y[3]**2*param.L*param.m**2*sin(Y[1])**2)/(param.L*param.m*(param.M - param.m*(cos(Y[1])**2 - 1))), -(cos(Y[1])*dU[3])/(param.L*(param.M - param.m*(cos(Y[1])**2 - 1))), -(2*Y[3]*param.L*cos(Y[1])*sin(Y[1])*param.m**2 + cos(Y[1])*dU[4]*param.m)/(param.L*param.m*(param.M - param.m*(cos(Y[1])**2 - 1)))]])