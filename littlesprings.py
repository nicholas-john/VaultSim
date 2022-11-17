# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 00:10:59 2022

@author: Nicholas John
"""

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.animation as animation

class Vault:
    
    def __init__(
            self,
            L=4, # length of pole in meters
            poledens=.1, # kg / m
            initial_veloc=[17,2], # initial velocity at top of pole
            initial_theta=np.pi-.4, # angle of pole from ground
            npts=100, # number of length elements in discretization of pole
            dtime=.01, # length of time step in simulation
            ntime=6000, # number of time steps in simulation
            stopping_tol=.001, # how close to unbent to stop simulation
            plot_interval=100, # how many time steps between plot frames
            E=1000, # Elasticity constant
            I=.03, # radius of pole's cross sesction
            name="vault" # file name for the gif
                 ):
        
        self.L = L
        self.initial_veloc = np.array(initial_veloc)
        self.npts = npts
        self.dtime = dtime
        self.ntime = ntime
        self.stopping_tol = stopping_tol
        self.plot_interval=plot_interval
        self.ds = L / (npts - 1)
        self.E = E
        self.I = I
        self.name = name
        
        # vector for beam
        self.s = np.linspace(0, self.L, self.npts)
        self.posx = self.s * np.cos(initial_theta)
        self.posy = self.s * np.sin(initial_theta)
        self.velocity = self.initial_veloc 
        
        