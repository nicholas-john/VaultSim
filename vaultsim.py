# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 20:15:57 2021

@author: Nicholas John
"""

import numpy as np
import matplotlib.pyplot as plt

from bvp import circular_section

import matplotlib.animation as animation

LITTLE_G = 9.81 # acceleration due to gravity

def trapezoid(y, dx):
    sum = y[:-1] + y[1:]
    return sum / (2 * dx)

def simple_energy(mass, velocity, y):
    assert y > 0, "y-position must be higher than reference height of zero"
    kinetic = .5 * mass * velocity@velocity
    potential = mass * LITTLE_G * y
    return kinetic + potential

class Vault:
    
    def __init__(self, 
                 L=4, # length of pole in meters
                 pointmass=70, # kilograms
                 initial_veloc=[10,5], # (x,y) m/s of point mass
                 initial_theta=np.pi-.4, # angle of pole from ground
                 npts=100,
                 dtime=.01,
                 ntime=200,
                 pole_dens=.5, # kilograms per meter
                 E=200, # Elasticity constant 
                 I=.01, # radius of pole's cross section
                 name="vault" # file name for the gif
                 ):
        
        self.L = L
        self.pointmass = pointmass
        self.initial_veloc = np.array(initial_veloc)
        self.npts = npts
        self.dtime = dtime
        self.ntime = ntime
        self.ds = L / (npts - 1)
        self.E = E
        self.I = I
        self.name = name
        
        # vector for beam
        self.s = np.linspace(0, self.L, self.npts)
        self.posx = self.s * np.cos(initial_theta)
        self.posy = self.s * np.sin(initial_theta)
        
        self.TOTAL_ENERGY = simple_energy(pointmass, 
                                          self.initial_veloc, 
                                          self.posy[-1])
        self.velocity = self.initial_veloc 
        
    def update_unit_tangent(self):
        """
        Finite difference calculation of unit tangent vector
        at index i of the pole position's discretization
        """
        gradx = np.gradient(self.posx, self.ds)
        grady = np.gradient(self.posy, self.ds)
        direction = np.array([gradx, grady])
        self.unit_tangent = direction / (gradx**2 + grady**2)**(1/2)
    
    def update_unit_normal(self):
        gradx = np.gradient(self.unit_tangent[0], self.ds)
        grady = np.gradient(self.unit_tangent[1], self.ds)
        direction = np.array([gradx, grady])
        self.unit_normal = direction / (gradx**2 + grady**2)**(1/2)
    
    def update_curvature(self):
        Tx = self.unit_tangent[:,0]
        Ty = self.unit_tangent[:,1]
        gradx = np.gradient(Tx, self.ds)
        grady = np.gradient(Ty, self.ds)
        self.curvature = (gradx**2 + grady**2)**(1/2)
    
    def update_chord(self):
        """
        The chord of the pole is the line connecting its two ends.
        The plug is always at the origin: (0, 0)
        This function computes the unit vector on the chord facing the cap.
        """
        chord = np.array([self.posx[-1], self.posy[-1]]) 
        self.chord_length = np.linalg.norm(chord)
        self.unit_chord = chord / self.chord_length 
        self.theta = np.arctan(self.unit_chord[1] / -self.unit_chord[0])
    
    def get_force(self):
        """
        Assuming the force is in the direction of the chord.
        Let's see if this assumption looks realistic.
        """
        integrand = self.E * self.I * self.curvature
        magnitude = trapezoid(integrand, self.ds)
        return magnitude * self.unit_chord
    
    def update_pole_position(self):
        """
        Solve a boundary value problem to find a position of the pole with 
        the desired energy and arc length that starts at (0, 0) and ends
        at (x_point, y_point)
        
        point_energy = simple_energy(self.pointmass, 
                                     self.velocity, 
                                     self.posy[-1])
        pole_energy = self.TOTAL_ENERGY - point_energy
        """
        # ^^^^^^^^ This is hard
        # so instead for now we assume the pole 
        # takes the arc of a circular section.
        # This assumption is probably reasonable for the 
        # point mass case, since it implies that the 
        # [stress -> strain -> curvature]
        # is evenly distributed over the length of the pole.
        # To build a model that works for multiple contact points, 
        # we need an energy-based formulation, which requires someone to
        # learn about mechanics of materials.
        
        # solve for circular arc of length L
        # starting at (0,0) and ending at (chord_length, 0)
        # to be reflected and rotated later
        arcx, arcy = circular_section(self.chord_length, self.L, self.npts)
        arcx = arcx - arcx[0]
        arcy = arcy - arcy[0]
        z = arcx + 1j * arcy
        z = z * np.exp( -1j * self.theta )
        #self.posx[:-1] = np.real(z)[:-1]
        #self.posy[:-1] = np.imag(z)[:-1]
        #self.posx[-1] = x_point
        #self.posy[-1] = y_point
        self.posx = np.real(z)
        self.posy = np.imag(z)
    
    def update_point_mass(self):
        """
        Euler's method --
        Could be upgraded to a Runge-Kutta method of higher order
        """
        force = self.get_force()
        # update position and velocity of point mass
        # with Euler's method:
        self.velocity = self.velocity + force * self.dtime / self.pointmass
        self.velocity[1] = self.velocity[1] - self.dtime * LITTLE_G
        x_point = self.posx[-1] + self.dtime * self.velocity[0]
        y_point = self.posy[-1] + self.dtime * self.velocity[1]
        self.posx[-1] = x_point
        self.posy[-1] = y_point
        
    def step(self):
        self.update_chord()
        if self.chord_length > self.L:
            return True # time to stop stepping
        self.update_unit_tangent()
        self.update_unit_normal()
        self.update_curvature()
        self.update_point_mass()
        self.update_chord()
        self.update_pole_position()
        return False
    
    def run(self):
        fig, ax = plt.subplots(figsize=(15, 12))
        ax.set_xlim([-self.L - 1, self.L + 1])
        ax.set_ylim([-1, self.L + 1])
        plt.axis('equal')
        plots = []
        time = 0
        stop = False
        while time < self.ntime and not stop:
            stop = self.step()
            plot = plt.plot(self.posx, self.posy,
                            color="grey",
                            linewidth=1,
                            )
            plot += plt.plot(self.posx + 2*self.I * self.unit_normal[0],
                            self.posy + 2*self.I * self.unit_normal[1],
                            color="grey",
                            linewidth=1)
            plot += plt.plot([self.posx[-1] + self.I * self.unit_normal[0,-1]],
                             [self.posy[-1] + self.I * self.unit_normal[1,-1]],
                             marker='o',
                             color="k")
            plot += plt.plot([self.posx[0] + self.I * self.unit_normal[0,0]], 
                             [self.posy[0] + self.I * self.unit_normal[1,0]],
                             marker='o',
                             color="grey")
            plots.append(plot)
            time += 1
        anim = animation.ArtistAnimation(fig, plots, 
                                         interval=1000*self.dtime, 
                                         blit=True,
                                         repeat_delay=3000)
        anim.save(self.name + ".gif")
        
def test():
    vault = Vault()
    vault.run()

test()