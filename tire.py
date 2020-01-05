#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 21:38:52 2018

@author: rariniello
"""

import numpy as np


# First lets define tire force function that returns the force the tire is producing
class Tire:
    """ A tire class meant to be extended for each tire. 
    
    All variables use the SAE tire axis system.
    """
    def __init__(self):
        """ Init """
    
    def magic_formula(self, x, B, C, D, E):
        return D*np.sin(C*np.arctan(B*(1-E)*x + E*np.arctan(B*x)))
    
    def get_force(self, Fz, SA, IA, SR):
        """ Calculates the forces acting on the tire. """
        # Fz is negative in this case
        muy = self.get_muy(Fz)
        mux = self.get_mux(Fz)
        cz = self.get_cz(Fz)
        SAbar = self.get_SAbar(muy, cz, SA)
        Fy = self.magic_formula(SAbar, self.By, self.Cy, muy*Fz, self.Ey)
        Fx = self.magic_formula(SR, self.Bx, self.Cx, mux*Fz, self.Ex)
        # Combined slip effects
        if Fx == 0.0:
            phi = np.pi/2
        else:
            phi = np.arctan(abs(Fy)*mux/abs(Fx)/muy)
        multiY = abs(np.sin(phi)) * 0.65 # 0.65 is de-rating to account for smoothness of the test belt
        multiX = abs(np.cos(phi)) * 0.65
        multi = np.sqrt(1-(SR/(Fz*mux))**2) * 0.65
        Fy *= multi
        #Fy *= multiY
        #Fx *= multiX
        return [Fx, Fy]
    
    def get_moment(self, Fz, SA, IA, Fx):
        """ Calculates the moments acting on the tire. """
        muy = self.get_muy(Fz)
        cz = self.get_cz(Fz)
        Tz = self.get_Tz(Fz)
        SAbar = self.get_SAbar(muy, cz, SA)
        Mz = self.magic_formula(SAbar, self.BMz, 2.0, self.DMz, self.EMz)
        multi = Tz*muy*Fz * 0.65
        Mz *= multi
        Mx = 0.0
        return [Mx, Mz]
    
    def get_muy(self, Fz):
        """ Calculates the lateral coefficient of friction. """
        return -self.my*Fz + self.y0
    
    def get_mux(self, Fz):
        """ Calculate the longitudinal coefficent of friction. """
        return self.mxa*np.exp(Fz/self.mxb)+self.mxc
    
    def get_cz(self, Fz):
        """ Calculates the cornering stiffness coefficient. """
        return self.c0*np.exp(Fz/self.cf)
    
    def get_Tz(self, Fz):
        """ Calculate the pneumatic trail. """
        return self.T0*(np.exp(Fz/self.Tf)-1)
    
    def get_SAbar(self, muy, cz, SA):
        """ Calculates the non-dimensionalized slip angle. """
        return np.rad2deg(cz * np.tan(np.deg2rad(np.array(SA))) / muy)
    
    def get_SA(self, muy, cz, SAbar):
        """ Calculate the slip angle from SAbar. """
        return np.rad2deg(np.arctan(np.deg2rad(muy * SAbar / cz)))


class Hoosier10X7(Tire):
    """ A tire class representing the Hoosier 18.0X7.5-10 tire. """
    def __init__(self):
        """ Set the tire fit parameters. """
        # Lateral coefficient of friction
        self.my = -7.573e-4
        self.y0 = 3.386
        # Longitudinal coefficient of friction
        self.mxa = 3.0819763
        self.mxb = 358.67327541
        self.mxc = 2.87227048
        # Cornering stiffness coefficient
        self.c0 = 1.613
        self.cf = 886.3
        # Lateral force magic formula
        self.By = 0.761
        self.Cy = 1.35
        self.Ey = 0.093
        # Longitudinal force magic formula
        self.Bx = 14.08078091
        self.Cx = 1.47615282
        self.Ex = 0.88199564
        # Tire parameters
        self.R = 0.216 #m
        # Useful information
        self.peakSAbar = 2
        # Pnuematic trail
        self.T0 = 7.43810407e-2
        self.Tf = 7.55013968e2
        # Self aligning torque magic formula
        self.BMz = 0.81205238
        self.DMz = 0.5371473
        self.EMz = -3.01115314
        