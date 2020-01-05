#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 15:44:32 2018

@author: robert
"""

import numpy as np
from scipy.optimize import newton


class Car:
    """ A car class to hold all the information about a car. """
    keys = [
            # Weights of parts of the car
            'W_uf',  #Front unsprung weight (N)
            'W_ur',  #Rear unsprung weight  (N)
            'W_s',   #Sprung weight         (N)
            # CG locations
            'z_wf',  #Front unsprung height (m)
            'z_wr',  #Rear unsprung height  (m)
            'h_s',   #Sprung height         (m)
            'a_s',   #Sprung distance from front axle (m)
            # Car geometry
            'l',     #Wheelbase   (m)
            't_f',   #Front track (m)
            't_r',   #Rear track  (m)
            # Roll center heights
            'z_rf',  #Front roll center height (m)
            'z_rr',  #Rear roll center height  (m)
            # Roll rates
            'K_f',   #Front roll rate (N-m/rad)
            'K_r',   #Rear roll rate  (N-m/rad)
            # Steering
            'ack',   #Steering Ackerman (%)
            'toe_f', #Front toe (deg)
            'toe_r', #Rear toe  (deg)
            'tire'   #Instance of a tire class
            ]
    
    # Initialization functions
    #--------------------------------------------------------------------------
    def __init__(self, params):
        self.params = params
        self.check_params(params)
        self.params_to_attrs(params)
        
        self.b_s = self.l - self.a_s
        self.theta, self.h_2 = self.calc_roll_axis()
        print('Roll axis inclination: %0.2f deg' % (self.theta*180/np.pi))
        print('Unsprung CG height from roll axis: %0.2f m' % self.h_2)
        
        self.phiA = self.calc_roll_sensitivity()
        print('Roll sensitivity: %0.2f deg/g' % (self.phiA*180/np.pi))
        
        self.WT_yf, self.WT_yr= self.calc_lateral_weight_transfer()
        print('Front lateral weight transfer: %0.2f N/g' % self.WT_yf)
        print('Rear lateral weight transfer:  %0.2f N/g' % self.WT_yr)
        
        self.WT_x = self.calc_longitudinal_weight_transfer()
        print('Longitudinal weight transfer:       %0.2f N/g' % self.WT_x)
        
        self.W_f0, self.W_r0 = self.calc_static_load()
        print('Front static wheel load:            %0.2f N' % self.W_f0)
        print('Rear static wheel load:             %0.2f N' % self.W_r0)
        
        self.W = 2*self.W_f0 + 2*self.W_r0
        print('Total weight:     %0.2f N' % self.W)
        
    
    def check_params(self, params):
        """ Check to ensure all required keys are in the params dictionary. """
        for key in self.keys:
            if key not in params:
                raise ValueError('The params object has no key: %s' % key)
    
    def params_to_attrs(self, params):
        """ Add all params as attributes of the class. """
        for key in params:
            setattr(self, key, params[key])
    
    # Calculation functions
    #--------------------------------------------------------------------------
    def calc_roll_axis(self):
        """ Calculate the roll axis angle (rad) and unsprung distance (m). """
        z_rr = self.z_rr
        z_rf = self.z_rf
        l    = self.l
        theta = np.arctan((z_rr-z_rf) / l)
        h_2 = (self.h_s - self.a_s*(z_rr-z_rf)/l - z_rf)*np.cos(theta)
        return theta, h_2
    
    def calc_roll_sensitivity(self):
        """ Calculate the roll sensitivity (rad/g). """
        h_2 = self.h_2
        W_s = self.W_s
        phiA = -W_s*h_2 / (self.K_f+self.K_r-W_s*h_2)
        return phiA
    
    def calc_lateral_weight_transfer(self):
        """ Calculate lateral weight transfer (N/g). """
        W_s =  self.W_s
        l   =  self.l
        h_2 =  self.h_2
        a_s =  self.a_s
        b_s =  self.b_s
        phiA = self.phiA
        # Front
        K_fp = self.K_f - b_s*W_s*h_2/l
        WT_yf = (-phiA*K_fp + W_s*b_s*self.z_rf/l+self.W_uf*self.z_wf)/self.t_f
        # Rear
        K_rp = self.K_r - a_s*W_s*h_2/l
        WT_yr = (-phiA*K_rp + W_s*a_s*self.z_rr/l+self.W_ur*self.z_wr)/self.t_r
        return WT_yf, WT_yr
    
    def calc_longitudinal_weight_transfer(self):
        """ Calculate longitudinal weight transfer (N/g). """
        Wh = self.h_s*self.W_s+self.z_wf*self.W_uf+self.z_wr*self.W_ur
        WT_x = Wh/self.l
        return WT_x
    
    def calc_static_load(self):
        """ Calculate static wheel loads (N). """
        W_s = self.W_s
        l   = self.l
        W_f0 = -(self.W_uf + W_s*self.b_s/l)/2
        W_r0 = -(self.W_ur + W_s*self.a_s/l)/2
        return W_f0, W_r0
    
    # Equation functions
    #--------------------------------------------------------------------------
    def load_from_A(self, Ax, Ay):
        """ Calculate the wheel loads (N) for a given acceleration (g). """
        W_f0  = self.W_f0
        W_r0  = self.W_r0
        WT_x  = self.WT_x
        WT_yf = self.WT_yf
        WT_yr = self.WT_yr
        # Calculate wheel loads
        W_fr = W_f0 + Ax*WT_x + Ay*WT_yf
        W_fl = W_f0 + Ax*WT_x - Ay*WT_yf
        W_rr = W_r0 - Ax*WT_x + Ay*WT_yr
        W_rl = W_r0 - Ax*WT_x - Ay*WT_yr
        return W_fr, W_fl, W_rr, W_rl
    
    def tire_angle_from_steer(self, delta):
        """ Calculate the steering angle of each tire from the steer angle. """
        l =   self.l
        t_f = self.t_f
        ack = self.ack
        # Calculate 100% Ackerman
        if delta == 0.0:
            return 0.0, 0.0
        delta_i = np.arctan(l / (l/np.tan(abs(delta))-t_f/2))
        delta_o = np.arctan(l / (l/np.tan(abs(delta))+t_f/2))
        # Apply Ackerman
        delta_i = ack*(delta_i - delta) + delta
        delta_o = ack*(delta_o - delta) + delta
        return delta_i, delta_o
    
    def radius_induced_slip(self, beta, r):
        """ Calculate the slip angles induced from going around a radius. """
        a_s = self.a_s
        b_s = self.b_s
        t_f = self.t_f
        t_r = self.t_r
        # Front
        xi_fr = np.arctan((a_s*np.cos(beta)+t_f*np.sin(beta)/2) 
                / (r+a_s*np.sin(beta)-t_f*np.cos(beta)/2))
        xi_fl = np.arctan((a_s*np.cos(beta)-t_f*np.sin(beta)/2) 
                / (r+a_s*np.sin(beta)+t_f*np.cos(beta)/2))
        # Rear
        xi_rr = -np.arctan((b_s*np.cos(beta)-t_r*np.sin(beta)/2) 
                / (r-b_s*np.sin(beta)-t_r*np.cos(beta)/2))
        xi_rl = -np.arctan((b_s*np.cos(beta)+t_r*np.sin(beta)/2) 
                / (r-b_s*np.sin(beta)+t_r*np.cos(beta)/2))
        return xi_fr, xi_fl, xi_rr, xi_rl
    
    def kinematic_slip(self, beta, delta):
        """ Calculate the kinematic angles of the tires. """
        toe_f = np.pi*self.toe_f/180
        toe_r = np.pi*self.toe_r/180
        # Steering contribution
        delta_i, delta_o = self.tire_angle_from_steer(delta)
        if delta >=0:
            delta_r = delta_i
            delta_l = delta_o
        else:
            delta_r = -delta_o
            delta_l = -delta_i
        KA_fr = beta + toe_f - delta_r
        KA_fl = beta - toe_f - delta_l
        KA_rr = beta + toe_r
        KA_rl = beta - toe_r
        return KA_fr, KA_fl, KA_rr, KA_rl
    
    def get_slip_angles(self, beta, r, KA):
        """ Calculate the slip angle for each tire. """
        # Radius contribution
        xi_fr, xi_fl, xi_rr, xi_rl = self.radius_induced_slip(beta, r)
        # Total slip angles
        SA_fr = KA[0] + xi_fr
        SA_fl = KA[1] + xi_fl
        SA_rr = KA[2] + xi_rr
        SA_rl = KA[3] + xi_rl
        return SA_fr, SA_fl, SA_rr, SA_rl
    
    def calc_peak_slip_angles(self, load):
        """ Calculate the optimal slip angle for each tire. """
        SApeak = []
        for i in range(4):
            muy = self.tire.get_muy(load[i])
            cz = self.tire.get_cz(load[i])
            SApeak.append(self.tire.get_SA(muy, cz, self.tire.peakSAbar))
        return SApeak[0], SApeak[1], SApeak[2], SApeak[3]
        
    def get_tire_forces(self, Fx, load, SA):
        """ Calculate the force produced by each tire. """
        tire = self.tire
        F_fr = tire.get_force(load[0], 180*SA[0]/np.pi, 0, 0)
        F_fl = tire.get_force(load[1], 180*SA[1]/np.pi, 0, 0)
        F_rr = tire.get_force(load[2], 180*SA[2]/np.pi, 0, Fx/2)
        F_rl = tire.get_force(load[3], 180*SA[3]/np.pi, 0, Fx/2)
        F_rr[0] = Fx/2
        F_rl[0] = Fx/2
        return F_fr, F_fl, F_rr, F_rl
    
    def get_tire_moments(self, Fx, load, SA):
        """ Calculate the moments produced by the tire. """
        tire = self.tire
        M_fr = tire.get_moment(load[0], 180*SA[0]/np.pi, 0, 0)
        M_fl = tire.get_moment(load[1], 180*SA[1]/np.pi, 0, 0)
        M_rr = tire.get_moment(load[2], 180*SA[2]/np.pi, 0, Fx/2)
        M_rl = tire.get_moment(load[3], 180*SA[3]/np.pi, 0, Fx/2)
        return M_fr, M_fl, M_rr, M_rl
    
    def get_effective_tire_forces(self, F, KA):
        """ Calculate the effective tire forces seen by the car. """
        FE = []
        for i in range(4):
            FE_x =  F[i][0]*np.cos(KA[i]) + F[i][1]*np.sin(KA[i])
            FE_y = -F[i][0]*np.sin(KA[i]) + F[i][1]*np.cos(KA[i])
            FE.append((FE_x, FE_y))
        return FE[0], FE[1], FE[2], FE[3]
    
    def get_chassis_tire_forces(self, F, KA, beta):
        """ Calculate the tire forces in the chassis frame. """
        FC = []
        for i in range(4):
            FC_x =  F[i][0]*np.cos(KA[i]-beta) + F[i][1]*np.sin(KA[i]-beta)
            FC_y = -F[i][0]*np.sin(KA[i]-beta) + F[i][1]*np.cos(KA[i]-beta)
            FC.append((FC_x, FC_y))
        return FC[0], FC[1], FC[2], FC[3]
    
    def calc_total_lateral_force(self, FE):
        """ Calculate the total effective lateral force. """
        return FE[0][1] + FE[1][1] + FE[2][1] + FE[3][1]
    
    def calc_total_longitudinal_force(self, FE):
        """ Calculate the total effective lateral force. """
        return FE[0][0] + FE[1][0] + FE[2][0] + FE[3][0]
    
    def calc_yaw_moment(self, FC):
        """ Calculate the total yaw moment about the CG. """
        t_f = self.t_f
        t_r = self.t_r
        a_s = self.a_s
        b_s = self.b_s
        Mz  = -FC[0][0]*t_f/2 + FC[0][1]*a_s
        Mz +=  FC[1][0]*t_f/2 + FC[1][1]*a_s
        Mz += -FC[2][0]*t_r/2 - FC[2][1]*b_s
        Mz +=  FC[3][0]*t_r/2 - FC[3][1]*b_s
        return Mz
    
    # Simulation functions
    #--------------------------------------------------------------------------
    def Ay_Fx_from_angle(self, r, beta, delta, N=100, Ay=1.0, Fx=0, extra=False):
        """ For the passed angles and radius, iteratively find Ay and Fx. """
        for k in range(100):
            load = self.load_from_A(-Ay*np.sin(beta), Ay*np.cos(beta))
            KA = self.kinematic_slip(beta, delta)
            SA = self.get_slip_angles(beta, r, KA)
            F = self.get_tire_forces(Fx, load, SA)
            FE = self.get_effective_tire_forces(F, KA)
            Fx += -self.calc_total_longitudinal_force(FE)
            Ay = -self.calc_total_lateral_force(FE) / self.W
        if extra:
            return Ay, Fx, load, KA, SA, F, FE
        return Ay, Fx
    
    def Mx_from_delta(self, delta, r, beta):
        """ Find the yaw moment for a given delta, r, and beta. """
        Ay, Fx, load, KA, SA, F, FE = self.Ay_Fx_from_angle(r, beta, delta,
                                                            extra=True)
        FC = self.get_chassis_tire_forces(F, KA, beta)
        Mz = self.calc_yaw_moment(FC)
        return Mz
    
    def Mx_from_beta(self, beta, r, delta):
        """ Find the yaw moment for a given delta, r, and beta. """
        Ay, Fx, load, KA, SA, F, FE = self.Ay_Fx_from_angle(r, beta, delta,
                                                            extra=True)
        FC = self.get_chassis_tire_forces(F, KA, beta)
        Mz = self.calc_yaw_moment(FC)
        return Mz
    
    def get_delta(self, r, beta, delta0=0.2):
        """ Calculate the steer angle needed for a given beta and radius. """
        try:
            delta = newton(self.Mx_from_delta, delta0, args=(r, beta))
        except RuntimeError:
            delta = float('nan')
        if abs(delta) > 0.5:
            delta = float('nan')
        return delta
    
    def get_beta(self, r, delta, beta0=0.2):
        """ Calculate the steer angle needed for a given beta and radius. """
        try:
            beta = newton(self.Mx_from_beta, beta0, args=(r, delta))
        except RuntimeError:
            beta = float('nan')
        if abs(delta) > 0.5:
            beta = float('nan')
        return beta
    
    # Output functions
    #--------------------------------------------------------------------------
    def print_force_info(self, Ay, beta, delta, Mz, KA, SA, SAP, load, F, FE, FC):
        """ Nicely display the important numbers. """
        con = 180/np.pi
        print('Values at the maximum')
        print('---------------------------------------------------------------------')
        print('Ay:              %0.2f g' % Ay)
        print('Beta:            %0.2f deg' % (beta*con))
        print('Delta:           %0.2f deg' % (delta*con))
        print('Mz:              %0.2f N-m' % Mz)
        print('')
        print('Wheel parameters             FR          FL          RR          RL')
        print('---------------------------------------------------------------------')
        print('Kin angle (deg):     %10.2f, %10.2f, %10.2f, %10.2f' % (KA[0]*con, KA[1]*con, KA[2]*con, KA[3]*con))
        print('Slip angle (deg):    %10.2f, %10.2f, %10.2f, %10.2f' % (SA[0]*con, SA[1]*con, SA[2]*con, SA[3]*con))
        print('Peak angle (deg):    %10.2f, %10.2f, %10.2f, %10.2f' % SAP)
        print('')
        print('Wheel loads (N):     %10.2f, %10.2f, %10.2f, %10.2f' % load)
        print('Tire lon forces (N): %10.2f, %10.2f, %10.2f, %10.2f' % (F[0][0], F[1][0], F[2][0], F[3][0]))
        print('Tire lat forces (N): %10.2f, %10.2f, %10.2f, %10.2f' % (F[0][1], F[1][1], F[2][1], F[3][1]))
        print('Eff lon forces (N):  %10.2f, %10.2f, %10.2f, %10.2f' % (FE[0][0], FE[1][0], FE[2][0], FE[3][0]))
        print('Eff lat forces (N):  %10.2f, %10.2f, %10.2f, %10.2f' % (FE[0][1], FE[1][1], FE[2][1], FE[3][1]))
        print('Cha lon forces (N):  %10.2f, %10.2f, %10.2f, %10.2f' % (FC[0][0], FC[1][0], FC[2][0], FC[3][0]))
        print('Cha lat forces (N):  %10.2f, %10.2f, %10.2f, %10.2f' % (FC[0][1], FC[1][1], FC[2][1], FC[3][1]))
        print('')
        print('Front lat force (N): %10.2f' % (FE[0][1] + FE[1][1]))
        print('Rear lat force (N):  %10.2f' % (FE[2][1] + FE[3][1]))
    
