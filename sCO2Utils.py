from CoolProp.CoolProp import PropsSI
import CoolProp as CP
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Compute crtical pressure
CO2 = CP.AbstractState("HEOS", "CO2")
pc = CO2.keyed_output(CP.iP_critical)

class ThermodynamicState:
  def __init__(
      self,
      p = 1e5,
      t = 273,
      h = PropsSI('H','P',1e5,'T',273,'CO2')/1000,
      s = PropsSI('S','P',1e5,'T',273,'CO2')/1000,      
  ):
    """
    p (float) : Pressure of fluid (Pa)
    t (float) : Temperature of fluid (K)
    h (float) : Enthalpy of fluid (KJ/kg)
    s (float) : Entropy of fluid (KJ/kgK)
    """
    self.p = p
    self.t = t
    self.h = h
    self.s = s
  
  def compute_from_p_t(self):
    self.h = PropsSI('H','P',self.p,'T', self.t,'CO2')/1000 # to convert to KJ/kg
    self.s = PropsSI('S','P',self.p,'T', self.t,'CO2')/1000 # to convert to KJ/kg

  def compute_from_p_h(self):

    self.t = PropsSI('T','P',self.p,'H', self.h*1000,'CO2')
    self.s = PropsSI('S','P',self.p,'H', self.h*1000,'CO2')/1000 # to convert to KJ/kg
  
  def compute_from_p_s(self):
    
    self.t = PropsSI('T','P',self.p,'S', self.s*1000,'CO2')
    self.h = PropsSI('H','P',self.p,'S', self.s*1000,'CO2')/1000 # to convert to KJ/kg

  
  def to_dict(self):
    """
    Returns : Dictionary of the state
    """
    return {
        'pressure (Pa)' : self.p,
        'temperature (K)' : self.t,
        'enthalpy (KJ/kg)' : self.h,
        'entropy (KJ/kgK)' : self.s
    }


class Turbine:
  def __init__(
      self,
      isentropic_efficiency = .8,
      output_pressure = 2e6):
    """
    isentropic_efficiency (float): Isentropic efficiency of turbine
    output_pressure (float) : The output operating pressure (Pa) of the turbine 
    """
    self.isentropic_efficiency = isentropic_efficiency
    self.output_pressure = output_pressure

  def resolve_isentropic_output_state(
      self,
      input_state : ThermodynamicState):
    """
    input_state (ThermodynamicState) : The input state of the turbine
    Returns 
    output_state (ThermodynamicState) : The output state of the turbine if the process is fully isentropic
    """
    h_is = PropsSI('H','P',self.output_pressure,'S', input_state.s*1000,'CO2')/1000 # to convert to KJ/kg
    T_is = PropsSI('T','P',self.output_pressure,'S', input_state.s*1000,'CO2')

    return ThermodynamicState(
        p = self.output_pressure,
        t = T_is,
        h = h_is,
        s = input_state.s
    )

  def resolve_actual_output_state(
      self,
      input_state : ThermodynamicState):
    """
    input_state (ThermodynamicState) : The input state of the turbine
    Returns 
    output_state (ThermodynamicState) : The actual output state of the turbine

    """

    h_is = PropsSI('H','P',self.output_pressure,'S', input_state.s*1000,'CO2')/1000 # to convert to KJ/kg
    h_actual = input_state.h + (h_is - input_state.h)*self.isentropic_efficiency
    s_actual = PropsSI('S','P',self.output_pressure,'H', h_actual*1000,'CO2')/1000 # to convert to KJ/kgK
    T_actual = PropsSI('T','P',self.output_pressure,'H', h_actual*1000,'CO2')
    
    return ThermodynamicState(
        p = self.output_pressure,
        t = T_actual,
        h = h_actual,
        s = s_actual
    )

class Pump:
  def __init__(
      self,
      isentropic_efficiency = .8,
      output_pressure = 2e6):
    """
    isentropic_efficiency (float): Isentropic efficiency of pump
    output_pressure (float) : The output operating pressure (Pa) of the pump
    """
    self.isentropic_efficiency = isentropic_efficiency
    self.output_pressure = output_pressure

  def resolve_isentropic_output_state(
      self,
      input_state : ThermodynamicState):
    """
    input_state (ThermodynamicState) : The input state of the pump
    Returns 
    output_state (ThermodynamicState) : The output state of the pump if the process is isentropic

    """
    h_is = PropsSI('H','P',self.output_pressure,'S', input_state.s*1000,'CO2')/1000 # to convert to KJ/kg
    T_is = PropsSI('T','P',self.output_pressure,'S', input_state.s*1000,'CO2')

    return ThermodynamicState(
        p = self.output_pressure,
        t = T_is,
        h = h_is,
        s = input_state.s
    )

  def resolve_actual_output_state(
      self,
      input_state : ThermodynamicState):
    """
    input_state (ThermodynamicState) : The input state of the pump
    Returns 
    output_state (ThermodynamicState) : The actual output state of the pump

    """

    h_is = PropsSI('H','P',self.output_pressure,'S', input_state.s*1000,'CO2')/1000 # to convert to KJ/kg
    h_actual = input_state.h + (h_is - input_state.h)/self.isentropic_efficiency
    s_actual = PropsSI('S','P',self.output_pressure,'H', h_actual*1000,'CO2')/1000 # to convert to KJ/kgK
    T_actual = PropsSI('T','P',self.output_pressure,'H', h_actual*1000,'CO2')
    
    return ThermodynamicState(
        p = self.output_pressure,
        t = T_actual,
        h = h_actual,
        s = s_actual
    )

def show_states_table(states_list):
  _l = [state.to_dict() for state in states_list]
  df = pd.DataFrame(_l)
  df.index = range(1,len(_l)+1)
  print(df)


def plot_tvss_phase_transition_region(starting_t, title):

  CO2 = CP.AbstractState("HEOS", "CO2")
  Tc = CO2.keyed_output(CP.iT_critical)  
  N = 100
  temp_space = np.linspace(starting_t,Tc,N)

  sf_array = np.array([PropsSI('S','T', temp_space[i],'Q',0,'CO2') for i in range(N)])/1000 # Conversion to KJ/kgK
  sg_array = np.array([PropsSI('S','T', temp_space[i],'Q',1,'CO2') for i in range(N)])/1000 # Conversion to KJ/kgK
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('Entropy (KJ/kgK)')
  ax.set_ylabel('T (K)')
  ax.set_title(title)
  saturation_line = np.concatenate((sf_array,np.flip(sg_array)), axis = 0)
  temps_vals = np.concatenate((temp_space, np.flip(temp_space)), axis = 0)
  ax.plot(saturation_line, temps_vals, color ='blue', label = 'saturation line')

  return fig,ax

