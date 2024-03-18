from CoolProp.CoolProp import PropsSI
from scipy import optimize
import CoolProp as CP
import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
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


class HeatExchanger:
  def __init__(
      self,
      stream1_entry_state :ThermodynamicState,
      stream2_entry_state :ThermodynamicState,
      stream1_mass_flow_rate = 1,
      stream2_mass_flow_rate = 1,
      grading = 1,
      efficiency = .85):
    """
    cold_fluid_entry_state (ThermodynamicState) : This is the state of the cold fluid which will be warm up by virtue of heat exchanger
    warm_fluid_entry_state (ThermodynamicState) : This is the state of the warm fluid which will warm the cold fluid
    grading (float) : A grading of 1 means the Heat Exchanger is ideal, meaning the maximum amount of heat is transferred and temperatures equalize. A 0 grading means there is no Heat Exchanger at all
    efficiency (float) : the efficiency dictates how much of the heat of the warm fluid is transferred to the cold fluid, the rest being lost to the environment
    """

    # The heat exchanger is taken to be the counterflow Heat Exchanger in its simplest form
    # The grading gives a linear interpolation between the ideal perfect case and the no heat transfer case

    self.stream1_entry_state = stream1_entry_state
    self.stream1_mass_flow_rate = stream1_mass_flow_rate
    self.stream2_entry_state = stream2_entry_state
    self.stream2_mass_flow_rate = stream2_mass_flow_rate
    self.grading = grading
    self.efficiency = efficiency 
  
  def heat_exchanged(self, T_eq):

    # Starting T1 = self.stream1_entry_state.t
    # We need to define the integration spacing
    # Break up T_eq - T1 into 10 pieces
    N = 10
    t1_space = np.linspace(self.stream1_entry_state.t, T_eq, N)
    cp1_avg = np.mean(np.array(PropsSI('CPMOLAR','P',self.stream1_entry_state.p,'T',t1_space,'CO2'))/44.01) # TO CONVERT IT TO KJ/kgK
    delta_Q1 = cp1_avg*(T_eq - self.stream1_entry_state.t)*self.stream1_mass_flow_rate

    # Starting T2 = self.stream2_entry_state.t
    # We need to define the integration spacing
    # Break up T_eq - T2 into 10 pieces
    t2_space = np.linspace(self.stream2_entry_state.t, T_eq, N)
    
    if self.stream2_entry_state.p < pc:
        # Make sure to avoid discontinuity at Ts
        # At the same time, avoid discontinuity if and only if p < p_crit, meaning CO2 can undergo a phase transition
        Ts = PropsSI('T','P',self.stream2_entry_state.p,'Q',0,'CO2')
        # If t2_space has a value that is too close to Ts it will raise an error into PropsSI when computing CPMOLAR
        # Thus we change those values to something that does not raise errors
        t2_space[(t2_space > Ts*.999)&(t2_space<Ts*1.001)] = Ts*1.001
        # Discontinuity avoided
    
    cp2_avg = np.mean(np.array(PropsSI('CPMOLAR','P',self.stream2_entry_state.p,'T',t2_space,'CO2'))/44.01) # TO CONVERT IT TO KJ/kgK
    delta_Q2 = cp2_avg*(T_eq - self.stream2_entry_state.t)*self.stream2_mass_flow_rate

    return delta_Q2 + delta_Q1


  def find_equilibrium_temperature(self):

    sol = optimize.root_scalar(self.heat_exchanged, bracket=[self.stream1_entry_state.t, self.stream2_entry_state.t])
    
    return sol.root

  def resolve_exit_states(self):
    """
    This function takes in all input state parameters and outputs the final states based on the grading and efficiency of the Heat Exchanger
    Returns 
    output_state : The output states of the Heat Exchanger, first the cold fluid exit state followed by warm fluid exit state

    """

    t_eq = self.find_equilibrium_temperature()
    
    # We need to know which stream is warmer at the start
    if self.stream1_entry_state.t > self.stream2_entry_state.t:
      warm_entry_state,cold_entry_state = self.stream1_entry_state,self.stream2_entry_state
      warm_mass_flow_rate, cold_mass_flow_rate = self.stream1_mass_flow_rate, self.stream2_mass_flow_rate
    else:
      warm_entry_state,cold_entry_state = self.stream2_entry_state,self.stream1_entry_state
      warm_mass_flow_rate, cold_mass_flow_rate = self.stream2_mass_flow_rate, self.stream1_mass_flow_rate


    ideal_warm_exit_state = ThermodynamicState(
        p = warm_entry_state.p,
        t = t_eq,
    )
    
    ideal_warm_exit_state.compute_from_p_t()

    # The ideal heat exchanged is
    ideal_delta_Q = warm_mass_flow_rate*(ideal_warm_exit_state.h - warm_entry_state.h)

    # The heat exchanged is based on grading and efficiency
    actual_delta_Q = ideal_delta_Q*self.grading*self.efficiency

    # The actual enthalpy change is thus
    actual_enthalpy_change = actual_delta_Q/warm_mass_flow_rate

    #The actual exit state of the warm stream is thus
    warm_exit_state = ThermodynamicState(
        p = warm_entry_state.p,
        h = warm_entry_state.h + actual_enthalpy_change,
    )

    warm_exit_state.compute_from_p_h()

    delta_Q_2 = actual_delta_Q

    enthalpy_change_2 = delta_Q_2/cold_mass_flow_rate # Keep in mind this enthalpy change is negative

    cold_exit_state = ThermodynamicState(
      p = cold_entry_state.p,
      h = cold_entry_state.h - enthalpy_change_2,
    )

    cold_exit_state.compute_from_p_h()

    if self.stream1_entry_state.t > self.stream2_entry_state.t:
      stream1_exit_state,stream2_exit_state = warm_exit_state, cold_exit_state      
    else:
      stream2_exit_state,stream1_exit_state = warm_exit_state, cold_exit_state      
      
    return stream1_exit_state, stream2_exit_state