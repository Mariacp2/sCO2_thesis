class Config:
    
    
    heat_exchanger_efficiency = .85
    turbine_isentropic_efficiency = .85
    pump_isentropic_efficiency = .80
    Pmax = 35*1e6
    T1 = 273 + 30
    Tmax = 273 + 200

    rankine_Pmin = 7.19*1e6
    rankine_rp_max = 4.9 # = 35/7.19
    rankine_rp_min = 1.03 # Pcrit/7.19

    brayton_Pmin = 7.5*1e6
    brayton_rp_max = 4.67  # = 35/7.5
    brayton_rp_min = 1.5     