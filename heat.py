import math
COLD_TEMP_IN = 25  # °C
COLD_TEMP_OUT = 35  # °C
W_COLD_MASS_FLOW = 1  #we are going to calculate kg/s
Cp_COLD = 4.1784*1000  # J/(kg·°K)
COLD_DENSITY = 995.8175  # kg/m³
COLD_VISCOSITY = 803.4/(10**6)  # N·s/m^2
COLD_THERMAL_CONDUCTIVITY = 617.2/(10**3)  # W/(m·K)
COLD_PRANDTL_NUMBER = 5.452
#COLD_HEAT_TRANSFER_COEFF = 1  # we are going to calculate W/(m²·K)
#hot => mix
HOT_TEMP_IN = 80  # °C
HOT_TEMP_OUT = 50  # °C
W_HOT_MASS_FLOW = 13000/3600  # kg/s
Cp_HOT = 3.70592*1000 # J/(kg·°K)
HOT_DENSITY = 1008.8649 # kg/m³
HOT_VISCOSITY = 5.94967/(10**4)  # N·s/m^2
HOT_THERMAL_CONDUCTIVITY = 0.485208  # W/(m·K)
HOT_PRANDTL_NUMBER = 4.54 # unitless
#HOT_HEAT_TRANSFER_COEFF = 1  # we are going to calculate W/(m²·K)

tube_diameters = [(0.0159, 0.0004572), (0.0191, 0.0004572),
                  (0.0222, 0.000381),(0.0254, 0.0004064),
                  (0.0318, 0.000381),(0.0381, 0.000381),
                  (0.0508,0.0003302 )]  # (outer diameter, inner diameter) in meter
lengths_and_shell_diameter = [(2.438, 0.3715), (3.658, 0.6128),
                              (4.877, 0.7874),(6.096, 0.889)]  # (lenght, shell diameter) in metres
def max_number_of_tubes(tube_outer_diameter,shell_diameter):
    shell_area = math.pi*(shell_diameter**2)/4
    tube_area = math.pi*(tube_outer_diameter**2)/4
    max_number_of_tubes = (shell_area/tube_area)/0.907
    max_number_of_tubes = math.floor(max_number_of_tubes)
    return max_number_of_tubes

delta_T_1 = (HOT_TEMP_IN-COLD_TEMP_OUT) # celcius
delta_T_2 = (HOT_TEMP_OUT-COLD_TEMP_IN) # celcius

Q_owerall = W_HOT_MASS_FLOW*Cp_HOT*(HOT_TEMP_IN-HOT_TEMP_OUT) # celcius
W_COLD_MASS_FLOW =  Q_owerall/(Cp_COLD*(COLD_TEMP_OUT-COLD_TEMP_IN)) #celcius

F_for_more_then_one_pass = 0.98
LMDT = (delta_T_1-delta_T_2)/ math.log((delta_T_1/delta_T_2), math.e) #celcius

average_temp_cold = (COLD_TEMP_IN+COLD_TEMP_OUT)/2 #celcius
average_temp_hot = (HOT_TEMP_IN+HOT_TEMP_OUT)/2 #celcius


def area_for_cross_flow_5(tube_outer_diameter,shell_diameter):
    tube_pitch = tube_outer_diameter*1.25
    buffle_space = shell_diameter*0.6
    cross_flow_area = (tube_pitch-tube_outer_diameter)*buffle_space*(shell_diameter/tube_pitch)
    return cross_flow_area

def Gc_calculator_6(tube_outer_diameter,shell_diameter):
    cross_flow_area = area_for_cross_flow_5(tube_outer_diameter,shell_diameter)
    Gc = W_HOT_MASS_FLOW/cross_flow_area
    return Gc
def flow_area_per_tube_7(tube_inner_diameter):
    a_t = (math.pi*(tube_inner_diameter**2))/4
    return a_t

def area_for_baffle_window_8(tube_number,shell_diameter,tube_outer_diameter):
    F_b = 0.1955 # optionel
    N_b = tube_number/4 # must be integer
    battle_window_area = F_b*((math.pi*(shell_diameter**2))/4)-N_b*((math.pi*(tube_outer_diameter**2))/4)
    return battle_window_area
def Gb_calculator_9( tube_number,shell_diameter,tube_outer_diameter):
    battle_window_area = area_for_baffle_window_8(tube_number,shell_diameter,tube_outer_diameter)
    Gb = W_HOT_MASS_FLOW/battle_window_area
    return Gb
def flow_area_per_pass_10(tube_number,tube_inner_diameter):
    a_t = flow_area_per_tube_7(tube_inner_diameter)
    n_t = 1 # number of tube passs
    flow_area_per_pass = (tube_number*a_t)/n_t
    return flow_area_per_pass
def mass_velocity_Ge_calculator_11(tube_outer_diameter,shell_diameter,tube_number):
    Gc = Gc_calculator_6(tube_outer_diameter,shell_diameter)
    Gb= Gb_calculator_9( tube_number,shell_diameter,tube_outer_diameter)
    Ge = (Gc*Gb)**(1/2)
    return Ge
def average_velocity_calculator_12(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    A_t = flow_area_per_pass_10(tube_number,tube_inner_diameter)
    average_velocity = W_COLD_MASS_FLOW/(COLD_DENSITY*A_t)
    return average_velocity
def Reynolds_s_calculator_13(tube_outer_diameter):
    Ge = mass_velocity_Ge_calculator_11(tube_outer_diameter,shell_diameter,tube_number)
    Reynolds_s = (tube_outer_diameter*Ge)/HOT_VISCOSITY
    return Reynolds_s
def Reynolds_t_calculator_14(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    v = average_velocity_calculator_12(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    Reynolds_t = (COLD_DENSITY*v*tube_inner_diameter)/COLD_VISCOSITY
    return Reynolds_t
def Pr_s_calculator_15():
    Pr_s = (Cp_HOT*HOT_VISCOSITY)/HOT_THERMAL_CONDUCTIVITY
    return Pr_s
def Pr_t_calculator_16():
    Pr_t = (Cp_COLD*COLD_VISCOSITY)/COLD_THERMAL_CONDUCTIVITY
    return Pr_t
def fi_s_calculator_17():
    fi_s = (HOT_VISCOSITY/COLD_VISCOSITY)**(0.14)
    return fi_s
def Hs_calculator_with_donohue_18(tube_outer_diameter,shell_diameter,tube_number):
    Ge = mass_velocity_Ge_calculator_11(tube_outer_diameter,shell_diameter,tube_number)
    fi_s = fi_s_calculator_17()
    Hs = (0.2*(((tube_outer_diameter*Ge)/HOT_VISCOSITY)**(0.6))*((Cp_HOT*HOT_VISCOSITY)/HOT_THERMAL_CONDUCTIVITY)**(0.33))*((HOT_THERMAL_CONDUCTIVITY)/tube_outer_diameter*fi_s)
    return Hs
def Ho_corrected_coef_calculator_19(tube_outer_diameter,shell_diameter,tube_number):
    Ho = Hs_calculator_with_donohue_18(tube_outer_diameter,shell_diameter,tube_number)
    return Ho
def design_overall_coeff_calculator_21(tube_number,tube_lenght,tube_outer_diameter):
    area = math.pi*tube_lenght*tube_outer_diameter*tube_number
    U_d = Q_owerall/(area*LMDT)
    return U_d
def Gnielinski_A_calculator_23(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    #epsilone/D = zero => clean pipe
    Re_t = Reynolds_t_calculator_14(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    A = (7.149 / Re_t) ** (0.8981)
    if 2300<=Re_t<=5000000:
        if Re_t <= 10000:
            Re_d = Re_t -1000
            return A,Re_t, Re_d
        else:
            Re_d = Re_t
            return A, Re_t, Re_d
    else:
        return -1,-1,-1
def Gnielinski_f_calculator_24(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    A, Re_t, Re_d = Gnielinski_A_calculator_23(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    if A==-1:
        return -1,-1
    else:
        f = (1/(-4*math.log(((-5.0452/Re_t)*math.log(A,10)),10)))**2
        return f, Re_d
def Gnielinski_Nu_calculator_25(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    Pr_t = Pr_t_calculator_16()
    if 0.5 <= Pr_t <= 2000:
        f, Re_d = Gnielinski_f_calculator_24(tube_number, tube_inner_diameter, tube_outer_diameter, shell_diameter)
        if f != -1:
            Nu = ((f / 2) * Re_d * Pr_t) / (1 + 12.7 * ((f / 2) ** (1 / 2)) * (Pr_t ** (2 / 3)) - 1)
            return Nu
        else:
            return -1
    else:
        return -1
def Hi_calculator_with_donohue_26(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    Nu = Gnielinski_Nu_calculator_25(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    if Nu != -1:
        Hi = (Nu*COLD_THERMAL_CONDUCTIVITY)/tube_inner_diameter
        return Hi
    else:
        return -1
def H_io_calculator_27(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter):
    Hi = Hi_calculator_with_donohue_26(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    if Hi !=-1:
        h_io = Hi*(tube_inner_diameter/tube_outer_diameter)
        return h_io
    else:
        return -1
def clean_overall_coeff_calculator_20(tube_inner_diameter, tube_outer_diameter,shell_diameter,tube_number):
    ho = Ho_corrected_coef_calculator_19(tube_outer_diameter,shell_diameter,tube_number)
    h_io = H_io_calculator_27(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    if h_io != -1:
        U_c = (h_io*ho)/(h_io+ho)
        return U_c
    else:
        return -1
def dirt_factor_calculator_22(tube_number,tube_lenght,tube_outer_diameter,tube_inner_diameter,shell_diameter):
    U_d = design_overall_coeff_calculator_21(tube_number,tube_lenght,tube_outer_diameter)
    U_c = clean_overall_coeff_calculator_20(tube_inner_diameter, tube_outer_diameter,shell_diameter,tube_number)
    if U_c != -1:
        R_d = (1/U_c) - (1/U_d)
        return R_d
    else:
        return -1
#pressure drop calculations
def Reynolds_number_calculator_28(tube_outer_diameter,shell_diameter):
    tube_pitch = tube_outer_diameter*1.25
    D_h = 4*((tube_pitch/2)*0.86*tube_pitch-0.5*math.pi*(tube_outer_diameter**2)/4)/(0.5*math.pi*tube_outer_diameter)
    cross_flow_area = area_for_cross_flow_5(tube_outer_diameter,shell_diameter)
    Gc = Gc_calculator_6(tube_outer_diameter,shell_diameter)
    Re = (D_h*Gc)/HOT_VISCOSITY
    return Re ,D_h , Gc
def f_calculator_with_Re_29(tube_outer_diameter,shell_diameter):
    Re, D_h , Gc= Reynolds_number_calculator_28(tube_outer_diameter,shell_diameter)
    if Re >= 400:
        f = 1.7789*(Re**(-0.195868))
        return f , D_h , Gc
    else :
        return -1 , -1 , -1
def number_of_crosses_30(shell_diameter,lenght):
    Buffle_spacing = 0.6*shell_diameter
    N_plus_1 = lenght/Buffle_spacing
    return N_plus_1
def inside_tube_t_delta_P_calc_31(lenght,tube_inner_diameter,tube_number,tube_outer_diameter,shell_diameter):
    v = average_velocity_calculator_12(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    f , D_h , Gc = f_calculator_with_Re_29(tube_outer_diameter,shell_diameter)
    if f != -1:
        delta_P_t = ((4*f*lenght*tube_number)/tube_inner_diameter)*COLD_DENSITY*((v**2)/2)
        return delta_P_t
    else:
        return -1
def tube_side_r_return_delta_P_calc_32( tube_number, tube_inner_diameter,tube_outer_diameter, shell_diameter):
    v = average_velocity_calculator_12(tube_number,tube_inner_diameter,tube_outer_diameter,shell_diameter)
    delta_P_r = 4*tube_number*(COLD_DENSITY*(v**2)/2)
    return delta_P_r
def total_pressure_drop_calculator(lenght,tube_inner_diameter,tube_number,tube_outer_diameter,shell_diameter):
    delta_P_t = inside_tube_t_delta_P_calc_31(lenght,tube_inner_diameter,tube_number,tube_outer_diameter,shell_diameter)
    if delta_P_t == -1:
        return -1
    else:
        delta_P_r = tube_side_r_return_delta_P_calc_32( tube_number, tube_inner_diameter,tube_outer_diameter, shell_diameter)
        total_drop = delta_P_t + delta_P_r
        return total_drop
def delta_P_s_drop_calculator_34(tube_outer_diameter,shell_diameter,lenght):
    f ,D_h ,Gc = f_calculator_with_Re_29(tube_outer_diameter,shell_diameter)
    if f == 1 :
        return -1
    else:
        N_plus_one = number_of_crosses_30(shell_diameter,lenght)
        fi_s =  fi_s_calculator_17()
        delta_P_s_drop = (f*(Gc**2)*shell_diameter*N_plus_one)/(2*D_h*fi_s*HOT_DENSITY)
        return delta_P_s_drop
combinations = []
for length, shell_diameter in lengths_and_shell_diameter:
    for outer_diameter, thickness in tube_diameters:
        inner_diameter = outer_diameter - thickness
        max_of_tubes = max_number_of_tubes(outer_diameter,shell_diameter)
        combinations.append((outer_diameter, inner_diameter, length, shell_diameter,max_of_tubes))
results = []
for outer_diameter, inner_diameter, length, shell_diameter,max_of_tubes in combinations:
    for i in range(1,max_of_tubes+1):
        if i % 4 == 0:
            dirt_factor = dirt_factor_calculator_22(i, length, outer_diameter, inner_diameter, shell_diameter)
            if dirt_factor != -1:
                delta_P_s = delta_P_s_drop_calculator_34(outer_diameter, shell_diameter, length)
                total_P_drop = total_pressure_drop_calculator(length, inner_diameter, i, outer_diameter, shell_diameter)
                if total_P_drop < 172368.925:
                    results.append((outer_diameter, shell_diameter, i, length, total_P_drop, delta_P_s))
            else:
                pass
results.sort(key=lambda x: x[4])
for outer_diameter, shell_diameter, i, length, total_P_drop, delta_P_s in results:
   print("Outer dia: {} m , Shell dia: {} m , Total Tube Number: {}, Lenght: {} m , delta P total: {} Pa, Delta Ps: {} Pa".format(outer_diameter, shell_diameter, i, length, total_P_drop, delta_P_s))
print(len(results))
