function calc_comp_thermodat(comp, T; gen_dat=gen_dat, Cp_dat=Cp_dat)
    """Function to compute the isobaric heat capacity, and the enthalpy of formation for a component. Correlation
    and data from the NIST chemistry webbook.

    Args:
        comp (str): Name of component.
        T (float): Temperature at which to compute the thermo-chemistry data [K].

    Returns:
        Cp (float): Molar heat capacity [J / mol K].
        dHf (float): Molar enthalpy of formation [J / mol].
    """

    A,B,C,D,E,F,G,H = (values(Cp_dat[findfirst(==(comp), Cp_dat.Component), 2:end]))
    t = T/1000                                                                                       # [K]
    Cp = (A + B*t + C*(t^2) + D*(t^3) + E/(t^2))                                                     # [J / mol K]
    dHf0 = gen_dat[findfirst(==(comp), gen_dat.Component), 2:end]["dHf, J / mol"]
    dHf = dHf0 + (A*t + (B*t^2)/2 + (C*t^3)/3 + (D*t^4)/4 - E/t + F - H) * 1000                      # [J / mol]
    
    return Cp, dHf
end

function get_thermo_coeffs(comp_vec; gen_dat=gen_dat, Cp_dat=Cp_dat)
    """Unpacks the relevant coefficients for thermochemical calculations for the components in `comp_vec`.

    Args:
        comp_vec (array): Names of the components in the system.

    Returns:
        Cpcoef_ik_mat (array): Temperature coefficients of the components for calculating the molar heat capacity [varies].
        dHf0i_vec (array): Standard heat of formation for the components [J / mol].
    """

    Ni = length(comp_vec)
    Cpcoef_ik_mat = zeros(Ni, 8)
    dHf0i_vec = zeros(Ni, 1)

    for (i, comp) in enumerate(comp_vec)
        Cpcoef_ik_mat[i, 1:end] = Array(Cp_dat[findfirst(==(comp), Cp_dat.Component), 2:end])
        dHf0i_vec[i] = gen_dat[findfirst(==(comp), gen_dat.Component), 2:end]["dHf, J / mol"]
    end
    
    return Cpcoef_ik_mat, dHf0i_vec
end

function calc_mix_thermodat!(Cpi_vec, dHfi_vec, Cpcoef_ik_mat, dHf0i_vec, T)
    """Calculates the molar heat capacity and heat of formation for the species i along the reactor coordinate z of a PFR.

    Args:
        Cpiz_mat (array): Storage for result.
        dHfiz_mat (array): Storage for result.
        Cpcoef_ik_mat (array): Temperature coefficients of the components for calculating the molar heat capacity [varies].
        dHf0i_vec (array): Standard heat of formation for the components [J / mol].
        Tz_vec (array): Temperature along the reactor coordinate [K].

    Returns:
        Cpiz_mat (array): Molar heat capacity [J / mol K].
        dHfiz_mat (array): Molar enthalpy of formation [J / mol]
    """

    t = T/1000
    t2 = t^2; t3 = t^3; t4 = t^4
    Ni = length(dHf0i_vec)

    for i in 1:Ni
        A,B,C,D,E,F,G,H = Cpcoef_ik_mat[i,:]
        Cpi_vec[i] = A + B*t + C*t2 + D*t3 + E/t2       
        dHfi_vec[i] = dHf0i_vec[i] + A*t + B*t2/2 + C*t3/3 + C*t4/4 - E/t + F - H
    end

    return nothing
end