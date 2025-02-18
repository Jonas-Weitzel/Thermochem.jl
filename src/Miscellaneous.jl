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

function calc_PFR_thermodat(Cpcoef_ik_mat, dHf0i_vec, Tz_vec)
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

    Nz = length(Tz_vec)
    tz_vec = Tz_vec/1000

    Cpiz_mat = Cpcoef_ik_mat[1:end, 1:5] * (@. [ones(Nz) tz_vec tz_vec^2 tz_vec^3 tz_vec^(-2)])'         
    dHfiz_mat = dHf0i_vec .+ Cpcoef_ik_mat * (@. [tz_vec (tz_vec^2 / 2) (tz_vec^3 / 3) (tz_vec^4 / 4) -tz_vec ones(Nz) zeros(Nz) -ones(Nz)])'

    return Cpiz_mat, dHfiz_mat
end