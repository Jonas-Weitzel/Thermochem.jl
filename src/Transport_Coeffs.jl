function calc_D12(T, p, comp1, comp2; fuller_df=fuller_df)
    """Calculates the binary diffusion coefficient from the Fuller equation.

    Args:
        T (float): Temperature [K]
        p (float): Pressure [Pa]
        comp1 (str): Formula of component 1
        comp2 (str): Formula of component 2

    Returns:
        D12 (float): Binary diffusion coefficient [m2 / s].
    """

    comp1_dat = fuller_df[findfirst(==(comp1), fuller_df.Component), :]
    comp2_dat = fuller_df[findfirst(==(comp2), fuller_df.Component), :]
    D12 = 1e-3 * ( T^1.75 * ( 1/comp1_dat["M"] + 1/comp2_dat["M"] )^0.5 ) / 
                ((p/101325) * (comp1_dat["V"]^(1.0/3.0) + comp2_dat["V"]^(1.0/3.0))^2)      # [cm2 / s]
    D12 = D12 * (1e-4)                       # [m2 / s]

end