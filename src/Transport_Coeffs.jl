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

function calc_bnry_Dij_mat(T, p, comp_vec)
    """Calculates the binary diffusion coefficient matrix with the Fuller equation.

    Args:
        T (float): Temperature [K]
        p (float): Pressure [Pa]
        comp_vec (list): Components in the system (e.g., ["N2", "H2", "CO"])

    Returns:
        Dij_mat (array): Matrix with the binary diffusion coefficients for the components [m2 / s]
    """
    
    Ni = length(comp_vec)
    Dij_mat = zeros((Ni, Ni))
    for (i, compi) in enumerate(comp_vec)
            for (j, compj) in enumerate(comp_vec)
                Dij_mat[i,j] = calc_D12(T, p, compi, compj)
                Dij_mat[j,i] = Dij_mat[i,j]
            end
        end
    return Dij_mat
end

function calc_Dmixi_vec(T, p, comp_vec, yi_vec)
    """Calculates the mixture-averaged diffusion coefficient from binary diffusion coefficients from the Chapman-Enskog equation.

    Args:
        T (float): Temperature [K]
        p (float): Pressure [Pa]
        comp_vec (list): Components in the system (e.g., ["N2", "H2", "CO"])
        yi_vec (list): Mole fractions of the components in the system (in the same order as comp_vec)

    Returns:
        Dmixi_vec (array): Mixture diffusion coefficients for the components [m2 / s]
    """
    
    Dmixi_vec = zeros(length(comp_vec))
    for (i, compi) in enumerate(comp_vec)
            for (j, compj) in enumerate(comp_vec)
                if compj == compi
                     continue
                else
                    Dij = calc_D12(T, p, compi, compj)
                    Dmixi_vec[i] = Dmixi_vec[i] + yi_vec[j] * Dij
                end
            end
            Dmixi_vec[i] *= 1 / (1-yi_vec[i])
        end
    return Dmixi_vec
end

function calc_Deffi_vec(T, p, comp_vec, yi_vec, Morph_par)
    """Calculates the effective diffusion coefficient from mixture averaged and Knudsen diffusion coefficients.

    Args:
        T (float): Temperature [K]
        p (float): Pressure [Pa]
        comp_vec (list): Components in the system (e.g., ["N2", "H2", "CO"])
        yi_vec (list): Mole fractions of the components in the system (in the same order as comp_vec)
        Morph_par (dict): Morphological parameters of the solid medium where diffusion takes place 
                          (e.g., Dict("d_pore" => 25e-9, "eps" => 0.5, "tau" => 4.0))

    Returns:
        Deffi_vec (array): Effective diffusion coefficients for the components [m2 / s]
    """

    Dmixi_vec = calc_Dmixi_vec(T, p, comp_vec, yi_vec)
    Deffi_vec = zeros(length(comp_vec))

    for (i, comp) in enumerate(comp_vec)
        Mi = fuller_df[findfirst(==(comp), fuller_df.Component), :]["M"]
        DKi = 97/2 * Morph_par["d_pore"] * sqrt.(T / Mi)
        Deffi_vec[i] = 1.0 / (1.0/Dmixi_vec[i] + 1.0/DKi) * Morph_par["eps"]/Morph_par["tau"]
    end

    return Deffi_vec
end

export calc_D12