import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch

def calculate_factor_loading(
    input_df: pd.DataFrame,
    factors: list[str],
    assets: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, np.array, np.array]: 
    """
    Given DataFrame of (non-excess) asset returns 
    and factor returns, 
    returns panel data of factor loadings

    Args:
        input_df (pd.DataFrame): DataFrame indexed on date,
            with column names corresponding to assets
        factors (list[str]): list of factors
        assets (list[str]): list of risky assets

    Returns:
        pd.DataFrame: panel data of factor loadings
        pd.DataFrame: modified returns dataframe with excess returns
        np.array: realized covariance matrices of shape (N, K, T)
    """

    assert type(input_df.index) == pd.DatetimeIndex, "input_df has wrong index"
    for factor in factors:
        assert factor in input_df.columns, f"missing factor {factor}"
    for asset in assets:
        assert asset in input_df.columns, f"missing asset {asset}"
    assert "RF" in input_df.columns, f"Missing risk free"
    assert "Quarter" in input_df.columns, f"Missing quarter"
    
    input_df.sort_index(inplace=True)
    N = len(assets)
    K = len(factors)
    T = input_df['Quarter'].nunique()

    for col in assets:
        input_df.loc[:, col] = input_df[col] - input_df["RF"]
    
    cols = list(assets) + list(factors)

    realized_covariance_matrices = np.zeros((N, K, T))
    residuals = np.zeros((N, K, T))

    quarters = sorted(input_df['Quarter'].unique())
    for i, quarter in enumerate(quarters):
        returns = (
            input_df.loc[
                input_df['Quarter'] == quarter,
                cols
            ]
            .values
        )
        Omega_hat_t = returns.T @ returns
        if (i == T):
            print(f"DEBUG ERROR FOUND")
        realized_covariance_matrices[:, :,i] = Omega_hat_t[:N, N:N+K]
    
    beta_loading = pd.DataFrame(
        index = pd.MultiIndex.from_product([assets, factors]),
        columns = input_df['Quarter'].unique(),
    )

    for i, asset in enumerate(assets):
        for j, factor in enumerate(factors):
            omega_i_j_series = realized_covariance_matrices[i, j, :]
            Y = omega_i_j_series[1:]
            X = (
                np.column_stack([
                    np.ones(len(Y)),
                    omega_i_j_series[:-1]
                ])
            )
            b = np.linalg.lstsq(X, Y, rcond=None)[0]
            delta0, delta1 = b
            beta_loading.loc[(asset, factor)] = delta0 + delta1 * omega_i_j_series

            fitted = delta0 + delta1 * omega_i_j_series[:-1]
            resid = omega_i_j_series[:-1] - fitted
            residuals[i, j, 1:] = resid
            # residuals[i, j, 0] = np.nan


    return beta_loading, input_df, realized_covariance_matrices, residuals

def penalty_based_minimization(
    beta_hat: np.array,
    excess_returns: np.array,
    N: int,
    K: int,
    R: int, 
    T: int,
    lam: float = 1.0,
    lr: float= 1e-2,
    n_iter: int = 2000,
    device: str = "cpu",
    seed: int = 0
) -> tuple[np.array, np.array, np.array, np.array]:
    """
    Solves unconstrained version of equation (24) 
    With penalty

    Args:
        beta_hat (np.array): estimated beta loadings, N * K * T
        excess_returns (np.array): excess returns, N * T
        N (int): number of assets
        K (int): number of observed assets
        R (int): number of unobserved factors
        T (int): number of time periods
        lam (float): penalty weight on deviation from identity
        lr (float): learning rate
        n_iter (int): number of iterations
        device (str): cpu
        seed (int): for reproducibility

    Returns:
        tuple[np.array, np.array, np.array]:
            eta: N * (1 + K)
            G: T * R
            beta^*: N * R
            objective: np.array of dimensions num_iter
    """

    assert beta_hat.ndim == 3, f"beta_hat must be 3D, got {beta_hat.ndim}"
    assert beta_hat.shape == (N, K, T), f"beta_hat.shape {beta_hat.shape} != ({N}, {K}, {T})"
    assert excess_returns.ndim == 2, f"excess_returns must be 2D, got {excess_returns.ndim}"
    assert excess_returns.shape == (N, T), f"excess_returns.shape {excess_returns.shape} != ({N}, {T})"

    torch.manual_seed(seed)

    beta_hat_t = torch.from_numpy(beta_hat).float().to(device)
    beta_hat_t = beta_hat_t.permute(0, 2, 1)  #(N, T, K)
    r = torch.from_numpy(excess_returns).float().to(device) #(N, T)

    ones = torch.ones((N, T, 1), device=device)
    X = torch.cat([ones, beta_hat_t], dim=2) # (N, T, 1 + K)

    # parameters of optimization problem 
    eta = torch.nn.Parameter(torch.zeros(N, 1 + K, device=device))
    beta_star = torch.nn.Parameter(torch.zeros(N, R, device=device))
    G = torch.nn.Parameter(torch.zeros(T, R, device=device))

    torch.nn.init.normal_(eta, mean=0.0, std=0.1)
    torch.nn.init.normal_(beta_star, mean=0.0, std=0.1)
    torch.nn.init.normal_(G, mean=0.0, std=0.1)

    optimizer = torch.optim.Adam([eta, beta_star, G], lr=lr)
    I_R = torch.eye(R, device=device)

    objective = np.empty(shape=(n_iter))

    for it in range(n_iter):
        optimizer.zero_grad()
        obs_part = (X * eta[:, None, :]).sum(dim=2)
        latent_part = (G @ beta_star.t()).t()
        pred = obs_part + latent_part
        mse_loss = torch.mean((r - pred) ** 2)
        GTG = G.t() @ G / T 
        penalty = torch.norm(GTG - I_R, p='fro')**2 
        loss = mse_loss + lam * penalty
        loss.backward()
        optimizer.step()

        objective[it] = mse_loss.item()

        if (it + 1) % 100 == 0:
            log_str = (
                f"Iter {it + 1}/{n_iter}, "
                f"objective={mse_loss.item():.6f}, "
                f"loss={loss.item():.6f}, "
                f"pen={penalty.item():.6f}"
            )
            print(log_str)
    
    eta_np = eta.detach().cpu().numpy()
    G_np = G.detach().cpu().numpy()
    beta_star_np = beta_star.detach().cpu().numpy()

    return eta_np, G_np, beta_star_np, objective

def iterative_convergence(
    beta_hat: np.array,
    excess_returns: np.array,
    N: int,
    K: int,
    R: int, 
    T: int,
    rtol: float = 1e-05,
    atol: float = 1e-08,
    n_iter: int = 2000,
    seed: int = 0,
    verbose: bool = False,
) -> tuple[np.array, np.array, np.array]:
    """
    Solves constrained optimization by iterating
    until convergence

    Args:
        beta_hat (np.array): estimated beta loadings, N * K * T
        excess_returns (np.array): excess returns, N * T
        N (int): number of assets
        K (int): number of observed factors
        R (int): number of unobserved factors
        T (int): number of time periods
        rtol (float): relative tolerance for convergence, refer to
            numpy.allclose documentation
        atol (float): absolute tolerance for convgernce
        n_iter (int): number of iterations
        seed (int): for reproducibility

    Returns:
        tuple[np.array, np.array, np.array, np.array]:
            eta: N * (1 + K)
            G: T * R
            beta^*: N * R
            objective: n_iter
    """

    assert beta_hat.ndim == 3, f"beta_hat must be 3D, got {beta_hat.ndim}"
    assert beta_hat.shape == (N, K, T), f"beta_hat.shape {beta_hat.shape} != ({N}, {K}, {T})"
    assert excess_returns.ndim == 2, f"excess_returns must be 2D, got {excess_returns.ndim}"
    assert excess_returns.shape == (N, T), f"excess_returns.shape {excess_returns.shape} != ({N}, {T})"

    np.random.seed(seed)

    beta_hat_t = beta_hat.transpose(0, 2, 1)  #(N, T, K)
    r = excess_returns                      #(N, T)
    ones = np.ones((N, T, 1))
    X = np.concatenate([ones, beta_hat_t], axis=2) # (N, T, 1 + K)

    eta = np.random.normal(0, 0.1, size = (N, 1 + K))
    beta_star = np.random.normal(0, 0.1, size = (N, R))
    G = np.random.normal(0, 0.1, size = (T, R))

    objective = np.empty(shape=(n_iter))

    for i in range(N):
        Xi = X[i]              # (T, 1+K)
        ri = r[i]              # (T,)

        # (X'X)^{-1} X'r
        A = Xi.T @ Xi          # (1+K, 1+K)
        b = Xi.T @ ri          # (1+K,)

        A = A + 1e-8 * np.eye(1 + K)
        eta[i] = np.linalg.solve(A, b)

    # U_i = r_i - X_i η_i
    U = np.zeros((T, N))
    for i in range(N):
        Xi = X[i]               # (T, 1+K)
        ri = r[i]               # (T,)
        ui = ri - Xi @ eta[i]   # (T,)
        U[:, i] = ui

    # (1/NT) \sum_{i=1}^{N} (r_i - X eta) (r - X eta)'
    S = (U @ U.T) / (N * T)      # (T, T)
    eigvals, eigvecs = np.linalg.eigh(S)
    G = eigvecs[:, -R:] * np.sqrt(T)            # (T, R)
    I_T = np.eye(T)

    for it in range(n_iter):
        G_old = G.copy()

        # M_G = I - G(G'G)^{-1}G'
        GtG = G.T @ G                    # (R, R)
        GtG_inv = np.linalg.inv(GtG)  # (R, R)
        Proj_G = G @ GtG_inv @ G.T       # (T, T)
        M_G = I_T - Proj_G

        # update eta and beta_star given G
        for i in range(N):
            Xi = X[i]                # (T, 1+K)
            ri = r[i]                # (T,)

            # eta = (X' M_G X)^{-1} X' M_G r
            A = Xi.T @ M_G @ Xi      # (1+K, 1+K)
            b = Xi.T @ M_G @ ri      # (1+K,)

            A = A + 1e-8 * np.eye(1 + K)
            eta[i] = np.linalg.solve(A, b)

            # calculate new residuals and beta_star
            vi = ri - Xi @ eta[i]              # (T,)
            beta_star[i] = GtG_inv @ (G.T @ vi)  # (R,)

        # update G
        U = np.zeros((T, N))
        for i in range(N):
            Xi = X[i]                  # (T, 1+K)
            ri = r[i]                  # (T,)
            ui = ri - Xi @ eta[i]
            U[:, i] = ui

        S = (U @ U.T) / (N * T)
        eigvals, eigvecs = np.linalg.eigh(S)
        G = eigvecs[:, -R:] * np.sqrt(T)

        # objective value
        U = np.zeros((T, N))
        for i in range(N):
            Xi = X[i]
            ri = r[i]
            U[:, i] = ri - Xi @ eta[i]  

        # term_i[t] = u_i[t] - G[t] @ beta_star[i]
        loss_matrix = np.zeros((T, N))
        for i in range(N):
            loss_matrix[:, i] = U[:, i] - G @ beta_star[i]

        obj_value = np.mean(loss_matrix ** 2)

        objective[it] = obj_value.item()

        if (it + 1) % 100 == 0 and verbose:
            loss = np.linalg.norm(G - G_old, ord="fro")
            log_str = (
                f"Iter {it + 1}/{n_iter}, "
                f"frobenius_norm(G - G_old)={loss:.4f}, "
                f"objective={obj_value:.6f}"
            )
            print(log_str)

        if np.allclose(G, G_old, rtol=rtol, atol=atol):
            if verbose:
                print(f"Converged at iteration {it+1}")
            break

    return eta, G, beta_star, objective

def estimate_avar(
    beta_hat: np.array,
    excess_returns: np.array,
    eta: np.array,
    G: np.array,
    beta_star: np.array,
    realized_covariance: np.array,
    residuals: np.array,
    N: int,
    K: int,
    R: int, 
    T: int,
) -> np.array:
    """
    Estimates Avar matrix using equation (30)

    Args:
        beta_hat (np.array): estimated time-varying betas, N * K * T
        excess_returns (np.array): excess returns, N * T
        eta (np.array): estimated eta,  N * (1 + K)
        G (np.array): estimated G matrix, T * R
        beta_star (np.array): estimated beta^*, N * R
        realized_covariance (np.array): realized covariance matrix, N * K * T
        residuals (np.array): residuals from AR(1) regression, N * K * T
        N (int): number of assets
        K (int): number of observed factors
        R (int): number of unobserved factors
        T (int): number of time periods

    Returns:
        np.array: estimated asymptotic variance N(K + 1) * N(K + 1)
    """

    assert beta_hat.ndim == 3, f"beta_hat must be 3D, got {beta_hat.ndim}"
    assert beta_hat.shape == (N, K, T), f"beta_hat.shape {beta_hat.shape} != ({N}, {K}, {T})"
    assert excess_returns.ndim == 2, f"excess_returns must be 2D, got {excess_returns.ndim}"
    assert excess_returns.shape == (N, T), f"excess_returns.shape {excess_returns.shape} != ({N}, {T})"
    assert eta.ndim == 2, f"eta must be 2D, got {eta.ndim}"
    assert eta.shape == (N, (1 + K)), f"eta.shape {eta.shape} != ({N}, {1 + K})"
    assert G.ndim == 2, f"G must be 2D, got {G.ndim}"
    assert G.shape == (T, R), f"G.shape {G.shape} != ({T}, {R})"
    assert beta_star.ndim == 2, f"beta_star must be 2D, got {beta_star.ndim}"
    assert beta_star.shape == (N, R), f"beta_star.shape {beta_star.shape} != ({N}, {R})"
    assert realized_covariance.ndim == 3, f"realized_covariance must be 3D, got {realized_covariance.ndim}"
    assert realized_covariance.shape == (N, K, T), f"realized_covariance.shape {realized_covariance.shape} != ({T}, {R})"
    assert residuals.ndim == 3, f"residuals must be 3D, got {residuals.ndim}"
    assert residuals.shape == (N, K, T), f"residuals.shape {residuals.shape} != ({T}, {R})"

    beta_hat_t = beta_hat.transpose(0, 2, 1)           # (N, T, K)
    ones = np.ones((N, T, 1))
    X = np.concatenate([ones, beta_hat_t], axis=2)     # (N, T, K+1)
    r = excess_returns                                 # (N, T)

    # projection matrix G
    I_T = np.eye(T)                                    # (T, T)
    GtG = G.T @ G                                      # (R, R)
    M_G = I_T - G @ np.linalg.inv(GtG) @ G.T           # (T, T)

    # Build block diagonal S
    Kp1 = K + 1
    S_hat = np.zeros((N * Kp1, N * Kp1))                # (N(K+1), N(K+1))

    for i in range(N):
        Xi = X[i]                                       # (T, K+1)
        Sii = Xi.T @ M_G @ Xi * (1.0/T)                       # (K+1, K+1)
        r0 = i*Kp1
        S_hat[r0:r0+Kp1, r0:r0+Kp1] = Sii

    # Build L
    L_hat = np.zeros((N * Kp1, N * Kp1))                # (N(K+1), N(K+1))
    GtG_over_N_inv = np.linalg.inv(GtG / N)

    for i in range(N):
        Xi = X[i]
        for j in range(N):
            Xj = X[j]

            a_ij = beta_star[i].T @ GtG_over_N_inv @ beta_star[j]
            Lij = (Xi.T @ M_G @ Xj) * (a_ij/T)
            r0 = i*Kp1
            c0 = j*Kp1
            L_hat[r0:r0+Kp1, c0:c0+Kp1] = Lij
    
    # sigma2_hat
    eps_sum = 0
    for i in range(N):
        for t in range(T):
            pred = X[i,t] @ eta[i] + G[t] @ beta_star[i]
            eps_sum += (r[i,t] - pred)**2

    df = N*T - N*Kp1 - (N+T)*R
    sigma2_hat = eps_sum / df

    W = np.zeros((N * Kp1, N * Kp1))
    for i in range(N):
        H_i = np.zeros((T, K))
        for j in range(K):
            Z_ij = np.column_stack([np.ones(T), realized_covariance[i, j, :]])
            ZTZ_over_T = (Z_ij.T @ Z_ij) * (1.0/T)
            v_ij = residuals[i, j, :]
            h_ij = Z_ij @ np.linalg.inv(ZTZ_over_T) @ (Z_ij.T @ v_ij) / np.sqrt(T)
            H_i[:, j] = h_ij

        W_i = (Xi.T @ M_G @ Xi) * (1.0 / T) * sigma2_hat

        lambda_i = eta[i][1:]
        local_vec = lambda_i * (M_G @ Xi)[:, 1:] * (1.0 / T)

        diag_HTH_over_T = np.diag(np.diag(
            H_i.T @ H_i
        )) * (1.0 / T)

        for j in range(K):
            weight = diag_HTH_over_T[j, j]
            col = local_vec[:, j]
            W_i[1 + j, 1 + j] += weight * (col @ col)

        W[
            i * (K + 1): (i + 1)* (K + 1),
            i * (K + 1): (i + 1)* (K + 1)
        ] = W_i

    avar = clean((
        np.linalg.inv(
            S_hat - L_hat.T / N
        ) 
        @ W 
        @ np.linalg.inv(
            S_hat - L_hat / N
        )
    ))
    avar 
    return avar

def full_homogeneity_test(
    eta: np.array,
    avar: np.array,
    N: int, 
    K: int,
    T: int
) -> float:
    """
    Joint hypothesis testing, under H_0
        - all intercepts are 0 AND 
        - all slope coefficients equal across assets

    Returns:
        float: gamma_ad, asymptotically standard normal
    """
    p = N * (K + 1)
    assert avar.shape == (p, p)
    assert eta.shape == (N, (K + 1))

    eta_mean = eta.mean(axis=0)
    eta_centered = eta - eta_mean
    d_vec = eta_centered.reshape(p)

    avar_inv = np.linalg.inv(avar)
    avar_inv = clean(avar_inv)

    W = T * d_vec.T @ avar_inv @ d_vec
    q = (N - 1) * (K + 1)

    gamma_ad = (W - q) / np.sqrt(2 * q)

    return gamma_ad


def intercept_homogeneity_test(
    eta: np.array,
    avar: np.array,
    N: int, 
    K: int,
    T: int
) -> float:
    """
    Under H_0
        - all intercepts are 0 AND 

    Returns:
        float: gamma_a, asymptotically standard normal
    """
    p = N * (K + 1)
    assert avar.shape == (p, p)
    assert eta.shape == (N, (K + 1))

    alpha = eta[:, 0]

    idx = np.arange(0, p, K + 1)


    avar_inv = np.linalg.inv(avar)
    avar_inv = clean(avar_inv)

    V_alpha = avar_inv[np.ix_(idx, idx)]


    W = alpha.T @ V_alpha @ alpha

    gamma_a = (W - (N-1)*K) / np.sqrt(2 * (N-1)*K)

    return gamma_a

def slope_homogeneity_test(
    eta: np.array,
    avar: np.array,
    N: int, 
    K: int,
    T: int
) -> float:
    """
    Under H_0
        - all slopes are equal

    Returns:
        float: gamma_a, asymptotically standard normal
    """
    p = N * (K + 1)
    assert avar.shape == (p, p)
    assert eta.shape == (N, (K + 1))

    slopes = eta[:,1:]
    slopes_centered = slopes - slopes.mean(axis = 0)
    slopes_vec = slopes_centered.reshape(N * K)

    avar_inv = np.linalg.inv(avar)
    avar_inv = clean(avar_inv)
    idx = []
    for i in range(N):
        for j in range(K):
            idx.append(i * (K + 1) + 1 + j)
    idx = np.array(idx)

    V_lambda = avar_inv[np.ix_(idx, idx)]

    W = slopes_vec.T @ V_lambda @ slopes_vec
    
    q = (N - 1) * K

    gamma_lambda = (W - q) / np.sqrt(2 * q)

    return gamma_lambda

def clean(x: np.array) -> np.array:
    return np.clip(
        x,
        a_min = np.percentile(x, 5),
        a_max = np.percentile(x, 95),
    )

