import numpy as np
import scipy.special as sp


def compute_phi(theta, q):
    """
    Compute phi_k values for k = -q to q.
    """
    phi = np.zeros(q + 1)  # Stores phi_k for k = -q to q

    for k in range(0, q + 1):
        phi[k] = sum(theta[s] * theta[s - abs(k)]
                     for s in range(abs(k), q + 1))

    return phi


def compute_gamma_j(theta, d, sigma_sq, j):
    """
    Compute gamma_j using the given theta coefficients and variance sigma_sq.
    """
    q = len(theta) - 1  # q is determined by theta size (q+1 elements in theta)


    phi = compute_phi(theta, q)  # Compute phi_k values
    print(f"phi: {phi}")
    gamma_j = np.zeros(j)
    print(f"gamma_j: {gamma_j}")
    for j in range(0, j):
        for k in range(-q, q+1):
            print(f"j: {j}, k: {k}")
            gamma_j[j] = gamma_j[j] + sigma_sq * phi[abs(k)] \
                * (sp.gamma(1 - 2 * d) / sp.gamma(1 - d)**2) \
                * (pochhammer(d, k - j) / pochhammer(1 - d, k - j))
    print(f"gamma_j: {gamma_j}")
    return gamma_j


def pochhammer(d, i):
    if (i == 0):
        return 1
    elif (i < 0):
        return pochhammer(d, i + 1) / (d + i)
    else:
        return (d + i - 1) * pochhammer(d, i - 1)


if __name__ == "__main__":
    # Example usage
    # Example theta values for q=2
    theta = np.array([1, 0.8330467201, 0.2358215965, 0.2533662929])
    sigma_sq = 0.8682117  # Example variance
    d = -0.175258
    j = 49  # Compute gamma_1

    gamma_j = compute_gamma_j(theta, d, sigma_sq, j)
