import numpy as np

def sfactor(dPML, N, omega_0, m = 4):
    lnR = -16;  ## final desired reflection coefficient

    PmlBound = [dPML, N - dPML];
    sigma_max = (m + 1) * lnR / 2
    # need start and stop coordinates of the PML
    # abs(wrange - wpml);
    sfactor_domain = np.ones((N - 2 * dPML))
    d = np.linspace(0, dPML, dPML)
    pmlPolynomial = 1/(1 + (1j) * (sigma_max * d / omega_0) ** m);
    sfactor_domain = np.concatenate((np.flip(pmlPolynomial, 0),\
                                     sfactor_domain, pmlPolynomial), axis=0);

    sfactorMatrix = np.repeat(sfactor_domain, N);
    sfactorMatrix = np.reshape(sfactorMatrix, [N,N])
    return sfactorMatrix;
