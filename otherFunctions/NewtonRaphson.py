import numpy as np

def newtonRaphson_poly(a:float, b:float, c:float, n:float, x0:float, tol=1e-6, max_iter=100):

    """
    Function solves the polynomial equation ax^n + bx + c = 0 using the Newton-Raphson method.

    Parameters:
    a: Coefficient of x^n
    b: Coefficient of x
    c: Constant term
    n: Degree of the polynomial
    x0: Initial guess for the root
    tol: Tolerance for convergence (default: 1e-6)
    max_iter: Maximum number of iterations (default: 100)

    Returns:
    float: The root of the polynomial (if found)
    """
    # Define the polynomial and its derivative
    def poly(x):
        return a * x**n + b * x + c

    def derivative(x):
        return n * a * x**(n - 1) + b

    # Newton-Raphson iteration
    x = x0
    for i in range(max_iter):
        f_x = poly(x)
        f_prime_x = derivative(x)

        if abs(f_x) < tol:
            #print(f"Converged in {i} iterations.")
            return x

        if f_prime_x == 0:
            #raise ValueError("Derivative is zero. Newton-Raphson method fails.")
            return x

        # Update the guess
        x -= f_x / f_prime_x

    #raise ValueError("Newton-Raphson method did not converge.")
    #in case the solution is not converged, take the latest the latest x
    return x

# Example Usage
'''
a, b, c, n = 1, -3, 2, 2  # Coefficients for x^2 - 3x + 2 = 0
x0 = 1.3  # Initial guess
root = newton_raphson_poly(a, b, c, n, x0)
print(f"Root: {root}")
'''