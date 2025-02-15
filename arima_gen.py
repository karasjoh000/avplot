import numpy as np
import pandas as pd


def gen_csv(args):
    # Parameters
    b0, b1, b2 = args.beta0, args.beta1, args.beta2
    phi = 0.6
    theta1, theta2, theta3 = 0.4, -0.2, 0.1
    n = args.n  # Number of data points

    # Generate exogenous variables
    np.random.seed(42)  # For reproducibility
    x1 = np.random.normal(0, args.x1_std, n)
    x2 = np.random.normal(0, args.x2_std, n)

    # Generate ARIMA(1, 0, 3) process
    # Extra points for lagged errors
    epsilon = np.random.normal(0, args.epsilon_std, n + 3)
    m = np.zeros(n + 3)

    for t in range(3, n + 3):
        m[t] = phi * m[t - 1] + epsilon[t] + theta1 * epsilon[t - 1] + \
            theta2 * epsilon[t - 2] + theta3 * epsilon[t - 3]

    m = m[3:]  # Discard the first 3 values (initialization)

    # Generate y_t
    y = b0 + x1 * b1 + x2 * b2 + m

    # Create a DataFrame to store the data
    data = pd.DataFrame({
        't': np.arange(1, n + 1),
        'x1': x1,
        'x2': x2,
        'm_t': m,
        'y_t': y
    })

    # save data to data dir
    import os
    path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(path, 'data')
    path = os.path.join(path, f'ma{args.ma}ar{args.ar}e{args.epsilon_std}n{n}.csv')

    # Save the data to a CSV file
    data.to_csv(path, index=False)

    # Display the first 10 rows of the data
    print(data.head(10))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate ARIMA data and save to CSV.")
    parser.add_argument("--ma", type=int, default=3, help="MA order")
    parser.add_argument("--ar", type=int, default=1, help="AR order")
    parser.add_argument("--beta0", type=float, default=2, help="Intercept")
    parser.add_argument("--beta1", type=float, default=1.5,
                        help="Coefficient for x1")
    parser.add_argument("--beta2", type=float, default=-
                        0.8, help="Coefficient for x2")
    parser.add_argument("--x1_std", type=float, default=1,
                        help="Standard deviation for x1")
    parser.add_argument("--x2_std", type=float, default=1,
                        help="Standard deviation for x2")
    parser.add_argument("--epsilon_std", type=float, default=1,
                        help="Standard deviation for epsilon")
    parser.add_argument("--n", type=int, default=50,
                        help="Number of data points")
    args = parser.parse_args()
    gen_csv(args)
    # accept MA, AR integers and beta0, beta1, beta2 coefficients, x1, x2 stddev and epsilon stdde
