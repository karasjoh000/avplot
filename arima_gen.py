import numpy as np
import pandas as pd

# Parameters
b0, b1, b2 = 2, 1.5, -0.8
phi = 0.6
theta1, theta2, theta3 = 0.4, -0.2, 0.1
n = 50  # Number of data points


# Generate exogenous variables
np.random.seed(42)  # For reproducibility
x1 = np.random.normal(0, 1, n)
x2 = np.random.normal(0, 1, n)

# Generate ARIMA(1, 0, 3) process
epsilon = np.random.normal(0, 1, n + 3)  # Extra points for lagged errors
m = np.zeros(n + 3)


for t in range(3, n + 3):
    m[t] = phi * m[t - 1] + epsilon[t] + theta1 * epsilon[t - 1] + theta2 * epsilon[t - 2] + theta3 * epsilon[t - 3]

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

# Save the data to a CSV file
data.to_csv('arima_data.csv', index=False)

# Display the first 10 rows of the data
print(data.head(10))
