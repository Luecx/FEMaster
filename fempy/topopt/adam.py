import numpy as np

class Adam:
    def __init__(self, initial_values, learning_rate=0.001, beta1=0.9, beta2=0.999, epsilon=1e-8):
        """
        Initializes the Adam optimizer.

        Parameters:
        - initial_values (np.ndarray): Initial values for optimization.
        - learning_rate (float): Step size for the updates.
        - beta1 (float): Exponential decay rate for the first moment estimates.
        - beta2 (float): Exponential decay rate for the second moment estimates.
        - epsilon (float): Small constant for numerical stability.
        """
        self.values = np.array(initial_values, dtype=np.float64)
        self.learning_rate = learning_rate
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon

        # Initialize first moment vector, second moment vector, and timestep
        self.m = np.zeros_like(self.values)
        self.v = np.zeros_like(self.values)
        self.t = 0

    def step(self, gradients):
        """
        Perform a single optimization step using Adam.

        Parameters:
        - gradients (np.ndarray): Gradients for the current step.

        Returns:
        - np.ndarray: Updated values after the optimization step.
        """
        self.t += 1  # Increment timestep

        # Update biased first moment estimate
        self.m = self.beta1 * self.m + (1 - self.beta1) * gradients

        # Update biased second raw moment estimate
        self.v = self.beta2 * self.v + (1 - self.beta2) * (gradients ** 2)

        # Compute bias-corrected first moment estimate
        m_hat = self.m / (1 - self.beta1 ** self.t)

        # Compute bias-corrected second moment estimate
        v_hat = self.v / (1 - self.beta2 ** self.t)

        # Update the parameters
        self.values -= self.learning_rate * m_hat / (np.sqrt(v_hat) + self.epsilon)

        return self.values

# Example usage
if __name__ == "__main__":
    # Initial orientations or parameters
    initial_orientations = np.random.rand(10, 3)  # Example: 10 elements with 3 orientation values each

    # Gradients (dummy example)
    gradients = np.random.randn(10, 3)

    # Instantiate the Adam optimizer
    optimizer = Adam(initial_orientations, learning_rate=0.01)

    # Perform a few steps
    for i in range(5):
        updated_orientations = optimizer.step(gradients)
        print(f"Step {i + 1}:")
        print(updated_orientations)
