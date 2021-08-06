import numpy as np
import matplotlib.pyplot as plt

#create observed data
mu, sigma = 0, 0.1 # mean and standard deviation
observed = np.random.normal(mu, sigma, 1000)


def find_distance(observed, simulated):
    distance = np.power(sum(observed), 2)+np.power(sum(simulated), 2)*(1/len(observed))
    return distance

def sampling(data, num_samples, threshold):
    posterior_distribution = []
    n = len(data)
    for i in range(0, num_samples):
        distance = threshold
        while distance >= threshold:
            variable = np.random.beta(1, 1, size=1)[0]
            simulated = np.random.normal(0, variable, n)
            distance = find_distance(simulated, data)
        posterior_distribution.append(variable)
    return posterior_distribution

posterior = sampling(observed, 1000, 1)
count, bins, ignored = plt.hist(posterior, 400, density=True)
plt.show()