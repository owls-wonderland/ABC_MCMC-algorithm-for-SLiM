import numpy as np
import math
import matplotlib.pyplot as plt

# create observed data
mu, sigma = 0, 0.1  # mean and standard deviation
observed = np.random.normal(mu, sigma, 100)


def find_distance(observed, simulated):
    distance = sum(np.power((observed - simulated), 2)) * (1 / len(observed))
    return distance


def sampling(data, num_samples):
    posterior_distribution_mean = []
    posterior_distribution_sd = []
    n = len(data)
    # defining the very first random variable and setting the threshold equal to distance of the first variable, so it would be accepted
    variable_mean = 1
    variable_sd = 1
    simulated = np.random.normal(variable_mean, variable_sd, n)
    threshold = find_distance(simulated, data)
    # defining standard deviation and mean for sampling s+1 variable from s variable
    variable_mean_st_dev = 0.1
    variable_sd_st_dev = 0.1
    # defining target acceptance rate
    target_acceptance_rate = 0.1
    for i in range(0, num_samples):
        running_mean_of_mean = variable_mean
        running_mean_of_sd = variable_sd
        new_variable_mean = np.random.normal(
            variable_mean, variable_mean_st_dev, size=1
        )[0]
        new_variable_sd = np.random.normal(variable_sd, variable_sd_st_dev, size=1)[0]
        if new_variable_sd < 0:
            new_variable_sd = -new_variable_sd
        simulated = np.random.normal(new_variable_mean, new_variable_sd, n)
        distance = find_distance(simulated, data)
        if distance <= threshold:
            # when new variable is accepted
            posterior_distribution_mean.append(new_variable_mean)
            posterior_distribution_sd.append(new_variable_sd)
            variable_mean = new_variable_mean
            variable_sd = new_variable_sd
            threshold = math.exp(
                math.log(threshold) + (target_acceptance_rate - 1) / (i + 1)
            )
        else:
            # when variable is rejected
            posterior_distribution_mean.append(variable_mean)
            posterior_distribution_sd.append(variable_sd)
            threshold = math.exp(
                math.log(threshold) + (target_acceptance_rate - 0) / (i + 1)
            )

        # Changing variables standard deviation
        if i > 1:
            running_mean_of_mean = running_mean_of_mean + (
                variable_mean - running_mean_of_mean
            ) / (i + 1)
            running_mean_of_sd = running_mean_of_sd + (
                variable_sd - running_mean_of_sd
            ) / (i + 1)
            variable_mean_st_dev = math.sqrt(
                variable_mean_st_dev ** 2
                + (
                    (variable_mean - running_mean_of_mean) ** 2
                    - variable_mean_st_dev ** 2
                )
                / (i + 1)
            )
            variable_sd_st_dev = math.sqrt(
                variable_sd_st_dev ** 2
                + ((variable_sd - running_mean_of_sd) ** 2 - variable_sd_st_dev ** 2)
                / (i + 1)
            )
    return posterior_distribution_mean, posterior_distribution_sd


posterior_mean, posterior_sd = sampling(observed, 1000000)
# posterior_sd = sampling(observed, 1000000)[1]
count, bins, ignored = plt.hist(posterior_mean, 100, density=True)
count, bins, ignored = plt.hist(posterior_sd, 100, density=True)
plt.show()
