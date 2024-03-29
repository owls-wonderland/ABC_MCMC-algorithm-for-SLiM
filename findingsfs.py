"""Code for finding parameters for background selection from SLiM"""
import msprime, pyslim
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import random as rd


def simulating_data(path, lenght_of_sfs, dominance_coefficient, mean, shape_parameter):
    """Function simulating_data runs Slim using terminal given number of times and generates .trees files.
    Extracts sfs from each .trees file. Later averages sfs and logs it."""
    samples = list(range(int(lenght_of_sfs)))
    random_seed_number = rd.randint(0, 2 ^ 62 - 1)
    try:
        # creating different .trees output files
        os.system(
            "slim -d seed="
            + str(random_seed_number)
            + " -d dominance_coefficient="
            + str(dominance_coefficient)
            + " -d mean="
            + str(mean)
            + " -d shape_parameter="
            + str(shape_parameter)
            + ' -d "path='
            + "'"
            + str(path)
            + "'\""
            + " background.slim"
        )

        # extracting sfs from each .trees file
        ts = pyslim.load(path + "/background" + str(random_seed_number) + ".trees")
        os.remove(
            path + "/background" + str(random_seed_number) + ".trees"
        )  # deleting used .trees file
        ts = ts.simplify(samples)
        sfs = ts.allele_frequency_spectrum(
            mode="branch", span_normalise=False, polarised=True
        )
        np.set_printoptions(threshold=np.inf)

        # Normalising SFS
        sumsfs = np.sum(sfs)
        sfs_norm = sfs / sumsfs

        # deleting first and last values
        modified_sfs = np.delete(sfs_norm, [0, len(sfs_norm) - 1])

        # Log
        log_sfs = np.log(modified_sfs) - np.log(1 - modified_sfs)
        return log_sfs, True
    except FileNotFoundError:
        return np.zeros(lenght_of_sfs - 1), False


def find_distance(observed, simulated):
    """Function find_distance calculates square distance between simulated and observed data"""
    distance = sum(np.power((observed - simulated), 2)) * (1 / len(observed))
    return distance


def sampling(path_sam, data, num_samples):
    # Finding length of sfs
    length = np.size(data)

    posterior_distribution_dominance = []
    posterior_distribution_mean = []
    posterior_distribution_shape = []

    # defining the very first random variables
    variable_dominance = 0.5
    variable_mean = -0.01
    variable_shape = 0.1

    # setting the first threshold
    succesful = False
    simulated = np.empty(length)
    simulated[:] = np.NaN
    while (np.all(np.isfinite(simulated)) == False) or (succesful == False):
        simulated, succesful = simulating_data(
            path_sam, length + 1, variable_dominance, variable_mean, variable_shape
        )
    # Saving simulated sfs into the matrix
    all_simulated_sfs = np.array([simulated])
    posterior_distribution_dominance.append(variable_dominance)
    posterior_distribution_mean.append(variable_mean)
    posterior_distribution_shape.append(variable_shape)

    # defining the very first threshold
    threshold = find_distance(data, simulated)
    # savind distances
    all_calculated_distances = np.array([threshold])

    # defining standard deviation and mean for sampling s+1 variable from s variable
    variable_dominance_st_dev = 0.1
    variable_mean_st_dev = 0.01
    variable_shape_st_dev = 0.1

    # defining target acceptance rate
    target_acceptance_rate = 0.1

    for i in range(0, num_samples):
        running_mean_of_dominance = variable_dominance
        running_mean_of_mean = variable_mean
        running_mean_of_shape = variable_shape

        # creating new variable for parameters

        new_variable_dominance = np.random.normal(
            variable_dominance, variable_dominance_st_dev, size=1
        )[0]
        new_variable_mean = np.random.normal(
            variable_mean, variable_mean_st_dev, size=1
        )[0]
        new_variable_shape = np.random.normal(
            variable_shape, variable_shape_st_dev, size=1
        )[0]

        # checking for negative values in parameter shape
        if new_variable_shape < 0:
            new_variable_shape = -new_variable_shape

        # finding simulated data & distance
        simulated[:] = np.NaN
        while np.all(np.isfinite(simulated)) == False:
            simulated, succesful = simulating_data(
                path_sam,
                length + 1,
                new_variable_dominance,
                new_variable_mean,
                new_variable_shape,
            )
        distance = find_distance(data, simulated)

        # accepting or rejecting new parameters
        if (distance <= threshold) & (succesful == True):
            # saving data
            all_simulated_sfs = np.vstack((all_simulated_sfs, simulated))
            all_calculated_distances = np.vstack((all_calculated_distances, distance))
            # when new variable is accepted
            posterior_distribution_mean.append(new_variable_mean)
            posterior_distribution_shape.append(new_variable_shape)
            posterior_distribution_dominance.append(new_variable_dominance)
            variable_mean = new_variable_mean
            variable_shape = new_variable_shape
            variable_dominance = new_variable_dominance
            threshold = math.exp(
                math.log(threshold) + (target_acceptance_rate - 1) / (i + 1)
            )
        else:
            # saving data
            print(succesful)
            all_simulated_sfs = np.vstack((all_simulated_sfs, all_simulated_sfs[-1, :]))
            all_calculated_distances = np.vstack(
                (all_calculated_distances, all_calculated_distances[-1])
            )
            # when variable is rejected
            posterior_distribution_mean.append(variable_mean)
            posterior_distribution_shape.append(variable_shape)
            posterior_distribution_dominance.append(variable_dominance)
            threshold = math.exp(
                math.log(threshold) + (target_acceptance_rate - 0) / (i + 1)
            )

        # Changing variables standard deviation
        if i > 1:
            running_mean_of_mean = running_mean_of_mean + (
                variable_mean - running_mean_of_mean
            ) / (i + 1)
            running_mean_of_shape = running_mean_of_shape + (
                variable_shape - running_mean_of_shape
            ) / (i + 1)
            running_mean_of_dominance = running_mean_of_dominance + (
                variable_dominance - running_mean_of_dominance
            ) / (i + 1)

            variable_mean_st_dev = math.sqrt(
                variable_mean_st_dev ** 2
                + (
                    (variable_mean - running_mean_of_mean) ** 2
                    - variable_mean_st_dev ** 2
                )
                / (i + 1)
            )
            variable_shape_st_dev = math.sqrt(
                variable_shape_st_dev ** 2
                + (
                    (variable_shape - running_mean_of_shape) ** 2
                    - variable_shape_st_dev ** 2
                )
                / (i + 1)
            )
            variable_dominance_st_dev = math.sqrt(
                variable_dominance_st_dev ** 2
                + (
                    (variable_dominance - running_mean_of_dominance) ** 2
                    - variable_dominance_st_dev ** 2
                )
                / (i + 1)
            )
    return (
        posterior_distribution_dominance,
        posterior_distribution_mean,
        posterior_distribution_shape,
        all_simulated_sfs,
        all_calculated_distances,
    )


def main():
    given_sfs = np.array(
        [
            -1.18783432,
            -1.74609835,
            -1.68097617,
            -2.63070513,
            -2.4895227,
            -3.51028512,
            -3.48102585,
            -3.88919502,
            -3.90956936,
            -4.1422421,
            -4.67984267,
            -4.88248064,
            -4.5274008,
            -4.21196954,
            -4.59830333,
            -4.4236753,
            -5.19989411,
            -5.10046828,
            -4.82342245,
            -4.81271866,
            -4.59234316,
            -4.91520473,
            -4.76823268,
            -4.43999834,
            -3.81415027,
            -4.02637127,
            -3.52964729,
            -5.35993467,
            -6.3311511,
        ]
    )
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    current_path = "/".join(arg for arg in args)
    posterior_dominance, posterior_mean, posterior_shape, simulated_sfs, calculated_distances = sampling(
        current_path, given_sfs, 40000
    )
    # making graphs
    count, bins, ignored = plt.hist(posterior_mean, 100, density=True, label="mean")
    count, bins, ignored = plt.hist(
        posterior_shape, 100, density=True, label="shape parameter"
    )
    count, bins, ignored = plt.hist(
        posterior_dominance, 100, density=True, label="dominance parameter"
    )
    plt.title("Parameters frequency")
    plt.xlabel("value")
    plt.ylabel("frequency")
    plt.legend()
    plt.savefig("mcmc_parameters.png")

    # saving data in matrix and turning into dataframe
    simulated_sfs = np.column_stack((simulated_sfs, np.array(posterior_dominance)))
    simulated_sfs = np.column_stack((simulated_sfs, np.array(posterior_mean)))
    simulated_sfs = np.column_stack((simulated_sfs, np.array(posterior_shape)))
    simulated_sfs = np.column_stack((simulated_sfs, calculated_distances))
    data_about_simulations = pd.DataFrame(simulated_sfs)
    pd.DataFrame(data_about_simulations).to_csv("All_used_data.csv")


if __name__ == "__main__":
    main()
