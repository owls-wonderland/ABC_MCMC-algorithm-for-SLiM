"""Code for finding parameters for background selection from SLiM"""
import msprime, pyslim
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys


def simulating_data(
    path, lenght_of_sfs, number_of_sfs, dominance_coefficient, mean, shape_parameter
):
    """Function simulating_data runs Slim using terminal given number of times and generates .trees files.
    Extracts sfs from each .trees file. Later averages sfs and logs it."""
    samples = list(range(int(lenght_of_sfs)))
    for i in range(number_of_sfs):
        # creating different .trees output files
        os.system(
            "slim -d seed="
            + str(i)
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
        ts = pyslim.load(path + "/background" + str(i) + ".trees")
        os.remove(path + "/background" + str(i) + ".trees") # deleting used .trees file
        ts = ts.simplify(samples)
        sfs = ts.allele_frequency_spectrum(
            mode="branch", span_normalise=False, polarised=True
        )
        np.set_printoptions(threshold=np.inf)

        # Normalising SFS
        sumsfs = np.sum(sfs)
        sfs_norm = sfs / sumsfs
        if i == 0:
            all_sfs = sfs_norm
        else:
            all_sfs = all_sfs + sfs_norm

    # average sfs
    modified_sfs = np.delete(all_sfs, [0, len(all_sfs) - 1])
    modified_sfs = modified_sfs / number_of_sfs

    # Log
    log_sfs = np.log(modified_sfs) - np.log(1 - modified_sfs)
    return log_sfs


def find_distance(observed, simulated):
    """Function find_distance calculates square distance between simulated and observed data"""
    distance = sum(np.power((observed - simulated), 2)) * (1 / len(observed))
    return distance


def sampling(path_sam, data, num_samples):
    posterior_distribution_dominance = []
    posterior_distribution_mean = []
    posterior_distribution_shape = []

    # defining the very first random variables
    variable_dominance = 0.5
    variable_mean = 0.3
    variable_shape = 0.3

    # setting the first threshold
    simulated = simulating_data(
        path_sam, 30, 10, variable_dominance, variable_mean, variable_shape
    )
    print("very first one simulated")
    threshold = find_distance(data, simulated)

    # defining standard deviation and mean for sampling s+1 variable from s variable
    variable_dominance_st_dev = 0.1
    variable_mean_st_dev = 0.1
    variable_shape_st_dev = 0.1

    # defining target acceptance rate
    target_acceptance_rate = 0.1

    for i in range(0, num_samples):
        running_mean_of_dominance = variable_dominance
        running_mean_of_mean = variable_mean
        running_mean_of_shape = variable_shape

        # creating new variable for parameters

        # new_variable_dominance = np.random.normal(
        #    variable_dominance, variable_dominance, size=1
        # )[0]
        new_variable_dominance = 0.5
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
        simulated = simulating_data(
            path_sam,
            30,
            10,
            new_variable_dominance,
            new_variable_mean,
            new_variable_shape,
        )
        print(i)
        print(
            "simulated once------------------------------------------------------------------"
        )
        distance = find_distance(data, simulated)

        # accepting or rejecting new parameters
        if distance <= threshold:
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
    )


"""
# This code was used to create given_sfs array
sfs = simulating_data(30, 20, 0.5, -0.01, 0.1)
"""
given_sfs = np.array(
    [
        -0.97521618
        - 1.87011656
        - 2.24700429
        - 2.65951118
        - 3.05000166
        - 3.07489556
        - 3.04522172
        - 3.47843442
        - 3.69398409
        - 4.06276335
        - 3.85851225
        - 4.29132519
        - 4.07425373
        - 4.09466838
        - 4.16694236
        - 4.25458452
        - 4.24621578
        - 4.51276235
        - 4.36592735
        - 4.62343833
        - 4.57562728
        - 4.64440639
        - 4.41582723
        - 4.57696382
        - 4.9021774
        - 4.92645268
        - 4.44071992
        - 4.94018856
        - 4.98841237
    ]
)


def main():
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    current_path = "/".join(arg for arg in args)
    posterior_dominance, posterior_mean, posterior_shape = sampling(
        current_path, given_sfs, 1
    )
    count, bins, ignored = plt.hist(posterior_mean, 100, density=True)
    count, bins, ignored = plt.hist(posterior_shape, 100, density=True)
    # count, bins, ignored = plt.hist(posterior_dominance, 100, density=True)
    plt.title("Parameters frequency")
    plt.xlabel("value")
    plt.ylabel("frequency")
    plt.savefig("mcmc_parameters.png")


if __name__ == "__main__":
    main()
