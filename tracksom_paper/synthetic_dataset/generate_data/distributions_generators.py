import numpy as np


def random_datapoints_generator(datapoints_to_exclude, num_datapoints, min_data, max_data):
    # ugly as hell, but otherwise how do u do a random sample excluding the cluster?
    noise_data = []
    for i in range(0, num_datapoints):
        n = None
        while n is None or n in noise_data or n in datapoints_to_exclude:
            n = [np.random.uniform(min_data[0], max_data[0], 1)[0],
                 np.random.uniform(min_data[1], max_data[1], 1)[0],
                 np.random.uniform(min_data[2], max_data[2], 1)[0]]
        noise_data.append(n)
    return np.array(noise_data)


def gaussian_datapoints_generator(mean, covariance, num_datapoints):
    """
    Generate random data points based on Gaussian Distribution

    Parameters
    ----------
    mean
    covariance
    num_datapoints

    Returns
    -------

    """
    data = np.random.multivariate_normal(mean, covariance, num_datapoints)
    return data