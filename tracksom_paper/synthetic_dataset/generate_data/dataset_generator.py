# Code to generate the synthetic dataset used in the ChronoClust paper and subsequently adapted for the TrackSOM paper.

import numpy as np
import sys
import pandas as pd

from distributions_generators import gaussian_datapoints_generator, random_datapoints_generator


"""
This is essentially version 4 but with the final day kinks all added.
"""

TOTAL_NOISE_PTS = 100
TOTAL_NUM_PTS_CLUSTER1 = 5000
TOTAL_NUM_PTS_CLUSTER2 = 2000

CLUSTER1_ID = 1
CLUSTER2_ID = 2
NOISE_ID = -1


def run(csv_file_directory):
    generate_dataset_for_day0(csv_file_directory)
    generate_dataset_for_day1(csv_file_directory)
    generate_dataset_for_day2(csv_file_directory)
    generate_dataset_for_day3(csv_file_directory)
    generate_dataset_for_day4(csv_file_directory)


def generate_dataset_for_day0(csv_file_directory):
    """
    Day 0
    2 clusters very far away from one another.
    Both clusters are of spherical shape.
    Cluster 1 contains 1000 datapoints are located at bottom left.
    Cluster 2 contains 2000 datapoints are located at top right.
    """
    covariance = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cluster1 = gaussian_datapoints_generator(mean=[10, 10, 10],
                                             covariance=covariance,
                                             num_datapoints=TOTAL_NUM_PTS_CLUSTER1)

    cluster2 = gaussian_datapoints_generator(mean=[30, 30, 30],
                                                        covariance=covariance,
                                                        num_datapoints=TOTAL_NUM_PTS_CLUSTER2)
    # Get noise datapoints
    all_datapoints = np.concatenate((cluster1, cluster2), axis=0)
    min_datapoint_per_dimension = np.amin(all_datapoints, axis=0)
    max_datapoint_per_dimension = np.amax(all_datapoints, axis=0)
    noise = random_datapoints_generator(datapoints_to_exclude=all_datapoints, num_datapoints=TOTAL_NOISE_PTS,
                                                   min_data=min_datapoint_per_dimension,
                                                   max_data=max_datapoint_per_dimension)

    # Label the datapoints with cluster it belongs to
    cluster1 = np.insert(cluster1, 3, CLUSTER1_ID, axis=1)
    cluster2 = np.insert(cluster2, 3, CLUSTER2_ID, axis=1)
    noise = np.insert(noise, 3, NOISE_ID, axis=1)

    all_datapoints = np.concatenate((cluster1, cluster2, noise), axis=0)

    # Shuffle the points
    np.random.shuffle(all_datapoints)

    csv_file = '{}/synthetic_d0.csv'.format(csv_file_directory)

    # Add the header and write out
    df = pd.DataFrame(all_datapoints, columns=['x', 'y', 'z', 'PopName'])
    df.to_csv(csv_file, index=False)


def generate_dataset_for_day1(csv_file_directory):
    """
    Day 1:
    2 clusters. One at bottom left (cluster 1), the other at top right (cluster 2).
    See data_points for day 0 for previous day.

    Cluster 1 move and kind of split in positive direction for each axes. So think of it like a dumbbell shape with very short handle in every direction.
    One end of the dumbbell will be shared among the 3.
    We call the non-shared dumbbell plate branch_off_cluster
    Datapoints need to be maintained at 5000 across all the dumbbell. Since one plate will be shared among the 3, we can have about 1000 for each other plate.

    Cluster 2 changes shape to be more elongated in x direction. Maybe a little bit in y and z as well. Datapoints need to be maintained at 2000 still.
    """

    def generate_cluster1():
        covariance_cluster1 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        # Work out the points allocation for cluster 1

        num_points_branched_off_cluster = 1000
        num_points_original_cluster = TOTAL_NUM_PTS_CLUSTER1 - (3 * num_points_branched_off_cluster)
        # Generate the dumbbells' plates
        cluster1_original = gaussian_datapoints_generator(mean=[10, 10, 10],
                                                                         covariance=covariance_cluster1,
                                                                         num_datapoints=num_points_original_cluster)
        cluster1_branched_off_x = gaussian_datapoints_generator(mean=[12, 10, 10],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_y = gaussian_datapoints_generator(mean=[10, 12, 10],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_z = gaussian_datapoints_generator(mean=[10, 10, 12],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)

        cluster1_without_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)

        # Label the datapoints with cluster it belongs to
        cluster1_original = np.insert(cluster1_original, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_x = np.insert(cluster1_branched_off_x, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_y = np.insert(cluster1_branched_off_y, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_z = np.insert(cluster1_branched_off_z, 3, CLUSTER1_ID, axis=1)

        cluster1_with_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)
        return cluster1_with_label, cluster1_without_label

    cluster1, cluster1_without_label = generate_cluster1()

    covariance_cluster2 = [[10.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cluster2 = gaussian_datapoints_generator(mean=[30, 30, 30],
                                                        covariance=covariance_cluster2,
                                                        num_datapoints=TOTAL_NUM_PTS_CLUSTER2)

    # Get noise datapoints
    all_datapoints = np.concatenate((cluster1_without_label, cluster2), axis=0)
    min_datapoint_per_dimension = np.amin(all_datapoints, axis=0)
    max_datapoint_per_dimension = np.amax(all_datapoints, axis=0)
    noise = random_datapoints_generator(datapoints_to_exclude=all_datapoints, num_datapoints=TOTAL_NOISE_PTS,
                                                   min_data=min_datapoint_per_dimension,
                                                   max_data=max_datapoint_per_dimension)

    # Label the datapoints with cluster it belongs to
    cluster2 = np.insert(cluster2, 3, CLUSTER2_ID, axis=1)
    noise = np.insert(noise, 3, NOISE_ID, axis=1)

    all_datapoints = np.concatenate((cluster1, cluster2, noise), axis=0)

    # Shuffle the points
    np.random.shuffle(all_datapoints)

    csv_file = '{}/synthetic_d1.csv'.format(csv_file_directory)

    # Add the header and write out
    df = pd.DataFrame(all_datapoints, columns=['x', 'y', 'z', 'PopName'])
    df.to_csv(csv_file, index=False)


def generate_dataset_for_day2(csv_file_directory):
    """
    Day 2:
    2 clusters. One at bottom left (cluster 1), the other at top right (cluster 2).
    See data_points for day 0 for previous day.

    Cluster 1 move and kind of split in positive direction for each axes. So think of it like a dumbbell shape with very short handle in every direction.
    One end of the dumbbell will be shared among the 3.
    We call the non-shared dumbbell plate branch_off_cluster
    Datapoints need to be maintained at 5000 across all the dumbbell. Since one plate will be shared among the 3, we can have about 1000 for each other plate.

    Cluster 2 changes shape to be more elongated in x direction. Maybe a little bit in y and z as well. Datapoints need to be maintained at 2000 still.
    """

    def generate_cluster1():
        covariance_cluster1 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        # Work out the points allocation for cluster 1

        num_points_branched_off_cluster = 1000
        num_points_original_cluster = TOTAL_NUM_PTS_CLUSTER1 - (3 * num_points_branched_off_cluster)
        # Generate the dumbbells' plates
        cluster1_original = gaussian_datapoints_generator(mean=[10, 10, 10],
                                                                         covariance=covariance_cluster1,
                                                                         num_datapoints=num_points_original_cluster)
        cluster1_branched_off_x = gaussian_datapoints_generator(mean=[13, 10, 10],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_y = gaussian_datapoints_generator(mean=[10, 13, 10],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_z = gaussian_datapoints_generator(mean=[10, 10, 13],
                                                                    covariance=covariance_cluster1,
                                                                    num_datapoints=num_points_branched_off_cluster)

        cluster1_without_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)

        # Label the datapoints with cluster it belongs to
        cluster1_original = np.insert(cluster1_original, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_x = np.insert(cluster1_branched_off_x, 3, CLUSTER1_ID + 0.1, axis=1)
        cluster1_branched_off_y = np.insert(cluster1_branched_off_y, 3, CLUSTER1_ID + 0.2, axis=1)
        cluster1_branched_off_z = np.insert(cluster1_branched_off_z, 3, CLUSTER1_ID + 0.3, axis=1)

        cluster1_with_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)
        return cluster1_with_label, cluster1_without_label

    cluster1, cluster1_without_label = generate_cluster1()

    covariance_cluster2 = [[15.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cluster2 = gaussian_datapoints_generator(mean=[30, 30, 30],
                                                        covariance=covariance_cluster2,
                                                        num_datapoints=TOTAL_NUM_PTS_CLUSTER2)

    # Get noise datapoints
    all_datapoints = np.concatenate((cluster1_without_label, cluster2), axis=0)
    min_datapoint_per_dimension = np.amin(all_datapoints, axis=0)
    max_datapoint_per_dimension = np.amax(all_datapoints, axis=0)
    noise = random_datapoints_generator(datapoints_to_exclude=all_datapoints, num_datapoints=TOTAL_NOISE_PTS,
                                                   min_data=min_datapoint_per_dimension,
                                                   max_data=max_datapoint_per_dimension)

    # Label the datapoints with cluster it belongs to
    cluster2 = np.insert(cluster2, 3, CLUSTER2_ID, axis=1)
    noise = np.insert(noise, 3, NOISE_ID, axis=1)

    all_datapoints = np.concatenate((cluster1, cluster2, noise), axis=0)

    # Shuffle the points
    np.random.shuffle(all_datapoints)

    csv_file = '{}/synthetic_d2.csv'.format(csv_file_directory)

    # Add the header and write out
    df = pd.DataFrame(all_datapoints, columns=['x', 'y', 'z', 'PopName'])
    df.to_csv(csv_file, index=False)


def generate_dataset_for_day3(csv_file_directory):
    """
    Day 3:
    2 clusters. One at bottom left (cluster 1), the other at top right (cluster 2).
    See data_points for day 1 for previous day.

    Cluster 1 moved further away, to the point where they are very separated.

    Cluster 2 split into cluster 2a and 2b each containing 1000 data_points points
    """

    def generate_cluster1():
        # Work out the points allocation for cluster 1
        num_points_branched_off_cluster = 1000
        num_points_original_cluster = TOTAL_NUM_PTS_CLUSTER1 - (3 * num_points_branched_off_cluster)

        # Generate the dumbbells' plates
        cluster1_original = gaussian_datapoints_generator(mean=[10, 10, 10],
                                                                         covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                     [0.0, 0.0, 1.0]],
                                                                         num_datapoints=num_points_original_cluster)
        cluster1_branched_off_x = gaussian_datapoints_generator(mean=[15, 10, 10],
                                                                    covariance=[[3.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_y = gaussian_datapoints_generator(mean=[10, 15, 10],
                                                                    covariance=[[1.0, 0.0, 0.0], [0.0, 3.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=num_points_branched_off_cluster)
        cluster1_branched_off_z = gaussian_datapoints_generator(mean=[10, 10, 15],
                                                                    covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                [0.0, 0.0, 3.0]],
                                                                    num_datapoints=num_points_branched_off_cluster)

        cluster1_without_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)

        # Label the datapoints with cluster it belongs to
        cluster1_original = np.insert(cluster1_original, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_x = np.insert(cluster1_branched_off_x, 3, CLUSTER1_ID + 0.1, axis=1)
        cluster1_branched_off_y = np.insert(cluster1_branched_off_y, 3, CLUSTER1_ID + 0.2, axis=1)
        cluster1_branched_off_z = np.insert(cluster1_branched_off_z, 3, CLUSTER1_ID + 0.3, axis=1)

        cluster1_with_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x, cluster1_branched_off_y, cluster1_branched_off_z), axis=0)
        return cluster1_with_label, cluster1_without_label

    cluster1_with_label, cluster1_without_label = generate_cluster1()

    # Cluster 2 split into 2.
    covariance = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cluster2a = gaussian_datapoints_generator(mean=[33, 30, 30],
                                                         covariance=covariance,
                                                         num_datapoints=int(TOTAL_NUM_PTS_CLUSTER2 / 2))
    cluster2b = gaussian_datapoints_generator(mean=[23, 30, 30],
                                                         covariance=covariance,
                                                         num_datapoints=int(TOTAL_NUM_PTS_CLUSTER2 / 2))

    # Get noise datapoints
    all_datapoints = np.concatenate((cluster1_without_label, cluster2a, cluster2b), axis=0)
    min_datapoint_per_dimension = np.amin(all_datapoints, axis=0)
    max_datapoint_per_dimension = np.amax(all_datapoints, axis=0)
    noise = random_datapoints_generator(datapoints_to_exclude=all_datapoints, num_datapoints=TOTAL_NOISE_PTS,
                                                   min_data=min_datapoint_per_dimension,
                                                   max_data=max_datapoint_per_dimension)

    # Label the datapoints with cluster it belongs to
    cluster2a = np.insert(cluster2a, 3, CLUSTER2_ID + 0.1, axis=1)
    cluster2b = np.insert(cluster2b, 3, CLUSTER2_ID + 0.2, axis=1)
    noise = np.insert(noise, 3, NOISE_ID, axis=1)

    all_datapoints = np.concatenate((cluster1_with_label, cluster2a, cluster2b, noise), axis=0)

    # Shuffle the points
    np.random.shuffle(all_datapoints)

    csv_file = '{}/synthetic_d3.csv'.format(csv_file_directory)
    # Add the header and write out
    df = pd.DataFrame(all_datapoints, columns=['x', 'y', 'z', 'PopName'])
    df.to_csv(csv_file, index=False)


def generate_dataset_for_day4(csv_file_directory):
    """
    Day 3:
    2 clusters. One at bottom left (cluster 1), the other at top right (cluster 2).
    See data_points for day 1 for previous day.

    Cluster 1 moved slightly further away. The tip of each extension (plate in day 1) has is starting to split and grow in different direction.

    Cluster 2a and 2b moved further apart.
    """

    def generate_cluster1():
        # Work out the points allocation for cluster 1
        num_points_branched_off_cluster = 1000
        num_points_original_cluster = TOTAL_NUM_PTS_CLUSTER1 - (3 * num_points_branched_off_cluster)

        # Generate the dumbbells' plates
        cluster1_original = gaussian_datapoints_generator(mean=[10, 10, 10],
                                                                         covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                     [0.0, 0.0, 1.0]],
                                                                         num_datapoints=num_points_original_cluster)
        cluster1_branched_off_x_1 = gaussian_datapoints_generator(mean=[15, 10, 10],
                                                                    covariance=[[3.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=int(num_points_branched_off_cluster/2))
        cluster1_branched_off_x_2 = gaussian_datapoints_generator(mean=[18, 7, 10],
                                                                    covariance=[[1.0, 0.0, 0.0], [0.0, 5.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=int(num_points_branched_off_cluster/2))
        cluster1_branched_off_y_1 = gaussian_datapoints_generator(mean=[10, 15, 10],
                                                                    covariance=[[1.0, 0.0, 0.0], [0.0, 3.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=int(num_points_branched_off_cluster/2))
        cluster1_branched_off_y_2 = gaussian_datapoints_generator(mean=[10, 18, 7],
                                                                      covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                  [0.0, 0.0, 5.0]],
                                                                      num_datapoints=int(
                                                                          num_points_branched_off_cluster / 2))
        cluster1_branched_off_z_1 = gaussian_datapoints_generator(mean=[10, 10, 15],
                                                                    covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                [0.0, 0.0, 3.0]],
                                                                    num_datapoints=int(num_points_branched_off_cluster/2))
        cluster1_branched_off_z_2 = gaussian_datapoints_generator(mean=[7, 10, 18],
                                                                    covariance=[[5.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                                                                [0.0, 0.0, 1.0]],
                                                                    num_datapoints=int(num_points_branched_off_cluster/2))

        cluster1_without_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x_1, cluster1_branched_off_x_2,
             cluster1_branched_off_y_1, cluster1_branched_off_y_2,
             cluster1_branched_off_z_1, cluster1_branched_off_z_2), axis=0)

        # Label the datapoints with cluster it belongs to
        cluster1_original = np.insert(cluster1_original, 3, CLUSTER1_ID, axis=1)
        cluster1_branched_off_x_1 = np.insert(cluster1_branched_off_x_1, 3, CLUSTER1_ID + 0.11, axis=1)
        cluster1_branched_off_x_2 = np.insert(cluster1_branched_off_x_2, 3, CLUSTER1_ID + 0.12, axis=1)
        cluster1_branched_off_y_1 = np.insert(cluster1_branched_off_y_1, 3, CLUSTER1_ID + 0.21, axis=1)
        cluster1_branched_off_y_2 = np.insert(cluster1_branched_off_y_2, 3, CLUSTER1_ID + 0.22, axis=1)
        cluster1_branched_off_z_1 = np.insert(cluster1_branched_off_z_1, 3, CLUSTER1_ID + 0.31, axis=1)
        cluster1_branched_off_z_2 = np.insert(cluster1_branched_off_z_2, 3, CLUSTER1_ID + 0.32, axis=1)

        cluster1_with_label = np.concatenate(
            (cluster1_original, cluster1_branched_off_x_1, cluster1_branched_off_x_2,
             cluster1_branched_off_y_1, cluster1_branched_off_y_2,
             cluster1_branched_off_z_1, cluster1_branched_off_z_2), axis=0)
        return cluster1_with_label, cluster1_without_label

    cluster1_with_label, cluster1_without_label = generate_cluster1()

    # Cluster 2 split into 2.
    covariance = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    cluster2a = gaussian_datapoints_generator(mean=[35, 30, 30],
                                                         covariance=covariance,
                                                         num_datapoints=int(TOTAL_NUM_PTS_CLUSTER2 / 2))
    cluster2b = gaussian_datapoints_generator(mean=[21, 30, 30],
                                                         covariance=covariance,
                                                         num_datapoints=int(TOTAL_NUM_PTS_CLUSTER2 / 2))

    # Get noise datapoints
    all_datapoints = np.concatenate((cluster1_without_label, cluster2a, cluster2b), axis=0)
    min_datapoint_per_dimension = np.amin(all_datapoints, axis=0)
    max_datapoint_per_dimension = np.amax(all_datapoints, axis=0)
    noise = random_datapoints_generator(datapoints_to_exclude=all_datapoints, num_datapoints=TOTAL_NOISE_PTS,
                                                   min_data=min_datapoint_per_dimension,
                                                   max_data=max_datapoint_per_dimension)

    # Label the datapoints with cluster it belongs to
    cluster2a = np.insert(cluster2a, 3, CLUSTER2_ID + 0.1, axis=1)
    cluster2b = np.insert(cluster2b, 3, CLUSTER2_ID + 0.2, axis=1)
    noise = np.insert(noise, 3, NOISE_ID, axis=1)

    all_datapoints = np.concatenate((cluster1_with_label, cluster2a, cluster2b, noise), axis=0)

    # Shuffle the points
    np.random.shuffle(all_datapoints)

    # Save the points
    csv_file = '{}/synthetic_d4.csv'.format(csv_file_directory)

    # Add the header and write out
    df = pd.DataFrame(all_datapoints, columns=['x','y','z','PopName'])
    df.to_csv(csv_file, index=False)


if __name__ == '__main__':
    # Version 4 is where cluster 1 (blue) is version 3 and cluster 2 (red) is version 2.
    if len(sys.argv) < 2:
        sys.exit("No csv files directory to store the dataset is given.")

    datafile_directory = sys.argv[1]
    run(datafile_directory)
