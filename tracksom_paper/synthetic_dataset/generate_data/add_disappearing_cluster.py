"""
This file will add a new cluster that appear in day 2, stay in day 3, but disappeared in day 4.
"""

import pandas as pd
import numpy as np
import os

from pathlib import Path
from distributions_generators import gaussian_datapoints_generator, random_datapoints_generator


def create_new_population():
    # create the new population, just round in shape to make things simpler.
    new_population = gaussian_datapoints_generator(mean=[50, 50, 50],
                                                covariance=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                                                num_datapoints=100)
    return new_population


def write_gating_data(gating_data_dir, outdir, new_population_name, new_population,
                      days_to_add_population_to, noise_data):

    df_cluster_new = pd.DataFrame(new_population, columns=['x', 'y', 'z'])
    df_cluster_new = df_cluster_new.assign(PopName=[new_population_name] * df_cluster_new.shape[0])

    for d in days_to_add_population_to:
        df = pd.read_csv(gating_data_dir / 'synthetic_d{}.csv'.format(d))

        # remove noise data and add new ones
        df = df[(df['PopName'] != 'Noise')]
        df_noise = pd.DataFrame(noise_data[d], columns=['x', 'y', 'z'])
        df_noise = df_noise.assign(PopName=['Noise'] * df_noise.shape[0])
        df_new = df.append(df_noise)

        # append new population
        df_new = df_new.append(df_cluster_new)

        # shuffle rows
        df_new = df_new.sample(frac=1)
        df_new.to_csv(outdir / 'synthetic_d{}.csv'.format(d), index=False)


if __name__ == '__main__':
    outdir = Path('/Users/givanna/Documents/phd/dataset/synthetic_dataset/disappearing_cluster')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # read existing data
    orig_dataset_dir = Path('/Users/givanna/Documents/phd/dataset/synthetic_dataset/original_cc_paper')

    new_population = create_new_population()

    days_to_add_population_to = [2, 3]

    # we need to generate new noise data that encompasses the new region
    # first, read the existing data points in both days, excluding noise data as we don't want them.
    # use any gating file
    gating_data_dir = orig_dataset_dir / 'gating_coarse'
    existing_points = {}
    for d in days_to_add_population_to:
        df = pd.read_csv(gating_data_dir / 'synthetic_d{}.csv'.format(d))
        # remove noise
        df = df[(df['PopName'] != 'Noise')]
        points = df[['x','y','z']].to_numpy()
        existing_points[d] = points

    # generate new noise data
    noise_data = {}
    total_noise_count = 100
    for d in days_to_add_population_to:

        existing_points_for_this_day = np.concatenate([existing_points[d], new_population])
        min_datapoint_per_dimension = np.amin(existing_points_for_this_day, axis=0)
        max_datapoint_per_dimension = np.amax(existing_points_for_this_day, axis=0)
        noise_data_for_this_day = random_datapoints_generator(datapoints_to_exclude=existing_points_for_this_day,
                                                 num_datapoints=total_noise_count,
                                                 min_data=min_datapoint_per_dimension,
                                                 max_data=max_datapoint_per_dimension)
        noise_data[d] = noise_data_for_this_day

    # create the output directory for gating
    gating_outdir = outdir / 'gating_coarse'
    if not os.path.exists(gating_outdir):
        os.makedirs(gating_outdir)

    write_gating_data(gating_data_dir=gating_data_dir,
                      outdir=gating_outdir,
                      new_population_name='E',
                      new_population=new_population,
                      days_to_add_population_to=days_to_add_population_to,
                      noise_data=noise_data)

    gating_outdir = outdir / 'gating_fine'
    if not os.path.exists(gating_outdir):
        os.makedirs(gating_outdir)
    write_gating_data(gating_data_dir=orig_dataset_dir / 'gating_fine',
                      outdir=gating_outdir,
                      new_population_name='K',
                      new_population=new_population,
                      days_to_add_population_to=days_to_add_population_to,
                      noise_data=noise_data)





