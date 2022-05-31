import numpy as np
import pandas as pd
import random
import sys
import math


def outliers_of_df(df):
    """
    Return outliers in DataFrame based on interquartile range (IQR).
    https://careerfoundry.com/en/blog/data-analytics/how-to-find-outliers/
    """
    q1 = df.quantile(0.20)
    q3 = df.quantile(0.80)
    iqr = q3 - q1
    
    outliers = df[((df < (q1 - 1.5 * iqr)) | (df > (q3 + 1.5 * iqr)))]

    return outliers


def drop_outliers(df, column):
    """
    Return df with dropped outliers with respect to given column.
    """
    return df.drop(outliers_of_df(df[column]).index.tolist())


def mean_std_sp_to_hsc_cardinality(sticky_paper_df, high_speed_camera_df, sample_count):
    """
    Return mean and std of the random quotient of cardinalities of SP and HSC splashes.

    Since there are only 16 splashes recorded by HSC we take a random sample of 16 splashes
    recorder by SP and divide the sum of beads in these 16 splashes by the sum of beads in
    all HSC splashes.
    """
    
    high_speed_camera_sum = sum([len(high_speed_camera_df[high_speed_camera_df['no'] == i]) for i in range(16)])
    sticky_paper_cardinalities = [len(sticky_paper_df[sticky_paper_df['no'] == i]) for i in range(48)]
    sp_hsc_quotient = [sum(random.sample(sticky_paper_cardinalities, 16)) / high_speed_camera_sum for i in range(sample_count)]

    return np.mean(sp_hsc_quotient), np.std(sp_hsc_quotient)


def shift_splashes_no(splashes_df):
    """
    Index splashes from 0 and not from 1.
    """
    splashes_df['no'] -= 1


def read_arguments():
    """
    Return CSV DataFrames for HSC and SP data and the number of quants of energy
    """
    return pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2]), int(sys.argv[3])


class EnergyStats():
    def __init__(self, splashes_hsc_df, splashes_sp_df, number_of_quants_of_energy):
        self.sum = [splashes_hsc_df[splashes_hsc_df['no'] == i]['e'].sum() for i in range(16)]
        self.mean = np.mean(self.sum)
        self.std = np.std(self.sum)
        self.min = np.min(self.sum)
        self.max = np.max(self.sum)
        self.min_particle = np.min(high_speed_camera_df['e'])
        self.max_particle = np.max(high_speed_camera_df['e'])
        self.scaling_mean, self.scaling_std = mean_std_sp_to_hsc_cardinality(splashes_sp_df, splashes_hsc_df, 1000)
        self.scaled_mean = self.mean * self.scaling_mean
        self.scaled_std = self.std * self.scaling_mean
        self.scaled_min = self.min * self.scaling_mean
        self.scaled_max = self.max * self.scaling_mean
        self.quant_of_energy = self.scaled_mean / number_of_quants_of_energy
        self.quantized_mean = self.scaled_mean / self.quant_of_energy
        self.quantized_std = self.scaled_std / self.quant_of_energy
        self.quantized_min = self.scaled_min / self.quant_of_energy
        self.quantized_max = self.scaled_max / self.quant_of_energy
        self.min_of_one_part = math.floor(self.min_particle / self.quant_of_energy)
        self.max_of_one_part = math.ceil(self.max_particle / self.quant_of_energy)
        self.min_number_of_parts = math.ceil(self.quantized_min / self.max_of_one_part)
        self.max_number_of_parts = math.floor(self.quantized_max  / self.min_of_one_part)

    def print_stats():
        print(f"Energy mean = {self.mean}, std = {self.std}")
        print(f"Energy min = {self.min}, max = {self.max}")
        print(f"Particle energy min = {self.min_particle}, max = {self.max_particle}")
        print(f"Scaling mean = {self.scaling_mean}, std = {self.scaling_std}, coeff of var = {self.scaling_std/self.scaling_mean}")
        print(f"Scaled energy mean = {self.scaled_mean}, std = {self.scaled_std}")
        print(f"Scaled energy min = {self.scaled_min}, std = {self.scaled_max}")
        print(f"Quant energy mean = {self.quantized_mean}, std = {self.quantized_std}")
        print(f"Quant energy min = {self.quantized_min}, max = {self.quantized_max}")
        print(f"kmin = {self.min_of_one_part}, kmax = {self.max_of_one_part}")
        print(f"parts min = {self.min_number_of_parts}, max = {self.max_number_of_parts}")


high_speed_camera_df, sticky_paper_df, number_of_q_energy = read_arguments()

shift_splashes_no(high_speed_camera_df)
shift_splashes_no(sticky_paper_df)

# indices of sticky paper registered by hsc
# indices_of_hsc_splashes_in_sp = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18]


high_speed_camera_df = drop_outliers(high_speed_camera_df, 'v')

self = EnergyStats(high_speed_camera_df, sticky_paper_df, number_of_q_energy)



print(f"{number_of_q_energy} {self.quantized_std:.4f} {self.quantized_min:.4f} {self.quantized_max:.4f} {self.min_of_one_part} {self.max_of_one_part}")
