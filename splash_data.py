import numpy as np
import pandas as pd
import random
import sys
import math


# Just to remember: indices (indexing from 0) of sticky paper registered by hsc
# indices_of_hsc_splashes_in_sp = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18]


class Splashes:
    def __init__(self, argv):
        """
        Take CSV DataFrames for HSC and SP data and the number of quants of energy.
        Index splashes from 0 and not from 1.
        """
        self.high_speed_camera_df = pd.read_csv(argv[0])
        self.sticky_paper_df = pd.read_csv(argv[1])
        self.number_of_quants_of_energy = int(argv[2])
        self.__shift_splashes_numbers()

    def __shift_splashes_numbers(self):
        """
        Index splashes from 0 and not from 1.
        """
        self.high_speed_camera_df['no'] -= 1
        self.sticky_paper_df['no'] -= 1

    def __velocity_outliers_of_hsc(self, quantile: float):
        """
        Return (indices of) outliers in DataFrame based on inter-quartile range (IQR).
        https://careerfoundry.com/en/blog/data-analytics/how-to-find-outliers/
        """
        q1 = self.high_speed_camera_df['v'].quantile(quantile)
        q3 = self.high_speed_camera_df['v'].quantile(1 - quantile)
        iqr = q3 - q1

        outliers = self.high_speed_camera_df['v'][((self.high_speed_camera_df['v'] < (q1 - 1.5 * iqr))
                                                   | (self.high_speed_camera_df['v'] > (q3 + 1.5 * iqr)))]

        return outliers

    def drop_velocity_outliers_of_hsc(self, quantile: float):
        self.high_speed_camera_df = self.high_speed_camera_df.drop(self.__velocity_outliers_of_hsc(quantile).index.tolist())

    def mean_std_sp_to_hsc_cardinality(self, sample_count):
        """
        Return mean and std of the random quotient of cardinalities of SP and HSC splashes.

        Since there are only 16 splashes recorded by HSC we take a random sample of 16 splashes
        recorder by SP and divide the sum of beads in these 16 splashes by the sum of beads in
        all HSC splashes.
        """

        high_speed_camera_sum = sum(
            [len(self.high_speed_camera_df[self.high_speed_camera_df['no'] == i]) for i in range(16)])
        sticky_paper_cardinalities = [len(self.sticky_paper_df[self.sticky_paper_df['no'] == i]) for i in range(48)]
        sp_hsc_quotient = [sum(random.sample(sticky_paper_cardinalities, 16))
                           / high_speed_camera_sum for _ in range(sample_count)]

        return np.mean(sp_hsc_quotient), np.std(sp_hsc_quotient)


class EnergyStats:
    def __init__(self, splashes: Splashes):
        self.sum = [splashes.high_speed_camera_df[splashes.high_speed_camera_df['no'] == i]['e'].sum() for i in
                    range(16)]
        self.mean = np.mean(self.sum)
        self.std = np.std(self.sum)
        self.min = np.min(self.sum)
        self.max = np.max(self.sum)
        self.min_particle = np.min(splashes.high_speed_camera_df['e'])
        self.max_particle = np.max(splashes.high_speed_camera_df['e'])
        self.scaling_mean, self.scaling_std = splashes.mean_std_sp_to_hsc_cardinality(1000)
        self.scaled_mean = self.mean * self.scaling_mean
        self.scaled_std = self.std * self.scaling_mean
        self.scaled_min = self.min * self.scaling_mean
        self.scaled_max = self.max * self.scaling_mean
        self.quant_of_energy = self.scaled_mean / splashes.number_of_quants_of_energy
        self.quantized_mean = self.scaled_mean / self.quant_of_energy
        self.quantized_std = self.scaled_std / self.quant_of_energy
        self.quantized_min = self.scaled_min / self.quant_of_energy
        self.quantized_max = self.scaled_max / self.quant_of_energy
        self.min_of_one_part = math.floor(self.min_particle / self.quant_of_energy)
        self.max_of_one_part = math.ceil(self.max_particle / self.quant_of_energy)
        self.min_number_of_parts = math.ceil(self.quantized_min / self.max_of_one_part)
        self.max_number_of_parts = math.floor(self.quantized_max / self.min_of_one_part)

    def print_stats(self):
        print(f"Energy mean = {self.mean}, std = {self.std}")
        print(f"Energy min = {self.min}, max = {self.max}")
        print(f"Particle energy min = {self.min_particle}, max = {self.max_particle}")
        print(f"Scaling mean = {self.scaling_mean}, std = {self.scaling_std}, "
              f"coefficient of variance = {self.scaling_std / self.scaling_mean}")
        print(f"Scaled energy mean = {self.scaled_mean}, std = {self.scaled_std}")
        print(f"Scaled energy min = {self.scaled_min}, std = {self.scaled_max}")
        print(f"Quant energy mean = {self.quantized_mean}, std = {self.quantized_std}")
        print(f"Quant energy min = {self.quantized_min}, max = {self.quantized_max}")
        print(f"Min part size = {self.min_of_one_part}, max = {self.max_of_one_part}")
        print(f"Number of parts min = {self.min_number_of_parts}, max = {self.max_number_of_parts}")


def main():
    hsc_filename = sys.argv[1]
    sp_filename = sys.argv[2]
    number_of_quants_of_energy = int(sys.argv[3])

    splashes = Splashes([hsc_filename, sp_filename, number_of_quants_of_energy])

    splashes.drop_velocity_outliers_of_hsc(0.25)

    energy_stats = EnergyStats(splashes)

    print(f"{splashes.number_of_quants_of_energy} "
          f"{energy_stats.quantized_std:.4f} "
          f"{energy_stats.quantized_min:.4f} "
          f"{energy_stats.quantized_max:.4f} "
          f"{energy_stats.min_of_one_part} "
          f"{energy_stats.max_of_one_part}")

main()
