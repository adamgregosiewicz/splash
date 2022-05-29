import numpy as np
import pandas as pd
import random
import sys
import math
import matplotlib.pyplot as plt
import scipy.stats as stats

def find_outliers_IQR(df):
    """
    Returns outliers in DataFrame based on interquartile range (IQR)
    https://careerfoundry.com/en/blog/data-analytics/how-to-find-outliers/
    """
    q1 = df.quantile(0.20)
    q3 = df.quantile(0.80)
    iqr = q3 - q1
    
    outliers = df[((df < (q1 - 1.5 * iqr)) | (df > (q3 + 1.5 * iqr)))]

    return outliers

high_speed_cam_df = pd.read_csv("high_speed_camera_1529.csv")
sticky_paper_df = pd.read_csv("sp_1529.csv")

#high_speed_cam_df = pd.read_csv("hsc_15.csv")
#sticky_paper_df = pd.read_csv("sp_15.csv")

# count splashes from 0
high_speed_cam_df['no'] -= 1
sticky_paper_df['no'] -= 1

# indicies of sticky paper registered by hsc
no_of_hsc = [1,2,3,4,5,7,8,9,10,11,12,13,14,16,17,18]


# outliers handling
high_speed_cam_df = high_speed_cam_df.drop(find_outliers_IQR(high_speed_cam_df['v']).index.tolist())



energy_sum_df = pd.DataFrame([high_speed_cam_df[high_speed_cam_df['no'] == i]['e'].sum() for i in range(16)], columns = ['energy_sum'])

energy_mean = energy_sum_df.mean()[0]
energy_std = energy_sum_df.std()[0]
energy_min = energy_sum_df.min()[0]
energy_max = energy_sum_df.max()[0]

energy_p_min = high_speed_cam_df['e'].min()
energy_p_max = high_speed_cam_df['e'].max()

print(f"Energy mean = {energy_mean}, std = {energy_std}")
print(f"Energy min = {energy_min}, max = {energy_max}\n")

print(f"Particle energy min = {energy_p_min}, max = {energy_p_max}\n")

# cardinality of each splash
high_speed_cam_len = [len(high_speed_cam_df[high_speed_cam_df['no'] == i]) for i in range(16)]
sticky_paper_len = [len(sticky_paper_df[sticky_paper_df['no'] == i]) for i in range(48)]

high_speed_cam_len_sum = sum(high_speed_cam_len)
sp_hsc_card = pd.DataFrame([sum(random.sample(sticky_paper_len, 16))/sum(high_speed_cam_len) for i in range(1000)], columns = ['z'])

scaling_mean = sp_hsc_card.mean()[0]
scaling_std = sp_hsc_card.std()[0]

scaling_mean = 2.95

energy_scaled_mean = energy_mean * scaling_mean
energy_scaled_std = energy_std * scaling_mean
energy_scaled_min = energy_min * scaling_mean
energy_scaled_max = energy_max * scaling_mean

print(f"Scaling mean = {scaling_mean}, std = {scaling_std}, coeff of var = {scaling_std/scaling_mean}")
print(f"Scaled energy mean = {energy_scaled_mean}, std = {energy_scaled_std}")
print(f"Scaled energy min = {energy_scaled_min}, std = {energy_scaled_max}\n")

Q = int(sys.argv[1])
q = energy_scaled_mean / Q

energy_q_mean = energy_scaled_mean / q
energy_q_std = energy_scaled_std / q
energy_q_min = energy_scaled_min / q
energy_q_max = energy_scaled_max / q

print(f"Quant energy mean = {energy_q_mean}, std = {energy_q_std}")
print(f"Quant energy min = {energy_q_min}, max = {energy_q_max}\n")

kmin = math.floor(energy_p_min / q)
kmax = math.ceil(energy_p_max / q)

numMin = math.ceil(energy_q_min / kmax)
numMax = math.floor(energy_q_max  / kmin)

print(f"kmin = {kmin}, kmax = {kmax}")
print(f"parts min = {numMin}, max = {numMax}\n")

print(f"{Q} {energy_q_std:.4f} {energy_q_min:.4f} {energy_q_max:.4f} {kmin} {kmax} {numMin} {numMax}")

