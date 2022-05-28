import numpy as np
import pandas as pd
import random
import sys
import math
import matplotlib.pyplot as plt
import scipy.stats as stats

def find_outliers_IQR(df):
    q1=df.quantile(0.20)
    q3=df.quantile(0.80)
    IQR=q3-q1
    
    outliers = df[((df<(q1-1.5*IQR)) | (df>(q3+1.5*IQR)))]

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

#sticky_paper_len = [len(sticky_paper_df[sticky_paper_df['no'] == i]) for i in range(48)]
#sticky_paper_len = [sticky_paper_len[i] for i in no_of_hsc]
#print(sticky_paper_len)

#high_speed_cam_len = [len(high_speed_cam_df[high_speed_cam_df['no'] == i]) for i in range(16)]

#print(sum(sticky_paper_len)/sum(high_speed_cam_len))

#sys.exit()


# outliers handling

# print(find_outliers_IQR(high_speed_cam_df['v']).index.tolist())
# print(high_speed_cam_df['v'].max())
high_speed_cam_df = high_speed_cam_df.drop(find_outliers_IQR(high_speed_cam_df['v']).index.tolist())
# print(high_speed_cam_df['v'].max())



energy_sum_df = pd.DataFrame([high_speed_cam_df[high_speed_cam_df['no'] == i]['e'].sum() for i in range(16)], columns = ['energy_sum'])
# print(energy_sum_df)

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

#numMin = math.ceil(energy_scaled_min / energy_p_max)
#numMax = math.floor(energy_scaled_max  / energy_p_min)

numMin = math.ceil(energy_q_min / kmax)
numMax = math.floor(energy_q_max  / kmin)

print(f"kmin = {kmin}, kmax = {kmax}")
print(f"parts min = {numMin}, max = {numMax}\n")

print(f"{Q} {energy_q_std:.4f} {energy_q_min:.4f} {energy_q_max:.4f} {kmin} {kmax} {numMin} {numMax}")

sys.exit()

sticky_paper_card_sorted = sorted(sticky_paper_len)
sticky_paper_card_sorted_prob = np.cumsum([i / sum(sticky_paper_len) for i in sticky_paper_card_sorted])
#print(sticky_paper_card_sorted)

model_df = pd.read_csv("data303.csv")
x = list(model_df['n'])
y = model_df['prob']
model_prob = model_df['prob'].to_list()
y = y.cumsum()
y = y.to_list()
# print(y)


mu = 200
sigma = 25
n_bins = 50
# print(x)
#fig, ax = plt.subplots(figsize=(8, 4))
fig, ax = plt.subplots()

# plot the cumulative histogram
bins = sticky_paper_card_sorted + [np.inf]
#n, bins, patches = ax.hist(sticky_paper_card_sorted, histtype='step',
#                           cumulative=True, label='Empirical', bins = bins)

def ecdf4plot(seq, assumeSorted = False):
    """
    In:
    seq - sorted-able object containing values
    assumeSorted - specifies whether seq is sorted or not
    Out:
    0. values of support at both points of jump discontinuities
    1. values of ECDF at both points of jump discontinuities
       ECDF's true value at a jump discontinuity is the higher one    """
    if not assumeSorted:
        seq = sorted(seq)
    prev = seq[0]
    n = len(seq)
    support = [prev]
    ECDF = [0.]
    for i in range(1, n):
        seqi = seq[i]
        if seqi != prev:
            preP = i/n
            support.append(prev)
            ECDF.append(preP)
            support.append(seqi)
            ECDF.append(preP)
            prev = seqi
    support.append(prev)
    ECDF.append(1.)
    return support, ECDF

x_sp, y_sp = ecdf4plot(sticky_paper_card_sorted)
plt.plot(x_sp, y_sp, label='Empirical')

#ax.plot(sticky_paper_card_sorted, sticky_paper_card_sorted_prob, 'ro')

ax.plot(x, y, 'r--', linewidth=1.5, label='Model')

# tidy up the figure
ax.grid(True)
plt.xlim(25, 105)
ax.legend(loc='right')
ax.set_title('Cumulative step histograms')
ax.set_xlabel('Annual rainfall (mm)')
ax.set_ylabel('Likelihood of occurrence')

plt.show()

#stats.chisquare(x_sp, x)

print(x[0],x[-1])

# sp_obs = [0] * (numMax - numMin + 1)
# for i in sticky_paper_card_sorted:
#     sp_obs[i-numMin] += 1
# #print(sp_obs)

# exp = [48 * p for p in model_prob]
# print(len(exp))

# chi2 = 0
# for i in range(numMax - numMin + 1):
#     chi2 += (sp_obs[i] - exp[i])**2 / exp[i]

# print(chi2)

# print(stats.chisquare(sp_obs, exp))

num_of_intervals = 9
p_for_intervals = 1 / num_of_intervals

# right ends of intervals
intervals = []

sum = 0
step = p_for_intervals
for i, p in enumerate(y):
    # print(i, p)
    if p > step:
        intervals += [i]
        print(p)
        step += p_for_intervals

print(intervals)

# theoretical distribution on [min,max]
def intervals_right_ends(cumulative_probabilities, num_of_intervals):
    step = 1 / num_of_intervals
    intervals = []
    probabilities = [step] * num_of_intervals

    threshold = step
    for i, p in enumerate(cumulative_probabilities):
        if p > threshold:
            intervals += [i]
            threshold += step
    
    if len(intervals) < num_of_intervals:
        intervals += [len(cumulative_probabilities) - 1]

    return intervals, probabilities

def uniform_distribution(min, max):
    length = max - min + 1
    return np.cumsum([1 / length] * length)

num_of_intervals = 5

# print(generate_uniform(34,81))
intervals, probabilities = intervals_right_ends(uniform_distribution(34,81),num_of_intervals)
print(intervals)
print(probabilities)

def idx_of_interval(min, intervals, empirical):
    for idx, i in enumerate(intervals):
        if empirical - min <= i:
            return idx
    return len(intervals) - 1

def intervals_cardinality(min, intervals, empirical):
    intervals_card = [0] * len(intervals)
    for e in empirical:
        intervals_card[idx_of_interval(min, intervals, e)] += 1
    
    return intervals_card

print(intervals_cardinality(34, intervals, [34, 36, 38, 32, 56, 61, 63, 70, 74, 80, 90]))

print(stats.chisquare(intervals_cardinality(34, intervals, [34, 36, 38, 32, 56, 61, 63, 70, 74, 80, 90]), np.multiply(probabilities, 11)))

cards = intervals_cardinality(34, intervals, sticky_paper_card_sorted)
print(sticky_paper_card_sorted)
print(cards)
cards_expected = np.multiply(probabilities, 48)
print(cards_expected)
print(stats.chisquare(cards, cards_expected))

intervals, probabilities = intervals_right_ends(y, num_of_intervals)
print(intervals)
print(probabilities)
cards = intervals_cardinality(20, intervals, sticky_paper_card_sorted)
print(sticky_paper_card_sorted)
print(cards)
cards_expected = np.multiply(probabilities,48)
print(cards_expected)
print(stats.chisquare(cards, cards_expected))