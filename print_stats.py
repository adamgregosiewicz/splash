import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats

def cdf_df(model_df):
    model_df.iloc[:,1] = np.cumsum(model_df.iloc[:,1])
    return model_df

def quantiles_list(cdf_df, num_of_intervals):
    """
    In:
        cdf_df - CDF DatFrame (args and values)
    Out:
        0. list of arguments i's such that cdf(i) < 1 / num_of_intevals <= cdf(i+1)
        1. [cdf(i)]
    """
    step = 1 / num_of_intervals
    quantiles = []
    probabilities = [step] * num_of_intervals

    threshold = step
    for arg, cdf in zip(cdf_df.iloc[:,0], cdf_df.iloc[:,1]):
        if cdf > threshold:
            quantiles += [arg]
            threshold += step
    
    if len(quantiles) < num_of_intervals:
        quantiles += [len(cdf_df) - 1]

    return quantiles, probabilities

def idx_of_interval(quantiles, single_empirical):
    for idx, i in enumerate(quantiles):
        if single_empirical <= i:
            return idx
    return len(quantiles) - 1

def intervals_cardinality(quantiles, empirical):
    intervals_card = [0] * len(quantiles)
    for e in empirical:
        intervals_card[idx_of_interval(quantiles, e)] += 1
    
    return intervals_card


def uniform_distribution_df(min, max):
    length = max - min + 1
    ud_cdf_df = pd.DataFrame()
    ud_cdf_df['no'] = range(min,max+1)
    ud_cdf_df['prob'] = [1 / length] * length
    return ud_cdf_df


def chisquare(empirical, model_df, num_of_intervals):
    quantiles, probabilities = quantiles_list(cdf_df(model_df), num_of_intervals)
    cardinalities = intervals_cardinality(quantiles, empirical)
    cardinalities_expected = np.multiply(probabilities,len(empirical))
    return stats.chisquare(cardinalities, cardinalities_expected)


empirical_filename = sys.argv[1]
model_filename = sys.argv[2]
num_of_intervals = int(sys.argv[3])

sticky_paper_df = pd.read_csv(empirical_filename)
model_df = pd.read_csv(model_filename)


# STICKY PAPER
# count splashes from 0
sticky_paper_df['no'] -= 1

sticky_paper_card = [len(sticky_paper_df[sticky_paper_df['no'] == i]) for i in range(48)]

sticky_paper_card_sorted = sorted(sticky_paper_card)

print(f'{num_of_intervals} intervals')
print(f'IP: {chisquare(sticky_paper_card_sorted, model_df, num_of_intervals)}')
print(f'U: {chisquare(sticky_paper_card_sorted, uniform_distribution_df(34,81), num_of_intervals)}')


sys.exit()

# MODEL
x = list(model_df['no'])
y = model_df['prob']
y = y.cumsum()
y = y.to_list()
# print(y)

# print(generate_uniform(34,81))
intervals, probabilities = quantiles_list(uniform_distribution_cdf_df(34,81),num_of_intervals)
print(intervals)
print(probabilities)
cards = intervals_cardinality(34, intervals, sticky_paper_card_sorted)
cards_expected = np.multiply(probabilities, 48)
# print(cards)
# print(cards_expected)
print(stats.chisquare(cards, cards_expected))

intervals, probabilities = quantiles_list(y, num_of_intervals)
# print(intervals)
# print(probabilities)
cards = intervals_cardinality(x[0], intervals, sticky_paper_card_sorted)
cards_expected = np.multiply(probabilities,48)
# print(cards)
# print(cards_expected)
print(stats.chisquare(cards, cards_expected))












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
plt.xlim(15, 115)
ax.legend(loc='right')
ax.set_title('Cumulative step histograms')
ax.set_xlabel('Annual rainfall (mm)')
ax.set_ylabel('Likelihood of occurrence')

plt.show()

#stats.chisquare(x_sp, x)

#print(x[0],x[-1])

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

