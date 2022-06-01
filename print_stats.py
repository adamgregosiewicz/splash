import numpy as np
import pandas as pd
import os
import sys
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def cdf_from_df_df(distribution_function_df):
    """
    In:
        distribution_function_df: probability distribution function DataFrame(arguments, probabilities)
    Out:
        cumulative_distribution_function_df: DataFrame(arguments, cumulative probabilities)
    """
    distribution_function_df.iloc[:, 1] = np.cumsum(distribution_function_df.iloc[:, 1])
    return distribution_function_df


def quantiles_list(cdf_df, num_of_intervals):
    """
    In:
        cdf_df - CDF DatFrame(arguments, probabilities)
    Out:
        quantiles: list of arguments i's such that
                   cdf(i) < k / num_of_intervals <= cdf(i+1) for k = 1,...,num_of_intervals
        probabilities: list of respective values of the cdf
    """
    step = 1.0 / num_of_intervals
    quantiles = []
    probabilities = [step] * num_of_intervals

    threshold = step
    for arg, cdf in zip(cdf_df.iloc[:, 0], cdf_df.iloc[:, 1]):
        if cdf >= threshold:
            quantiles += [arg]
            threshold += step
    
    # add the last quantile if needed (the last quantile is the maximum argument of the cdf)
    if len(quantiles) < num_of_intervals:
        quantiles += [cdf_df.iloc[-1, 0]]

    return quantiles, probabilities


def idx_of_interval(quantiles, single_empirical):
    """
    In:
        quantiles: list of quantiles
        single_empirical: single empirical data
    Out:
        Index of the first quantile that is as least single_empirical
    """
    for idx, i in enumerate(quantiles):
        if single_empirical <= i:
            return idx
    return len(quantiles) - 1


def intervals_cardinality(quantiles, empirical):
    """
    In:
        quantiles: list of quantiles
        empirical: list of empirical data
    Out:
        list of cardinalities of each interval determined by quantiles for empirical data
    """
    intervals_card = [0] * len(quantiles)
    for e in empirical:
        intervals_card[idx_of_interval(quantiles, e)] += 1
    
    return intervals_card


def uniform_distribution_df(min_value, max_value):
    """
    Return the distribution function (DataFrame(arguments, probabilities)) of the discrete uniform distribution
     on the interval [min_value, max_value].
    """
    length = max_value - min_value + 1
    ud_cdf_df = pd.DataFrame()
    ud_cdf_df['no'] = range(min_value, max_value + 1)
    ud_cdf_df['prob'] = [1 / length] * length
    return ud_cdf_df


def chi_square(empirical, cumulative_distribution_function_df, num_of_intervals):
    """
    Return chi-square statistics.

    In:
        empirical: list of empirical data
        cumulative_distribution_function_df: theoretical CDF (DataFrame(arguments, probabilities))
        num_of_intervals: number of intervals for chi-square test
    Out:
        chi-square statistics and p-value
    """
    quantiles, probabilities = quantiles_list(cumulative_distribution_function_df, num_of_intervals)
    cardinalities = intervals_cardinality(quantiles, empirical)
    cardinalities_expected = np.multiply(probabilities, len(empirical))
    return stats.chisquare(cardinalities, cardinalities_expected)


def cdf_for_plot_from_sequence(values):
    """
    In:
        values: sorted list containing values
    Out:
        cdf_support: values of support at both points of jump discontinuities
        cdf_values: values of ECDF at both points of jump discontinuities
        CDF's true value at a jump discontinuity is the higher one
    """
    prev = values[0]
    n = len(values)
    cdf_support = [prev]
    cdf_values = [0.]

    for i in range(1, n):
        value = values[i]
        if value != prev:
            probability = i / n
            cdf_support += [prev, value]
            cdf_values += [probability, probability]
            prev = value

    cdf_support.append(prev)
    cdf_values.append(1.)

    return cdf_support, cdf_values


def cdf_for_plot_from_cdf(args, values):
    """
    In:
        args: sorted list containing arguments
        values: sorted list containing probabilities (values of CDF)
    Out:
        cdf_support: values of support at both points of jump discontinuities
        cdf_values: values of ECDF at both points of jump discontinuities
        CDF's true value at a jump discontinuity is the higher one
    Example:
        [3, 5], [0.2, 1.0] -> [3, 3, 5, 5], [0, 0.2, 0.2, 1]
    """
    cdf_support = []
    cdf_values = [0.0]
    for i in range(len(args)):
        cdf_support += [args[i], args[i]]
        cdf_values += [values[i], values[i]]

    cdf_values.pop()

    return cdf_support, cdf_values


def cdf_for_plot(values, args=None):
    """
    Returns coordinates of points for plotting CDF.
    In:
        values: if args = None, then this is a list of empirical data.
                if args != None, then this is a list of values of CDF.
        args: list of arguments of CDF
    Out:
        cdf_support: support of CDF
        cdf_values: values of CDF
    Example:
        [3, 5], [0.2, 1.0] -> [3, 3, 5, 5], [0, 0.2, 0.2, 1]
    """
    if args is None:
        cdf_support, cdf_values = cdf_for_plot_from_sequence(values)
    else:
        cdf_support, cdf_values = cdf_for_plot_from_cdf(args, values)
    
    return cdf_support, cdf_values


# read data from files
empirical_filename = sys.argv[1]  # empirical data
model_filename = sys.argv[2]  # model data (distribution function)
number_of_intervals = int(sys.argv[3])  # number of intervals for chi-square test

sticky_paper_df = pd.read_csv(empirical_filename)
model_df = pd.read_csv(model_filename)
model_cdf_df = cdf_from_df_df(model_df)

# sticky paper
# count splashes from 0
sticky_paper_df.iloc[:, 0] -= 1
sticky_paper_card = [len(sticky_paper_df[sticky_paper_df.iloc[:, 0] == i]) for i in range(48)]
sticky_paper_card_sorted = sorted(sticky_paper_card)

print(f'{number_of_intervals} intervals')
print(f'Expected cardinalities: '
      f'{np.multiply(quantiles_list(model_df, number_of_intervals)[1], len(sticky_paper_card_sorted))}')
print(f'IP: {chi_square(sticky_paper_card_sorted, model_df, number_of_intervals)}')
print(f'Intervals cardinalities: '
      f'{intervals_cardinality(quantiles_list(model_df, number_of_intervals)[0], sticky_paper_card_sorted)}')

uniform_df = cdf_from_df_df(uniform_distribution_df(34, 81))
print(f'U: {chi_square(sticky_paper_card_sorted, uniform_df, number_of_intervals)}')
print(f'Intervals cardinalities: '
      f'{intervals_cardinality(quantiles_list(uniform_df, number_of_intervals)[0], sticky_paper_card_sorted)}')

# Print graphs

# MODEL

x_model = list(model_df.iloc[:, 0])
y_model = list(model_df.iloc[:, 1])

x_model, y_model = cdf_for_plot(y_model, x_model)

# extend CDF to the left
shift = 5
x_min = x_model[0] - shift
x_max = x_model[-1] - shift
x_model = [x_min] + x_model
y_model = [0.0] + y_model

# STICKY PAPER

x_sticky_paper, y_sticky_paper = cdf_for_plot(sticky_paper_card_sorted)

# extend CDF to the left and to the right
x_sticky_paper = [x_min] + x_sticky_paper + [x_max]
y_sticky_paper = [0.0] + y_sticky_paper + [1.0]

# PLOT CDFs

fig, ax = plt.subplots()

ax.plot(x_sticky_paper, y_sticky_paper, color='black', label='Experiment')
ax.plot(x_model, y_model, color='red', label='Model')

ax.set_yticks([0.0, 1.0], minor=False)
ax.set_yticks([0.2, 0.4, 0.6, 0.8], minor=True)
ax.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.grid(True, which='major')
ax.legend(loc='right')
ax.set_xlabel('Number of splashed beads')
ax.set_ylabel('Cumulative distribution function')

plt.xlim(x_min, x_max)
plt.savefig(os.path.splitext(model_filename)[0] + '.svg', format='svg', dpi=1200)


x_model = list(uniform_df.iloc[:, 0])
y_model = list(uniform_df.iloc[:, 1])

x_model, y_model = cdf_for_plot(y_model, x_model)

# extend CDF to the left
#shift = 20
#x_min = x_model[0] - shift
#x_max = x_model[-1] + shift
x_model = [x_min] + x_model
y_model = [0.0] + y_model
x_model = x_model + [x_max]
y_model = y_model + [1.0]

fig, ax = plt.subplots()

ax.plot(x_sticky_paper, y_sticky_paper, color='black', label='Experiment')
ax.plot(x_model, y_model, color='red', label='Model')

ax.set_yticks([0.0, 1.0], minor=False)
ax.set_yticks([0.2, 0.4, 0.6, 0.8], minor=True)
ax.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.grid(True, which='major')
ax.legend(loc='right')
ax.set_xlabel('Number of splashed beads')
ax.set_ylabel('Cumulative distribution function')

plt.xlim(x_min, x_max)
plt.savefig(os.path.splitext(model_filename)[0] + '_uniform.svg', format='svg', dpi=1200)