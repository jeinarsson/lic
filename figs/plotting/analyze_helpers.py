import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import time
import h5py




def plot_mean_and_deviation(ax, times, ydata, min_value = 1e-10):

	meanvalues = np.mean(ydata, axis=1)
	deviations = np.std(ydata, axis=1)

	try:
		ax.semilogy(times , meanvalues,label=label)
		lower_bound = np.clip(meanvalues-deviations,min_value,np.inf)
		upper_bound = meanvalues+deviations
		ax.fill_between(times,lower_bound,upper_bound,color=(.8,.8,.8,.5))
	except ValueError:
		ax.set_yscale('linear')
		ax.plot(times,meanvalues)



def plot_mean_and_deviation_nonlog(ax, times, ydata, min_value = 1e-10):

	meanvalues = np.mean(ydata, axis=1)
	deviations = np.std(ydata, axis=1)

	ax.plot(times , meanvalues)
	#lower_bound = np.clip(meanvalues-deviations,min_value,np.inf)
	#upper_bound = meanvalues+deviations
	ax.fill_between(times,meanvalues-deviations,meanvalues+deviations,color=(.8,.8,.8,.5))


def plot_mean_and_deviation_markers(ax, times, ydata, min_value = 1e-10,label=''):

	meanvalues = np.mean(ydata, axis=1)
	deviations = np.std(ydata, axis=1)

	#ax.plot(times , meanvalues, marker='o', color='white')
	ax.errorbar(times, meanvalues, yerr=deviations, ls='none',marker='o', label=label)
	

def offset_axes_color_cycle(axes, n):

	for ax in axes:
		[ax._get_lines.color_cycle.next() for k in xrange(n)]


def linear_to_square(ij):
	#return (np.ceil(ij/3)+1, 1+ij-3*(np.ceil(ij/3)))
	return (ij/3+1, 1+np.mod(ij,3))

def square_to_linear(ij):
	return (ij[0]-1)*3+ij[1]-1

def transpose_index(ij):
	return square_to_linear(linear_to_square(ij)[::-1])

def index_name(ij):
	tup = linear_to_square(ij)
	return "{}{}".format(tup[0],tup[1])

def get_Acorr2(output_file, ij, kl, time_slices=slice(None,None)):
	
	step = 1
	if ij > kl:
		temp = kl
		kl = ij
		ij = temp
		step = -1

	dataset = output_file['/A_{}_{}'.format(ij,kl)]
	data = dataset[time_slices,:]
	return data[::step,:]

def get_Acorr3_onetime(output_file, ij, kl, mn, time_slices=slice(None,None)):
	
	if ij > kl:
		temp = kl
		kl = ij
		ij = temp

	dataset = output_file['/A_{}_{}_{}'.format(ij,kl,mn)]
	data = dataset[time_slices,:]
	return data

def flatten(l):
	return [i for sub in l for i in sub]