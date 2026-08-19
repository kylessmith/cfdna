import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from lmfit import models
from scipy import stats


def peak_decov(x, y, peaks_cen, max_scale=None):
	'''
	x: array of sizes.
	y: array of numbers or frequencies.
	peaks_cen: array of roughtly locations.
	Without scale constrains.'''
	x, y = np.array(x), np.array(y)/sum(y)*100
	model_ar, model, params = [], None, None
	for peak in peaks_cen:
		mod = models.LorentzianModel(prefix='peak{}_'.format(peak))
		mod.set_param_hint('center', value=peak,)
		mod.set_param_hint('amplitude', min=0)
		mod.set_param_hint('sigma', max=max_scale)
		par = mod.make_params()
		model_ar.append(mod)
		if model is None:
			model = mod
			params = par
		else:
			model = model + mod
			params.update(par)
	result = model.fit(y, params, x=x)

	result_par = {}
	for peak in peaks_cen:
		loc = result.values['peak{}_center'.format(peak)]  # loc
		scale = result.values['peak{}_sigma'.format(peak)]  # width
		amp = result.values['peak{}_amplitude'.format(peak)]  # area
		result_par.setdefault(loc, (scale, amp))

	final_res = []
	for loc, (scale, amp) in sorted(result_par.items()):
		final_res.append((loc, scale, amp))
	return final_res


def peak_decov_constrain(fig_x, fig_y, peaks_cen, max_scale=None):
	'''With scale constrains.'''
	fig_x, fig_y = np.array(fig_x), np.array(fig_y)/sum(fig_y)*100
	model_ar, model, params = [], None, None
	for peak in peaks_cen:
		mod = models.LorentzianModel(prefix='peak{}_'.format(peak))
		mod.set_param_hint('center', value=peak)
		mod.set_param_hint('amplitude', min=0)

		if peak == 70:
			mod.set_param_hint('sigma', expr='peak60_sigma')
		elif peak == 90:
			mod.set_param_hint('sigma', expr='peak80_sigma')
		elif peak == 110:
			mod.set_param_hint('sigma', expr='peak100_sigma')
		elif peak == 130:
			mod.set_param_hint('sigma', expr='peak120_sigma')
		elif peak == 150:
			mod.set_param_hint('sigma', expr='peak140_sigma')
		elif peak >= 190:
			mod.set_param_hint('sigma', expr='peak180_sigma')
		else:
			mod.set_param_hint('sigma', max=max_scale)

		par = mod.make_params()
		model_ar.append(mod)
		if model is None:
			model = mod
			params = par
		else:
			model = model + mod
			params.update(par)
	result = model.fit(fig_y, params, x=fig_x)

	result_par = {}
	for peak in peaks_cen:
		loc = result.values['peak{}_center'.format(peak)]  # loc
		scale = result.values['peak{}_sigma'.format(peak)]  # width
		amp = result.values['peak{}_amplitude'.format(peak)]  # area
		result_par.setdefault(loc, (scale, amp))

	final_res = []
	for loc, (scale, amp) in sorted(result_par.items()):
		final_res.append((loc, scale, amp))
	return final_res


def color_iter(num, color_family, min_frac, max_frac):
	color_family = eval("cm.{}".format(color_family))
	return iter(color_family(np.linspace(min_frac, max_frac, num)))


def plot_peaks(x, y, pars, fig_file='size_decov.pdf', title='', ylim_max=None):
	x, y = np.array(x), np.array(y)/sum(y)*100
	model_ar, model, params = [], None, None
	x_high = np.linspace(min(x), max(x), num=1000)
	fig, axs = plt.subplots(
		2, 1,
		sharex=True,
		gridspec_kw={'height_ratios': [20, 1]}, 
		figsize=(6, 6))

	axs[0].plot(
		x,
		y,
		color='silver',
		linestyle='-',
		label='Size profile',
		linewidth=4,
		alpha=1)

	best_fit = np.array([0. for x in x])
	colors = color_iter(len(pars), 'spring', 0., 0.85)
	cols = list(colors)

	for par, col in zip(pars, cols):
		loc, scale, amp = par[:3]
		dist = stats.cauchy(loc=loc, scale=scale)
		fit = dist.pdf(x) * amp
		best_fit += fit
		y_high = dist.pdf(x_high)*amp 
		axs[0].plot(
			x_high,
			y_high,
			color=col,
			linestyle='-',
			linewidth=1.5)

	axs[0].plot(
		x,
		best_fit,
		color='k',
		linestyle='--',
		label='Best fit',
		linewidth=1)

	if ylim_max is None:
		axs[0].set_yticks([_ for _ in range(int(round(max(y)))+1)])
		axs[0].set_yticklabels([f'{_}' for _ in range(int(round(max(y)))+1)])
	else:
		axs[0].set_yticks([_ for _ in range(ylim_max+1)])
		axs[0].set_yticklabels([f'{_}' for _ in range(ylim_max+1)])
	axs[0].tick_params(axis='y', labelsize=15)
	axs[0].set_ylim(ymin=0)
	axs[0].set_ylabel('Frequency (%)', fontsize=20)
	axs[0].tick_params(axis="y", labelsize=15)

	for par, col in zip(pars, cols):
		loc, scale, amp = par[:3]
		dist = stats.cauchy(loc=loc, scale=scale)

		axs[1].scatter(
			[loc],
			[0],
			color=col,
			s=15 if 158 < loc < 161 else 10,
			marker='^' if 158 < loc < 161 else 'o',
			facecolor='none' if 158 < loc < 161 else col,
			)  # center

		axs[1].hlines(
			-0.5,
			loc-scale/2,
			loc+scale/2,
			color=col,
			lw=1)
	axs[1].set_ylim(-0.8, 0.5)

	axs[1].set_yticks([])
	axs[1].set_xlabel('Size (bp)', fontsize=20)
	axs[1].tick_params(axis="x", labelsize=15)

	fig.suptitle(title, fontsize=20)
	plt.subplots_adjust(left=0.15, right=0.85, top=0.92, bottom=0.15, hspace=0)
	plt.savefig(fig_file)


def decompose_frag_lengths(frags, plot_name=None):
	region_ld = frags.length_dist()
	y = region_ld
	x = np.arange(len(region_ld)) + 1
	y = y[np.logical_and(x >= 50, x <= 250)]
	x = x[np.logical_and(x >= 50, x <= 250)]

	result_pars = peak_decov_constrain(
		x,
		y,
		[60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200],
		max_scale=8,
		)
	
	if plot_name is not None:
		plot_peaks(x, y, result_pars, fig_file=plot_name)

	return result_pars

