import numpy as np
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.backends.backend_pdf
import time
import scipy.io
from os.path import exists
from os import mkdir
import h5py

import matplotlib.gridspec as gridspec
from itertools import cycle
from matplotlib.patches import Ellipse

def main():

	# poincare A
	# series_files = ["series10-1-1-{}.mat".format(k) for k in xrange(1,4)]
	# output = "A"



	# # poincare B
	# series_files = ["series10-1.1-1-{}.mat".format(k) for k in xrange(1,4)]
	# filename = "poincaredata/1o10.mat"
	# output = "B"

	# # poincare C
	series_files = ["series10-1.3-1-{}.mat".format(k) for k in xrange(1,4)]
	filename = "poincaredata/1o30.mat"
	output = "C"



	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')


	fig1 = plt.figure(figsize=(12./2.5, 16./2.5))
	
	gs = gridspec.GridSpec(6+6+1,1, hspace=.5,bottom=0.1,left=0.18,right=0.94,top=0.95)

	ax_poincare = plt.subplot(gs[0:6, 0])
	ax_series = [plt.subplot(gs[(7+2*k):(9+2*k), 0]) for k in xrange(3)]
	
	# plot poincare map
	ax = ax_poincare

	if output == "A":

		numlines = 51
		for yk in xrange(numlines):
			y=-1.0+2.0*float(yk)/(numlines-1.)
			ax.plot([-.5,.5],[y,y],linewidth=0.4,color='#444444')		

	else:

		datafile = h5py.File(filename, 'r')
		print datafile['/'].keys()
		
		
		thetapoints = datafile['/thetapoints'][:,:].flatten()
		psipoints = datafile['/psipoints'][:,:].flatten()

		
		ax.scatter(psipoints*1./np.pi,np.cos(thetapoints),s=.05, facecolor='0.5', lw = 0)
		ax.scatter(-psipoints*1./np.pi,-np.cos(thetapoints),s=.05, facecolor='0.5', lw = 0)

	ax.set_xlim(-.5,.5)
	ax.set_ylim(-1,1)
	ax.set_xlabel(r'$\psi$',labelpad=-1)
	ax.set_ylabel(r'$n_z$',rotation=0)
	ax.set_xticks([-.5,0,.5])
	ax.set_yticks([-1,-.5,0,.5,1])
	ax.set_xticklabels([r'$-\pi/2$',r'$0$',r'$\pi/2$'])

	pos = ax.get_position().get_points()
	print pos
	aspect_ratio = 2.1
	pointsize = .015
	
	colors=['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33']
	for (series_k, filename) in enumerate(series_files):

		datafile = h5py.File("poincaredata/series/{}".format(filename), 'r')
		print datafile['/'].keys()
		psis = datafile['/poincare_psi'][0,:]
		zs = datafile['/poincare_nz'][0,:]
		ts = datafile['/t'][:].transpose()
		nzs = datafile['/n_z'][:].transpose()
		for k in range(psis.shape[0]):
			x = psis[k]*1./np.pi;
			y = zs[k];
			ax.add_artist(Ellipse(xy=(x,y), width=pointsize, height=pointsize*aspect_ratio,facecolor=colors[series_k], linewidth=0.2))	


		#ax_series[series_k].plot(ts, nzs)
		ax_series[series_k].plot(ts, nzs, color=colors[series_k])
		ax_series[series_k].set_ylim(-1.1,1.1)
		ax_series[series_k].set_xlim(0,500.)
		ax_series[series_k].set_ylabel('$n_z$',rotation=0)
		ax_series[series_k].set_xticklabels([])
		ax_series[series_k].set_yticks([-1,0,1])
		
	ax_series[2].set_xlabel('$t$')
		
	#ax.add_artist(Ellipse(xy=(0,.35), width=pointsize, height=pointsize*aspect_ratio,facecolor='#800000', linewidth=0.5))	
	#ax.text(.03, .22, r'(a)')


	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	#pp = matplotlib.backends.backend_pdf.PdfPages('output/poincare.pdf')
	#pp.savefig(fig1)
	#pp.close()
	#ax.set_frame_on(False)
	#plt.subplots_adjust(left=0.2, right=.96, top=.96, bottom=0.2)

	fig1.savefig('output/poincare{}.png'.format(output),dpi=300)


	#plt.show()
	print "Done!"



if __name__ == '__main__':
	main()