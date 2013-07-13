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
import matplotlib.image as mpimg
def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')


	fig1 = plt.figure(figsize=(12./2.5, 7./2.5))
	
	gs = gridspec.GridSpec(1,3, hspace=.5,bottom=0.1,left=0.1,right=0.9,top=0.9)


	ax1 = plt.subplot(gs[0, 0])
	ax2 = plt.subplot(gs[0, 1])
	ax3 = plt.subplot(gs[0, 2])

	axes = [ax1, ax2, ax3]

	for k in xrange(3):

		img=mpimg.imread('orbit{}.png'.format(k+1))
		axes[k].imshow(img)
		#axes[k].axis('off')
		
	#for ax in axes:
	#ax.text(.03, .22, r'(a)')


	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/projection.pdf')
	pp.savefig(fig1, dpi=600)
	pp.close()
	#ax.set_frame_on(False)
	#plt.subplots_adjust(left=0.2, right=.96, top=.96, bottom=0.2)

	#fig1.savefig('output/poincare{}.png'.format(output),dpi=300)


	#plt.show()
	print "Done!"



if __name__ == '__main__':
	main()