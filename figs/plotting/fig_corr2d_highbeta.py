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

def approximation6(l, betasqr, ts):

	part_two = 1 + betasqr*(-1)*(1+l**2)**(-2)*(3+2*l**2+3*l**4)*((-1)+np.exp(-ts)+ts);
	part_four = betasqr**2*(float(1)/2)*(1+l**2)**(-4)*(3+2*l**2+3*l**4)*((-47)+ \
				np.exp(-2*ts)*(13+(-18)*l**2+13*l**4)+14*ts+3*ts**2+2* \
				l**2*(51+((-22)+ts)*ts)+l**4*((-47)+ts*(14+3*ts))+ \
				2*np.exp(-ts)*(17+23*ts+(-2)*l**2*(21+19*ts)+l**4*(17+ \
				23*ts)));
	part_six = betasqr**3*(float(1)/18)*(1+l**2)**(-6)*(3+2*l**2+3*l**4)*( \
				np.exp(-3*ts)*((-1191)+2668*l**2+(-3146)*l**4+2668* \
				l**6+(-1191)*l**8)+(-9)*np.exp(-2*ts)*(693+759*ts+33* \
				l**8*(21+23*ts)+(-4)*l**2*(621+487*ts)+(-4)* \
				l**6*(621+487*ts)+2*l**4*(1759+1221*ts))+(-4)* \
				l**2*(10657+3*ts*((-805)+3*((-23)+ts)*ts))+(-4)* \
				l**6*(10657+3*ts*((-805)+3*((-23)+ts)*ts))+(-3)*(( \
				-6067)+3*ts*(391+3*ts*(17+ts)))+(-3)*l**8*((-6067)+3* \
				ts*(391+3*ts*(17+ts)))+(-2)*l**4*((-24523)+3*ts*(2143+ \
				ts*(27+11*ts)))+(-9)*np.exp(-ts)*(3*(399+ts*(610+239*ts))+3* \
				l**8*(399+ts*(610+239*ts))+(-4)*l**2*(489+ts*( \
				1198+329*ts))+(-4)*l**6*(489+ts*(1198+329*ts))+2* \
				l**4*(791+ts*(2898+631*ts))));

	return part_two + part_four + part_six
			
def approximation_exp(l, betasqr, ts):
	return np.exp(betasqr*(3*l**4 + 2*l**2 +3)/(l**2+1)**2 * (-ts));

def approximation_flowdecorrelation(l, betasqr, ts):
	return np.exp((-ts));

def approximation2_renorm(l, betasqr, ts):
	return np.exp(2.*betasqr*(1.-ts-np.exp(-ts)));


def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')

	fig1 = plt.figure(1, figsize=(18/2.5, 8/2.5))
	ax1 = fig1.add_subplot(1,2,1)
	ax2 = fig1.add_subplot(1,2,2)


	all_parameters = []
	all_parameters.append({
		'lambda': 1,
		'betasqr': .9
	})


	for parameters in all_parameters:
		datafile_name = "paper_all_2d/" + "correlations_l%(lambda)e_betasqr%(betasqr)e.hdf5" % parameters
		datafile = h5py.File(datafile_name, 'r')
		correlations = datafile['/correlations']

		dt = float(datafile.attrs['dt'])
		steps_per_invocation = float(datafile.attrs['stepsPerInvocation'])
		num_invocations = float(datafile.attrs['numInvocations'])
		sample_interval = dt*steps_per_invocation
		end_time = dt*steps_per_invocation*num_invocations
		times = np.arange(0,num_invocations+1)*sample_interval

		# process data
		meanvalues = np.mean(correlations, axis=1)
		deviations = np.std(correlations, axis=1)

		# plot data
		sparse=slice(None,None,5)
		ax1.plot(times[sparse],meanvalues[sparse],ls='none',ms=5,marker='o',color='#ffffff',label='Numerical simulation')
		ax2.semilogy(times[sparse],meanvalues[sparse],ls='none',ms=5,marker='o',color='#ffffff',label='Numerical simulation')
		lower_bound = np.clip(meanvalues-deviations,1e-4,np.inf)
		upper_bound = meanvalues+deviations
		ax2.fill_between(times,lower_bound,upper_bound,color=(.8,.8,.8,.5))

		# setup theory
		approx6 = approximation6(parameters['lambda'],parameters['betasqr'],times)
		approx2_renorm = approximation2_renorm(parameters['lambda'],parameters['betasqr'],times)

		#approx_exp = approximation_exp(parameters['lambda'],parameters['betasqr'],times)
		approx_flowdecorrelation = approximation_flowdecorrelation(parameters['lambda'],parameters['betasqr'],times)


		#plot theory
		ax1.plot(times,approx6,'k--',label='6th standard')
		ax1.plot(times,approx2_renorm,'r:', label='2nd renormalized')
		ax2.semilogy(times,approx6,'k--',label='6th standard')
		ax2.semilogy(times,approx2_renorm,'r:', label='2nd renormalized')

	ax1.set_xlabel(r'$t/\tau$')
	ax1.set_ylabel('C(t)',rotation=0)
	ax1.set_ylim(0,1)
	ax1.set_xlim(0,2.5)

	ax2.set_xlabel(r'$t/\tau$')
	ax2.set_ylabel('C(t)',rotation=0)
	ax2.set_ylim(1e-3,1)
	ax2.set_xlim(0,5)

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/corr2d_highbeta.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()