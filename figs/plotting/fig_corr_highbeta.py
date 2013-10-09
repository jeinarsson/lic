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

	part_two = 1 + betasqr*(-4)*(1+l**2)**(-2)*(2+l**2+2*l**4)*((-1)+ np.exp(-ts)+ts)
	part_four = betasqr**2*2*(1+l**2)**(-4)*(2+l**2+2*l**4)*((-87)+  \
			np.exp(-2*ts)*(27+(-39)*l**2+27*l**4)+22*ts+8*ts**2+ \
			l**2*(219+(-94)*ts+4*ts**2)+l**4*((-87)+22*ts+8* \
			ts**2)+4*np.exp(-ts)*(15+23*ts+l**4*(15+23*ts)+(-1)* \
			l**2*(45+41*ts))) 
	part_six = betasqr**3*(float(1)/9)*(1+l**2)**(-6)*(2+l**2+2*l**4)*(np.exp(-3*ts)*( \
			(-9268)+22862*l**2+(-28188)*l**4+22862*l**6+( \
			-9268)*l**8)+4*(29713+(-81257)*l**2+107838* \
			l**4+(-81257)*l**6+29713*l**8)+(-3)*(6037+( \
			-20558)*l**2+36342*l**4+(-20558)*l**6+6037* \
			l**8)*ts+(-2160)*(2+(-5)*l**2+l**4+(-5)* \
			l**6+2*l**8)*ts**2+(-96)*(2+l**2+2*l**4) \
			**2*ts**3+(-9)*np.exp(-2*ts)*(5228+5833*ts+6*l**4*(4828+3533* \
			ts)+l**8*(5228+5833*ts)+(-2)*l**2*(10106+8191*ts)+( \
			-2)*l**6*(10106+8191*ts))+(-18)*np.exp(-ts)*(3474+(-6681)* \
			l**2+527*ts*(12+5*ts)+(-2)*l**2*ts*(9543+2800*ts) \
			+l**8*(3474+527*ts*(12+5*ts))+6*l**4*(1319+67* \
			ts*(62+15*ts))+(-1)*l**6*(6681+2*ts*(9543+2800*ts))));

	return part_two + part_four + part_six

def approximation2_renorm(l, betasqr, ts):

	b = (l**2-1.0)/(l**2+1.0)
	k2 = -(3.*b**2+5.0)*(np.exp(-ts) + (-1.+ts))

	part_two = np.exp(betasqr*k2)

	return part_two

def approximation4_renorm(l, betasqr, ts):

	b = (l**2-1.0)/(l**2+1.0)
	k2 = -(3.*b**2+5.0)*(np.exp(-ts) + (-1.+ts))
	k4 = 1./8. * (3.*b**2+5.0) *(-5.+81*b**2) * (np.exp(-2.*ts) + 4.*np.exp(-ts)*(1.+ts)+(-5.+2.*ts))
	
	return np.exp(betasqr*k2 + betasqr**2*k4)	

def approximation6_renorm(l, betasqr, ts):

	b = (l**2-1.0)/(l**2+1.0)
	k2 = -(3.*b**2+5.0)*(np.exp(-ts) + (-1.+ts))
	k4 = 1./8. * (3.*b**2+5.0) *(-5.+81*b**2) * (np.exp(-2.*ts) + 4.*np.exp(-ts)*(1.+ts)+(-5.+2.*ts))
	k6 = -1./144. * (3.*b**2+5.0)* ( \
		25.*(4.*np.exp(-3.*ts) + 9.*np.exp(-2.*ts)*(4.+3*ts) + 18.*np.exp(-ts)*(2.+4*ts+ts**2) +(-76.+21*ts) ) + \
		27.*b**4*(686.*np.exp(-3.*ts) + 9.*np.exp(-2.*ts)*(688.+515.*ts) + 18.*np.exp(-ts)*(341.+686*ts+171.*ts**2) + (-13016+3597*ts)) + \
		30.*b**2*(199.*np.exp(-3.*ts) + np.exp(-2.*ts)*(-270+828*ts) + 9.*np.exp(-ts)*(657+398*ts+214*ts**2) + 2.*(-2921+780*ts) ))
	return np.exp(betasqr*k2 + betasqr**2*k4 + betasqr**6*k6)		

			
def approximation_exp(l, betasqr, ts):
	return np.exp(-2*betasqr*(4*l**4+2*l**2+4)/(l**2+1)**2 * (ts));

def approximation_mw_slow(l, betasqr, ts):
	return np.exp(-(np.sqrt(17)-1)*0.5* (ts));

def approximation_ein_slow(l, betasqr, ts):
	return np.exp(-np.pi/2. * ts);

def approximation_flowdecorrelation(l, betasqr, ts):
	return 0.5*np.exp((-ts));

def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')

	fig1 = plt.figure(1, figsize=(13/2.5, 6/2.5))
	ax1 = fig1.add_subplot(1,2,1)
	ax2 = fig1.add_subplot(1,2,2)


	all_parameters = []
	all_parameters.append({
		'lambda': 1,
		'betasqr': .9
	})
	

	for parameters in all_parameters:
		datafile_name = "paper_all/" + "correlations_l%(lambda)e_betasqr%(betasqr)e.hdf5" % parameters
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
		sparse=slice(None,None,15)
		ax1.plot(times[sparse],meanvalues[sparse],ls='none',ms=5,marker='o',color='#ffffff',label='Numerical simulation')
		ax2.semilogy(times[sparse],meanvalues[sparse],ls='none',ms=5,marker='o',color='#ffffff',label='Numerical simulation')
		lower_bound = np.clip(meanvalues-deviations,1e-4,np.inf)
		upper_bound = meanvalues+deviations
		ax2.fill_between(times,lower_bound,upper_bound,color=(.8,.8,.8,.5))

		# setup theory
		approx6 = approximation6(parameters['lambda'],parameters['betasqr'],times)
		approx2_renorm = approximation2_renorm(parameters['lambda'],parameters['betasqr'],times)
		approx4_renorm = approximation4_renorm(parameters['lambda'],parameters['betasqr'],times)
		approx6_renorm = approximation6_renorm(parameters['lambda'],parameters['betasqr'],times)
		#approx_exp = approximation_exp(parameters['lambda'],parameters['betasqr'],times)
		approx_mw_slow = approximation_mw_slow(parameters['lambda'],parameters['betasqr'],times)
		approx_ein_slow = approximation_ein_slow(parameters['lambda'],parameters['betasqr'],times)
		approx_flowdecorrelation = approximation_flowdecorrelation(parameters['lambda'],parameters['betasqr'],times)

		#plot theory
		ax1.plot(times,approx6,'k--',label='6th standard')
		ax1.plot(times,approx2_renorm,'r:', label='2nd renormalized')
		ax1.plot(times,approx4_renorm,'g:', label='4th renormalized')
		ax1.plot(times,approx6_renorm,'m:', label='6th renormalized')
		#ax2.semilogy(times,approx_mw_slow,'k--')
		#ax2.semilogy(times,approx_flowdecorrelation,'k:')
		ax2.semilogy(times,approx6,'k--',label='6th standard')
		ax2.semilogy(times,approx2_renorm,'r:', label='2nd renormalized')
		ax2.semilogy(times,approx4_renorm,'g:', label='4th renormalized')
		ax2.semilogy(times,approx6_renorm,'m:', label='6th renormalized')
		#ax2.semilogy(times,.5*np.exp(-np.pi**2/6. * times),'c:', label='$-\pi^2/6$')
		ax2.semilogy(times,.5*approx_mw_slow,'y:', label='"MW Slow mode"')
		

	ax1.set_xlabel(r'$t/\tau$')
	ax1.set_ylabel(r'$\langle n(0)n(t)\rangle$',rotation=90)
	ax1.set_ylim(0,1)
	ax1.set_xlim(0,2)

	ax2.set_xlabel(r'$t/\tau$')
	ax2.set_ylabel(r'$\langle n(0)n(t)\rangle$',rotation=90)
	ax2.set_ylim(1e-3,1)
	ax2.set_xlim(0,4)
	ax2.set_xticks([0,1,2,3,4])

	#ax1.legend()

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/corr_highbeta.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()