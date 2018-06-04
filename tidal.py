from astropy import constants as const
from pylab import *
import scipy
from scipy import interpolate
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc,rcParams

rc('font', **{'family': 'Helvetica'})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size':20})
plt.rc('legend', **{'fontsize':7})

# Ticks to the outside:
rcParams['axes.linewidth'] = 3.0
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

def get_ecc_change(a,e):
	A = (63./4.)*np.sqrt(G*Ms**3)*(Rp**5/(Qp*Mp))
	B = (225./16.)*np.sqrt(G/Ms)*(Rs**5*Mp/Qs)
	dedt = -e*(A+B)*a**(-6.5)
	#dedt = -e*(63*Rp**5.*np.sqrt(G*Ms**3)/(4.*Qp*Mp) + 171.*Rs**5*Mp*np.sqrt(G/Ms)/(16.*Qs)) * a**(-13./2.)
	return dedt

def get_sma_change(a,e):
	A = e**2*(63./2.)*np.sqrt(G*Ms**3)*(Rp**5/(Qp*Mp))
	B = (9.*np.sqrt(G/Ms)*Rs**5*Mp/(2.*Qs)) * (1. + (57./4.)*e**2)
	dadt = -a*(A+B) * a**(-6.5)
	return dadt

def kepler3law(sma,Ms,Mp):
	period = np.sqrt(4.*np.pi**2*sma**3/(G*(Ms+Mp)))
	return period



forward = True
backward = True

G    = const.G.value
Msun = const.M_sun.value
Rsun = const.R_sun.value
Mjup = const.M_jup.value
Rjup = const.R_jup.value
AU   = const.au.value

GyrToSec = 10**9 * 365. * 24 * 3600.
SecToDay = 3600*24.

Mp0,MpE =     0.179*Mjup, 0.02*Mjup
Rp0,RpE =     0.840*Rjup, 0.011*Rjup
Ms0,MsE =     1.101*Msun, 0.015*Msun
Rs0,RsE =     1.669*Rsun, 0.021*Rsun
Qp = 3.2e4
Qs = 1e6

sma0,smaE = 0.10376 * AU, 0.0005 * AU
ecc0,eccE = 0.42, 0.03
age0,ageE = 8.51, 0.6

realizations = 500
itera = 0

d = np.loadtxt('CL002_14_Rs.txt')
tck = interpolate.splrep(d[:,0],d[:,1],k=1)
f, (ax2, ax3) = plt.subplots(2, 1, sharex=True,figsize=(10,12))
all_roches = []
rochetimes = []

while itera < realizations:

	if itera != realizations - 1:
		while True:
			Mp = np.random.normal(Mp0,MpE)
			if Mp>0:
				break
		Rp = np.random.normal(Rp0,RpE)
		Ms = np.random.normal(Ms0,MsE)
		Rs = np.random.normal(Rs0,RsE)
		sma = np.random.normal(sma0,smaE)
		while True:
			ecc = np.random.normal(ecc0,eccE)
			if ecc>=0 and ecc <= 1:
				break
		while True:
			age = np.random.normal(age0,ageE)
			if age < 9.5:
				break

	else:
		Mp = Mp0
		Rp = Rp0
		Ms = Ms0
		Rs = Rs0
		sma = sma0
		ecc = ecc0
		age = age0	

	A = (63./4.)*np.sqrt(G*Ms**3)*(Rp**5/(Qp*Mp))
	tau_circ = A * sma0**(-6.5)
	tau_circ = 1. / tau_circ
	aroche = 2.7 * Rp * (Ms/Mp)**(1./3.)

	print '\n'
	print 'iteracion:', itera
	print 'Mp =', Mp
	print 'Rp =', Rp
	print 'Ms =', Ms
	print 'Rs =', Rs
	print 'Age =', age


	print 'a =', sma
	print 'ecc =', ecc
	print 'tau_cir', tau_circ / GyrToSec 

	print Qp
	print Qs

	print '\n'

	per = kepler3law(sma,Ms,Mp)
	Qpp = 3*Qp/(2.*0.37)
	H0 = (63./4.)*sma**(-15./2.)*ecc**2*((G*Ms)**(3./2.)*Ms*Rp**5)/Qpp
	#print H0
	#print gfcd
	times,eccs,smas,pers = [],[],[],[]

	t = age * GyrToSec#10Gyr
	dt = 0.0001 * GyrToSec
	tf = 9.80 * GyrToSec
	oecc,osma = ecc, sma
	while t<tf:
		if t/GyrToSec < d[0,0]:
			Rs = d[0,1]*Rsun
		elif t/GyrToSec > d[-1,0]:
			Rs = d[-1,1]*Rsun
		else:
			Rs = interpolate.splev(t/GyrToSec,tck)*Rsun
		#Rs = Rs0
		dadt = get_sma_change(sma,ecc)
		dedt = get_ecc_change(sma,ecc)
		#print dadt,dt,dadt*dt
		#print dedt,dt,dedt*dt
		ecc = ecc + dedt * dt
		sma = sma + dadt * dt
		#print sma*(1-ecc)/AU, aroche/AU
		#print sma*(1-ecc) < aroche

		#print ecc, sma*(1-ecc) 
		if ecc < 0 or sma*(1-ecc) < Rs:
			break

		per = kepler3law(sma,Ms,Mp)

		#print ecc,sma
		#print gfd
		times.append(t)
		eccs.append(ecc)
		smas.append(sma)
		pers.append(per)
		#print t/GyrToSec,Rs/Rsun,ecc,sma,per/SecToDay

		#if len(pers)==2:
		#	print pers[0]/SecToDay
		#	print ((pers[1] - pers[0]) / dt) * pers[0]*100
		#	print gfd
		t += dt

	t = age * GyrToSec#10Gyr
	dt = -0.0001 * GyrToSec
	tf = .100 * GyrToSec
	ecc,sma = oecc,osma
	rochetime = -999
	rdone = False
	while t>tf:
		#print 'hola', t
		if t/GyrToSec < d[0,0]:
			Rs = d[0,1]*Rsun
		elif t/GyrToSec > d[-1,0]:
			Rs = d[-1,1]*Rsun
		else:
			Rs = interpolate.splev(t/GyrToSec,tck)*Rsun
		#Rs = Rs0

		dadt = get_sma_change(sma,ecc)
		dedt = get_ecc_change(sma,ecc)
		#print dadt,dt,dadt*dt
		#print dedt,dt,dedt*dt
		#print t,ecc,sma,dedt*dt,dadt*dt
		ecc = ecc + dedt * dt
		sma = sma + dadt * dt

		if ecc > 1 or sma < 0:
			break

		if sma*(1-ecc) < aroche and rdone == False:
			rochetime = t
			print 'Roche Time:', rochetime/GyrToSec
			rochetimes.append(rochetime)
			rdone = True
		per = kepler3law(sma,Ms,Mp)
		#print ecc,sma
		#print gfd
		times.append(t)
		eccs.append(ecc)
		smas.append(sma)
		pers.append(per)
		#print t/GyrToSec,Rs/Rsun,ecc,sma,per/SecToDay
		t += dt

	times,eccs,smas,pers = np.array(times),np.array(eccs),np.array(smas),np.array(pers)
	I = np.argsort(times)
	times = times[I]
	eccs = eccs[I]
	smas = smas[I]
	pers = pers[I]

	if itera != realizations - 1:
		ax2.plot(times/GyrToSec,eccs, 'k', linewidth=3.,alpha=0.1)
		ax2.set_ylabel('ecc')
		#ax2.set_xlabel('time [Gyr]')
		ax2.plot(age,oecc,'ro',alpha=0.3)
		#show()
		ax3.plot(times/GyrToSec,smas/AU, 'k', linewidth=3.,alpha=0.1)
		#ax2.axhline(aroche/AU,c='r',alpha=0.3)
		ax3.plot(times/GyrToSec,smas*(1-eccs)/AU,'b', linewidth=3.,alpha=0.1)
		ax3.plot(age,osma/AU,'ro',alpha=0.3)
	else:
		ax2.plot(times/GyrToSec,eccs, 'y', linewidth=2.,alpha=0.8)
		ax2.set_ylabel('ecc')
		ax2.plot(age,oecc,'ro',alpha=0.8)
		#show()
		ax3.plot(times/GyrToSec,smas*(1-eccs)/AU,'y--', linewidth=2.,alpha=0.7)
		ax3.plot(times/GyrToSec,smas/AU, 'y', linewidth=2.,alpha=0.8)
		#ax2.axhline(aroche/AU,c='r',alpha=0.8)
		ax3.plot(age,osma/AU,'ro',alpha=0.8)
	all_roches.append(aroche/AU)

	itera += 1

for roche in all_roches:
	ax3.axhline(roche,c='r',alpha=0.3)
ax3.axhline(all_roches[-1],c='r',alpha=0.8,linewidth=2.)

ax3.plot(d[:,0],d[:,1]*Rsun/AU,'g', linewidth=3.,alpha=0.7)
ax3.plot([0,1],[d[0,1]*Rsun/AU,d[0,1]*Rsun/AU],'g', linewidth=3.,alpha=0.7)
ax3.set_ylabel('a [AU]')
ax3.set_xlabel('time [Gyr]')

ax3.set_ylim([0.,0.25])
ax2.set_ylim([0.,1.])
ax2.set_xlim([0.1,10.])
ax3.set_xlim([0.1,10.])

#ax4.set_xlabel('time [Gyr]')
#ax4.set_ylabel('N')

#rochetimes = np.array(rochetimes)
#print rochetimes
#I = np.where(rochetimes!=-999)[0]
#ax4.invert_yaxis()
#ax4.hist(rochetimes[I]/GyrToSec,color='black',alpha=0.3)
#print float(len(I))/float(len(rochetimes))

#print np.mean(rochetimes[I])/GyrToSec

f.subplots_adjust(hspace=.0)

#show()
#ax1.plot(d[:,0],d[:,1])
#ax1.set_ylabel(r'R$_{\odot}}$')

#ax1.plot(times/GyrToSec,pers/SecToDay)
#ax1.set_ylabel(r'P [day]')
#ax1.plot(age,per0/SecToDay,'ro')
#suptitle('log(Qp)='+str(np.log10(Qp)), fontsize=16)

#H = (63./4.)*smas**(-15./2.)*eccs**2*((G*Ms)**(3./2.)*Ms*Rp**5)/Qpp
#ax1.plot(times/GyrToSec,np.log10(H))
#ax1.set_ylabel(r'H [W]')
#ax1.plot(age,np.log10(H0),'ro')
#suptitle('log(Qp)='+str(np.log10(Qp)), fontsize=16)


#print rochetime/GyrToSec
plt.tight_layout()
plt.savefig('tidal_dist.pdf', dpi=300)
