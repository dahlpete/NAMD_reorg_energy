import numpy as np
import matplotlib as mpl
mpl.use('tkagg')
from matplotlib import pyplot as plt


def free_energy_surface(data,temp,ev=None):
	if ev:
		data = data * 1.0/23.06035
		kB = 8.617e-5
	else:
		kB = 8.617e-5 * 23.06035
	pdf,bins = np.histogram(data,bins=500,density=True)
	bin_width = bins[-1] - bins[-2]
	bins = bins[0:-1] + bin_width

	free_energy = -kB * temp * np.log(pdf)

	return bins,free_energy

def reorganization_energy(data_ox,data_red,temp,ev=None):
	if ev:
		data_ox = data_ox * 1.0/23.06035 
		data_red = data_red * 1.0/23.06035 
		kB = 8.617e-5
	else:
		kB = 8.617e-5 * 23.06035
	r_st = 0.5 * (np.mean(data_red) - np.mean(data_ox))
	r_ox = np.std(data_ox)**2 / (2.0*kB*temp)
	r_red = np.std(data_red)**2 / (2.0*kB*temp)

	lambda_r = (r_st**2) / (0.5*(r_ox + r_red))

	r_st_std = np.sqrt((0.5 * np.std(data_ox))**2 + (0.5 * np.std(data_red)**2))

	error = (2.0 * r_st * r_st_std) / (0.5*(r_ox + r_red))
	

	#print(r'$\lambda^{St} =$ %.2d' % r_st
	#print(r'$\lambda^{OX} =$ %.2d' % r_ox)
	#print(r'$\lambda^{RED} =$ %.2d' % r_red)
	return r_st,r_ox,r_red,lambda_r,error

T = 310

redIN = open('data_red.txt','r').readlines()
oxIN = open('data_ox.txt','r').readlines()

shift = False

energy_red = np.zeros(len(redIN))
for i in range(len(redIN)):
	line = redIN[i]
	values = line.split()
	energy_red[i] = values[1]

energy_ox = np.zeros(len(oxIN))
for i in range(len(oxIN)):
	line = oxIN[i]
	values = line.split()
	energy_ox[i] = values[1]

[r_st,r_ox,r_red,lambda_r,error] = reorganization_energy(energy_ox,energy_red,T,ev=True)
print('r_st = %.3f eV' % r_st)
print('r_ox = %.3f eV' % r_ox)
print('r_red = %.3f eV\n' % r_red)
print('_________________________________\n')
print('lambda_r = %.3f eV\n' % lambda_r)
print('lambda_r_err = %.3f eV\n' % error)
print('_________________________________\n')

[bin_ox,oxidized_surface] = free_energy_surface(energy_ox,310,ev=True)
[bin_red,reduced_surface] = free_energy_surface(energy_red,310,ev=True)
mpr = 0.5*(np.max(bin_red)+np.min(bin_red)); mpo = 0.5*(np.max(bin_ox)+np.min(bin_ox))

idx = np.isfinite(bin_ox) & np.isfinite(oxidized_surface)
fit_ox = np.polyfit(bin_ox[idx],oxidized_surface[idx],2)
parabola_ox = lambda x: fit_ox[0] * x**2 + fit_ox[1] * x + fit_ox[2]
x_ox = np.linspace(1.75*mpo-0.75*mpr,0.25*mpo+0.75*mpr,100)

idx = np.isfinite(bin_red) & np.isfinite(reduced_surface)
fit_red = np.polyfit(bin_red[idx],reduced_surface[idx],2)
parabola_red = lambda x: fit_red[0] * x**2 + fit_red[1] * x + fit_red[2]
x_red = np.linspace(0.25*mpr+0.75*mpo,1.75*mpr-0.75*mpo,100)

if shift:
	solutions = np.roots(fit_red - fit_ox)
	mean_values = [np.mean((energy_red * 1.0/23.06035)),np.mean((energy_ox)*1.0/23.06035)]
	idx1 = solutions < np.max(mean_values)
	idx2 = solutions > np.min(mean_values)
	idx = idx1 * idx2
	zero_point = solutions[idx]; print(zero_point)
	bin_ox = bin_ox - zero_point; bin_red = bin_red - zero_point

	idx = np.isfinite(bin_ox) & np.isfinite(oxidized_surface)
	fit_ox = np.polyfit(bin_ox[idx],oxidized_surface[idx],2)
	parabola_ox = lambda x: fit_ox[0] * x**2 + fit_ox[1] * x + fit_ox[2]

	idx = np.isfinite(bin_red) & np.isfinite(reduced_surface)
	fit_red = np.polyfit(bin_red[idx],reduced_surface[idx],2)
	parabola_red = lambda x: fit_red[0] * x**2 + fit_red[1] * x + fit_red[2]

	x_ox = x_ox - zero_point; x_red = x_red - zero_point


fig = plt.figure(figsize=(5,5))
ax = fig.gca()

ax.scatter(bin_ox,oxidized_surface,marker='o',c='g',linewidths=0,s=5,label='OX')
ax.scatter(bin_red,reduced_surface,marker='o',c='r',linewidths=0,s=5,label='RED')
ax.plot(x_ox,parabola_ox(x_ox),color='g',linewidth=2)
ax.plot(x_red,parabola_red(x_red),color='r',linewidth=2)
plt.xlabel('energy gap (eV)',fontsize=14)
plt.ylabel('free energy (eV)',fontsize=14)
plt.legend(loc=1,fontsize=16,markerscale=3,framealpha=1)

#plt.show()
