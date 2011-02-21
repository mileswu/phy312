#!/usr/bin/python
import sys
import math
import numpy

def read_file(fn):
	f = open(fn, 'r')

	l1 = f.readline().strip() #v_min v_max v_incr
	l2 = f.readline().strip() #time
	time = float(l2)

	i = 0
	for l in f:
		i += 1
	
	f.seek(0)
	f.readline()
	f.readline()
	
	current = numpy.zeros(i)
	rate = numpy.zeros(i)

	i = 0
	for l in f:
		l_split = l.strip().split() #current count

		current[i] = float(l_split[0]) 
		rate[i] = float(l_split[1])/time
		
		i += 1

	err_rate = numpy.sqrt(rate)
	return current, rate, err_rate

def bk_analyse(bk):
	current, rate, err_rate = bk
	err = numpy.sqrt(numpy.sum(err_rate*err_rate)) / len(current) #error prop
	mean = numpy.mean(rate)
#	print "Background count - %f +/- %f" % (mean, err)
	return mean, err

def na_analyse(na, bk_mean, bk_err):
	current, rate, err_rate = na
	scaling = 0.0897
	me = 0.51091
	z = 11 #na
	sign = -1.0 #positrons

	momentum = current*scaling #in MeV/c
	rate = rate - bk_mean
	for i in range(len(rate)):
		if(rate[i] < 0):
			rate[i] = 0
	err_rate = numpy.sqrt(err_rate + bk_err) #error prop

	energy = numpy.sqrt(momentum*momentum + me*me) - me #KE in MeV

	eta = momentum/me
	alpha = 1.0/137.0
	SS = math.pow(1-(alpha*alpha*z*z), 0.5)
	delta = alpha*z*numpy.power(eta*eta+1, 0.5) / eta
	fermi = numpy.zeros(len(current))
	for i in range(len(current)):
		fermi[i] = math.pow(eta[i], 2*SS)*math.exp(sign * math.pi * delta[i])*abs_gamma(SS, delta[i])
	
	kurie = numpy.sqrt(rate / fermi / momentum / momentum)
	err_kurie = 0.5 * err_rate / numpy.sqrt(rate * fermi * momentum * momentum)
	for i in range(len(kurie)):
		if(kurie[i] == 0):
			err_kurie[i] = 0
	
	return energy, momentum, rate, err_rate, kurie, err_kurie


def gamma(x):
	g = 1.0;
	gamma = 1e99;
	
	while (x != 0.0):
		g = g*x;
		x += 1;
		if (x >= 3):
			gamma = (1.0 - 2.0/7.0/x/x*(1.0 - 2.0/3.0/x/x))/30.0/x/x;
			gamma = (1.0 - gamma)/12.0/x + x*(math.log(x) - 1.0);
			gamma = math.exp(gamma)/g*math.sqrt((2*math.pi)/x);
			return gamma;

def abs_gamma(x, y):
	f = 1.0;
	for j in range(1,700):
		f *= abs(x + j) / math.sqrt(y*y + (j+x)*(j+x));
	AbsGamma = f*f*gamma(x)*gamma(x);
	return AbsGamma;




if(len(sys.argv) != 3):
	print("Usage: %s <na_filename> <bk_filename>" % sys.argv[0])
	sys.exit()

bk_filename = sys.argv[2]
na_filename = sys.argv[1]

bk = read_file(bk_filename)
bk_mean, bk_err = bk_analyse(bk)

na = read_file(na_filename)
na = na_analyse(na, bk_mean, bk_err)
energy, momentum, rate, err_rate, kurie, err_kurie = na

for i in range(len(momentum)):
	print "%f %f %f %f %f %f" % (energy[i], momentum[i], rate[i], err_rate[i], kurie[i], err_kurie[i])

