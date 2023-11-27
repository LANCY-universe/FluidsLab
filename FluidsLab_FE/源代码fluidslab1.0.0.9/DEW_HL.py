import math

#Calculates the pressure of water as a function of density and temperature
def calculatePressure(density,temperature,equation):
	m = 18.01528
	if equation == 1:
		ZD05_R = 83.144
		ZD05_Vc = 55.9480373
		ZD05_Tc = 647.25

		TK = temperature + 273.15
		Vr = m / density / ZD05_Vc
		Tr = TK / ZD05_Tc
		B = 0.349824207 - 2.91046273 / (Tr * Tr) + 2.00914688 / (Tr * Tr * Tr)
		C = 0.112819964 + 0.748997714 / (Tr * Tr) - 0.87320704 / (Tr * Tr * Tr)
		D = 0.0170609505 - 0.0146355822 / (Tr * Tr) + 0.0579768283 / (Tr * Tr * Tr)
		E = -0.000841246372 + 0.00495186474 / (Tr * Tr) - 0.00916248538 / (Tr * Tr * Tr)
		f = -0.100358152 / Tr
		g = -0.00182674744 * Tr

		delta = 1 + B / Vr + C / (Vr * Vr) + D / pow(Vr, 4) + E / pow(Vr, 5) + (f / (Vr * Vr) + g / pow(Vr, 4)) * math.exp(-0.0105999998 / (Vr * Vr))
		calculatePressure = ZD05_R * TK * density * delta / m
	elif equation == 2:
		ZD09_R = 0.083145
		ZD09_c1 = 6.971118009
		dm = 475.05656886 * density
		Vm = 0.0021050125 * (m / density)
		Tm = 0.3019607843 * (temperature + 273.15)
		B = 0.029517729893 - 6337.56452413 / (Tm * Tm) - 275265.428882 / (Tm * Tm * Tm)
		C = 0.00129128089283 - 145.797416153 / (Tm * Tm) + 76593.8947237 / (Tm * Tm * Tm)
		D = 2.58661493537E-06 + 0.52126532146 / (Tm * Tm) - 139.839523753 / (Tm * Tm * Tm)
		E = -2.36335007175E-08 + 0.00535026383543 / (Tm * Tm) - 0.27110649951 / (Tm * Tm * Tm)
		f = 25038.7836486 / (Tm * Tm * Tm)
		delta = 1 + B / Vm + C / (Vm * Vm) + D / pow(Vm, 4) + E / pow(Vm, 5) + f / (Vm * Vm) * (0.73226726041 + 0.015483335997 / (Vm * Vm)) * math.exp(-0.015483335997 / (Vm * Vm))
		Pm = ZD09_R * Tm * delta / Vm
		calculatePressure = Pm * ZD09_c1
	else:
		calculatePressure = 0
	return(calculatePressure)

#calculate water density
def calculateDensity(pressure,temperature,error,equation,Psat):
	if Psat == 1:
		calculateDensity = -1.01023381581205E-104 * pow(temperature, 40) + -1.1368599785953E-27 * pow(temperature, 10) + -2.11689207168779E-11 * pow(temperature, 4) + 1.26878850169523E-08 * pow(temperature, 3) + -4.92010672693621E-06 * pow(temperature, 2) + -3.2666598612692E-05 * temperature + 1.00046144613017
	else:
		minGuess = 0.00001
		guess = 0.00001
		maxGuess = 7.5 * equation - 5
		calcP = 0
		for i in range(500):
			calcP = calculatePressure(guess,temperature,equation)
			if abs(calcP - pressure) > error:
				if calcP > pressure:
					maxGuess = guess
					guess = (guess + minGuess) / 2
				elif calcP < pressure:
					minGuess = guess
					guess = (guess + maxGuess) / 2
			else:
				break
		calculateDensity = guess
	return(calculateDensity)

#calculate epsilon
def calculateEpsilon(density,temperature,equation,Psat):
	if Psat == 1:
		calculateEpsilon = -1.66686763214295E-77 * pow(temperature, 30) + \
		-9.02887020379887E-07 * pow(temperature, 3) + 8.4590281449009E-04 * pow(temperature, 2) + \
		-0.396542037778945 * temperature + 87.605024245432
	else:
		if equation == 1:
			T_hat = (temperature + 273.15) / 298.15
			k0 = 1
			k1 = 14.70333593 / T_hat
			k2 = 212.8462733 / T_hat - 115.4445173 + 19.55210915 * T_hat
			k3 = -83.3034798 / T_hat + 32.13240048 * T_hat - 6.69409865 * (T_hat * T_hat)
			k4 = -37.86202045 / (T_hat * T_hat) + 68.87359646 / T_hat - 27.29401652
			calculateEpsilon = k0 + k1 * density + k2 * (density * density) + k3 * pow(density, 3) + k4 * pow(density, 4)
		if equation == 2:
			pi = 3.14159265358979
			omega = 0.0000000268
			k = 1.380648E-16
			Na = 6.022E+23
			mu = 2.33E-18
			rhostar = (density * 0.055508) * pow(omega, 3) * Na
			mustarsq = pow(mu, 2) / (k * (temperature + 273.15) * pow(omega, 3))
			y = (4 * pi / 9) * rhostar * mustarsq
			f1 = 0.4341 * pow(rhostar, 2)
			f2 = -(0.05 + 0.75 * pow(rhostar, 3))
			f3 = -0.026 * pow(rhostar, 2) + 0.173 * pow(rhostar, 4)
			calculateEpsilon = ((3 * y) / (1 - f1 * y)) * (1 + (1 - f1) * y + f2 * (y * y) + f3 * (y * y * y)) + 1
		if equation == 3:
			N_k = [0.978224486826,-0.957771379375,0.237511794148,0.714692224396,-0.298217036956,-0.108863472196,\
			0.0949327488264,-0.00980469816509,0.000016516763497,9.37359795772E-5,-1.2317921872E-10,0.00196096504426]
			i_k = [1,1,1,2,3,3,4,5,6,7,10]
			j_k = [0.25,1,2.5,1.5,1.5,2.5,2,2,5,0.5,10]
			avogadro = 6.0221367E+23
			dipole = 6.138E-30
			epsilon_o = 8.8541878176204E-12
			boltzmann = 1.380658E-23
			alpha = 1.636E-40
			density_c = 17873.728
			T_c = 647.096
			density_molm3 = density * 0.055508 * 1000000
			T_K = temperature + 273.15
			g = 1
			for x in range(12):
				if x == 11:
					g = g + N_k[11] * (density_molm3 / density_c) * pow(T_K / 228 - 1, -1.2)
				else:
					g = g + N_k[x] * pow((density_molm3 / density_c), i_k[x]) * pow((T_c / T_K), j_k[x])

			A = (avogadro * pow(dipole, 2) * density_molm3 * g) / (epsilon_o * boltzmann * T_K)
			B = (avogadro * alpha * density_molm3) / (3 * epsilon_o)
			C = 9 + 2 * A + 18 * B + A * A + 10 * A * B + 9 * B * B
			calculateEpsilon = (1 + A + 5 * B + math.sqrt(C)) / (4 - 4 * B)
		if equation == 4:
			a1 = -1.57637700752506E-03
			a2 = 6.81028783422197E-02
			a3 = 0.754875480393944
			b1 = -8.01665106535394E-05
			b2 = -6.87161761831994E-02
			b3 = 4.74797272182151
			A = a1 * temperature + a2 * math.sqrt(temperature) + a3
			B = b1 * temperature + b2 * math.sqrt(temperature) + b3
			calculateEpsilon = math.exp(B) * pow(density, A)
	return(calculateEpsilon)

# temperature in these after two equations in the form Kelvin

# def A_gamma(pressure,temperature,equation,Psat):
# 	density = calculateDensity(pressure,temperature,0.01,equation,Psat)
# 	Epsilon = calculateEpsilon(density,temperature,1,Psat)
# 	return((1.824829238*10E6*(density)**0.5)/((Epsilon*(temperature+273.15))**(1.5))*0.1)

# def B_gamma(pressure,temperature,equation,Psat):
# 	density = calculateDensity(pressure,temperature,0.01,equation,Psat)
# 	Epsilon = calculateEpsilon(density,temperature,1,Psat)
# 	return((50.29158649*density**0.5)/(Epsilon*(temperature+273.15))**0.5)	


def one_Debye_Huckel(pressure,temperature,equation,Psat,ironic_strength,charge,ion_size,Ar,Br):
	I = ironic_strength
	z_i = charge
	aio = ion_size
	# Ar = A_gamma(pressure,temperature,equation,Psat)
	# Br = B_gamma(pressure,temperature,equation,Psat)
	logri = - (z_i**2*Ar*I**0.5)/(1+aio*Br*I**0.5)
	if z_i == 0:
		logri = 0
	# print(logri,'德拜虎克单点计算已运行')
	return round(logri,6)

def ironic_strength_function(mole_list,Z_list):
	# print(mole_list,Z_list)
	ironic_strength = 0
	for i in range(len(mole_list)):
		ironic_strength += 0.5*mole_list[i]*Z_list[i]**2

	# print(ironic_strength,'离子强度已运行')
	return ironic_strength

def all_Debye_huckel(P,T,mole_list,Z_list,Ar,Br):
	I = ironic_strength_function(mole_list,Z_list)
	logr_list_calculation = []
	for ii in range(len(mole_list)):
		logr_list_calculation.append(one_Debye_Huckel(P,T,1,0,I,Z_list[ii],4,Ar,Br))
	return logr_list_calculation

# Ar = A_gamma(5000,300,0.01,1,0)
# Br = B_gamma(5000,300,0.01,1,0)
# print(A,B)
# I = 0.5*(1+1)
# logr = - (Ar*I**0.5)/(1+4*Br*I**0.5)

# print(logr)