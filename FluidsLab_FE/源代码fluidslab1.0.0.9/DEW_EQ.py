import math


def RationalTaylor(z0,A01,B01,B02,C02,A1,A2,B1,B2,C2,P,T):
	# P in kbar T in oC
	return (z0+A01*P+B01*T+B02*T**2+C02*P*T)/(1+A1*P+B1*T+A2*P**2+B2*T**2+C2*T*P)

def Rational2D(z0,A01,B01,B02,B03,A1,A2,A3,B1,B2,P,T):
	# P in kbar T in oC
	return (z0+A01*P+B01*T+B02*T**2+B03*T**3)/(1+A1*P+A2*P**2+A3*P**3+B1*T+B2*T**2)	

def ExtremeCum(z0,B,C,D,E,F,G,H,P,T):
	# P in kbar T in oC
	return z0+B*math.exp(-math.exp(-(P-C)/D))+E*math.exp(-math.exp(-(T-F)/G))+H*math.exp(-math.exp(-(P-C)/D) -math.exp(-(T-F)/G))

def gamma_a(P,T):
	# P in kbar T in oC
	return RationalTaylor(0.52716,0.0303,1.95906E-5,-3.02075E-7,3.0844E-6,0.11435,3.67273E-4,-0.00193,9.13888E-7,-7.31377E-5,P,T)

def gamma_b(P,T):
	# P in kbar T in oC
	return RationalTaylor(0.34155,0.02159,-3.89476E-4,9.4728E-8,-6.32174E-6,0.07003,2.95282E-5,-0.00143,5.20864E-7,-3.14848E-5,P,T)

def logk(P,T,species_name):
	# P in kbar T in oC	
	if species_name == '<OH->':
		return Rational2D(11.63286,0.15341,-0.02218,1.888E-5,-4.83456E-9,0.03213,4.55078E-6,1.55319E-6,-0.00148,7.61407E-7,P,T)
	if species_name == '<H2CO3>':
		return RationalTaylor(920.79375,353.02087,-67.95418,-0.28862,2.0126,115.46896,-2.98907,21.0776,0.00463,0.05991,P,T)
	if species_name == '<HCO3->':
		return Rational2D(-10.07275,0.17443,-0.00393,4.40638E-6,2.45698E-10,0.01615,-7.58418E-4,4.54301E-6,1.87832E-4,-5.66887E-7,P,T)
	if species_name == '<CaHCO3+>':
		return Rational2D(2.67072,0.07948,-0.03817,-2.22127E-5,2.31785E-8,0.18387,-0.00441,3.73307E-5,0.00631,-4.39946E-6,P,T)
	if species_name == '<CaCO3>':
		return Rational2D(-1.76354,0.15088,-0.02574,8.54172E-6,2.45234E-9,0.07085,-0.00184,1.40307E-5,0.00136,-1.35571E-6,P,T)
	if species_name == '<CO2>':
		return Rational2D(474.10783,-39.21757,2.85953,0.00632,-4.40326E-6,-7.03204,0.17529,-0.00142,-0.81054,4.59957E-4,P,T)
	if species_name == '<CaOH+>':
		return Rational2D(10.32759,0.02395,-0.02178,1.84795E-5,-6.6613E-9,0.00363,1.2977E-4,-1.14818E-6,-5.81161E-4,-2.92904E-7,P,T)
	if species_name == '<MgCO3>':
		return RationalTaylor(6.25792E12,3.97555E11,2.98238E10,2.52666E7,-3.73033E9,-2.98568E11,1.78316E9,-1.2325E10,8966545.88129,-4.62011E8,P,T)
	if species_name == '<MgHCO3+>':
		return Rational2D(407.69101,30.80209,-4.18957,-0.01221,7.09914E-6,79.95876,-1.77094,0.01397,2.25925,-0.00134,P,T)
	if species_name == '<MgOH+>':
		return Rational2D(11.1216,-0.00161,-0.01741,7.0778E-6,-5.68459E-10,-0.00522,1.34312E-4,-1.36343E-6,9.92649E-4,-1.82147E-6,P,T)
	if species_name == '<NaCO3->':
		return RationalTaylor(-23.47745,3.88005,0.39371,5.94669E-4,-0.03656,-3.68956,0.09496,-0.20431,1.73291E-4,-0.0198,P,T)
	if species_name == '<NaHCO3>':
		return RationalTaylor(-31.88857,-24.19592,0.26562,4.3703E-4,0.06347,-12.98957,0.01147,-0.19639,8.14603E-5,-0.00528,P,T)
	if species_name == '<NaOH>':
		return Rational2D(12.57105,-0.013,-0.0112,1.06059E-5,-3.78181E-9,0.00778,-7.88228E-5,4.57349E-7,2.48011E-4,-4.27896E-8,P,T)
	if species_name == '<CH4>':
		return RationalTaylor(-1471.60278,-2308.84166,-1.30023,0.00587,1.87201,-17.35997,0.00461,-0.11386,5.48557E-5,-0.05089,P,T)
	if species_name == '<ETHANE>':
		return RationalTaylor(-1642.44488,-3791.04089,2.26441,0.00437,3.14025,-15.45452,-0.03239,-0.07492,1.11339E-4,-0.04882,P,T)
	if species_name == '<CaCl+>':
		return Rational2D(63.65869,2.73882,-0.25727,-0.00121,7.49865E-7,8.96864,-0.23235,0.00209,0.17497,-1.2461E-4,P,T)
	if species_name == '<CaCl2>':
		return RationalTaylor(1.07984E7,4002612.64064,-50312.40183,-139.42982,-3905.95526,2109334.47387,-15539.07682,12895.10036,-5.05973,127.50845,P,T)
	if species_name == '<MgCl+>':
		return Rational2D(0.89433,0.02022,-0.00651,-1.2138E-5,6.08398E-9,0.22233,-0.00529,4.6987E-5,4.81664E-4,-5.48488E-7,P,T)
	if species_name == '<NaCl>':
		return Rational2D(0.87305,-0.29521,-0.00326,1.47561E-5,-5.16637E-9,-0.42898,0.00813,-7.55121E-5,-0.00276,1.64475E-6,P,T)
	if species_name == '<H8Si3O10>':
		return RationalTaylor(-20.67804,-11.31019,-0.35477,3.82429E-4,0.00658,-1.84879,0.00258,-0.08217,-8.28285E-5,-0.00221,P,T)
	if species_name == '<H6Si2O7>':
		return Rational2D(-1.46578,0.00459,0.00475,-5.27922E-6,1.9717E-9,0.00908,-1.70147E-4,1.68764E-6,-0.002,9.69127E-7,P,T)
	if species_name == '<H3SiO4->':
		return RationalTaylor(8.08455,0.09139,-0.0154,8.28532E-6,5.0784E-4,0.01979,7.02713E-4,-0.00194,9.44838E-7,5.43106E-5,P,T)
	if species_name == '<CaH3SiO4+>':
		return RationalTaylor(-127414.35503,-21478.35365,159.10317,0.1315,16.178,-2922.55136,27.0408,-80.7373,0.04302,-3.61827,P,T)
	if species_name == '<MgH3SiO4+>':
		return Rational2D(6.43314,-0.27621,0.00825,1.08933E-5,1.33074E-9,-0.02446,7.72832E-4,-5.30822E-6,0.00174,1.68165E-6,P,T)
	if species_name == '<HS->':
		return Rational2D(139.71308,0.05429,-0.22372,7.92523E-5,7.39762E-9,-0.00658,4.04162E-4,-5.46465E-6,0.00241,-2.94446E-6,P,T)
	if species_name == '<CaSO4>':
		return Rational2D(-1.90236,0.15346,-0.00359,-2.68086E-7,-2.35681E-11, 0.12991,-0.00243,1.95794E-5,-0.00131,3.95082E-7,P,T)
	if species_name == '<MgSO4>':
		return RationalTaylor(-134.37759,-11.02853,1.75022,6.2889E-5,-0.00752,-6.372,0.16797,-0.21007,1.32533E-4,-0.00393,P,T)
	if species_name == '<NaSO4->':
		return Rational2D(170.48522,0.42455,-0.23707,1.19492E-4,-4.18929E-8,0.00261,3.18514E-4,-3.69641E-6,0.0012,-1.88651E-6,P,T)
	if species_name == '<HSO4->':
		return Rational2D(-0.84284,0.0719,-0.0178,1.25009E-5,-1.08354E-9,0.03315,-0.00145,1.8041E-5,-9.69631E-5,-4.30549E-7,P,T)


def logQ_mineral(P,T,mineral_name):
	# P in kbar T in oC	
	if mineral_name == '[Calcite]':
		return Rational2D(-7.91017,0.26146,-0.02408,9.97885E-6,1.48765E-9,0.05315,-0.00142,9.18933E-6,7.40507E-4,-9.17496E-7,P,T)
	elif mineral_name == '[Aragonite]':
		return Rational2D(-7.77244,0.1975,-0.01722,8.65543E-6,1.77664E-10,0.04415,-0.0012,8.16899E-6,1.76256E-4,-5.30966E-7,P,T)
	elif mineral_name == '[Dolomite]':
		return Rational2D(-16.70306,0.35858,-0.03065,2.12664E-5,-2.53257E-9,0.03704,-0.00106,7.50027E-6,-7.55234E-5,-3.67991E-7,P,T)
	elif mineral_name == '[Magnesite]':
		return Rational2D(-8.18952,0.13426,-0.01203,1.17614E-5,-3.10989E-9,0.03081,-8.95701E-4,6.73259E-6,-4.72551E-4,-8.74112E-8,P,T)
	elif mineral_name == '[Graphite]':
		return RationalTaylor(-348.93099,-674.24621,-0.30631,0.00322,0.79254,-11.344,-0.05225,-0.10154,6.44717E-5,-0.03697,P,T)
	elif mineral_name == '[Quartz]':
		return Rational2D(-3.20946,0.00202,0.00997,-1.0537E-5,3.7638E-9,0.0136,-3.68528E-4,3.52205E-6,-0.00185,7.97698E-7,P,T)
	elif mineral_name == '[Coesite]':
		return ExtremeCum(-1.81511,-24.27469,0.09708,5.4945,0.82074,-891.97116,537.11567,26.47817,P,T)
	elif mineral_name == '[Wollastonite]':
		return RationalTaylor(-17037.92036,-2302.57336,-8.35923,0.03115,-1.1679,-193.87293,1.61551,-11.1354,0.00975,-0.40916,P,T)
# print(gamma_b(50,300))
# print(logk(50,500,'<CaH3SiO4+>'))
# print(logQ_mineral(50,500,'[Wollastonite]'))
