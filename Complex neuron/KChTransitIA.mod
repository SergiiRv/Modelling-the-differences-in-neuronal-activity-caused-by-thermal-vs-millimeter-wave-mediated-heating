TITLE transient potassium current (A-current)

COMMENT
AP back-prop. explains threshold variability and rapid rise (McCormick et al. 2007, Yu et al. 2008)
	*********************************************
	reference:	Huguenard & McCormick (1992) 
			J.Neurophysiology 68(4), 1373-1383
	found in:	thalamic relay neurons		 	
	*********************************************
	Original by Alain Destexhe
	Rewritten for MyFirstNEURON by Arthur Houweling
	Post-edited and Modified by S. Romanenko
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KChTransitIA
	USEION k READ ek WRITE ik 
    RANGE gkbar, m_inf1, tau_m, h_inf, tau_h1, ik, Q10KiA, gbQ10KiA, tauMmax, Vh1m, Vs1m, Vh2m, Vs2m, CMM, VhMinf, VsMinf, tauHmax, Vh1h, Vs1h, Vh2h, Vs2h, CHH, VhHinf, VsHinf, s_celsius
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v		(mV)
	celsius (degC)
	dt		(ms)
	ek		(mV)
	gkbar = 0.009	(mho/cm2) :0.00345
	Q10KiA = 3 (1)
	gbQ10KiA = 1.6 (1)
	tauMmax = 1 
	Vh1m = 35.82 
	Vs1m = 19.69 
	Vh2m = 79.69 
	Vs2m = 12.2 
	CMM= 0.37
	VhMinf = 60 
	VsMinf = 8.5
	tauHmax = 1 
	Vh1h = 46.05 
	Vs1h = 5 
	Vh2h = 238.4 
	Vs2h = 37.45 
	CHH= 0
	VhHinf = 78 
	VsHinf = 6
	s_celsius = 0 (degC)
}

STATE {
	m1 h1
}

ASSIGNED {
	ik		(mA/cm2)
	m_inf1
	tau_m	(ms)
	h_inf
	tau_h1	(ms)
	tadj
	gbtadj
}

BREAKPOINT { 
	SOLVE states :METHOD euler
 	ik = gbtadj*gkbar * m1^2*h1 * (v-ek) :S_power 4
}

:DERIVATIVE states { 
:	evaluate_fct(v)
:
:	m1'= (m_inf1-m1) / tau_m
:	h1'= (h_inf-h1) / tau_h1
:}

PROCEDURE states() {
        evaluate_fct(v)

	m1= m1 + (1-exp(-dt/tau_m))*(m_inf1-m1)
	h1= h1 + (1-exp(-dt/tau_h1))*(h_inf-h1)
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m1 = m_inf1
    h1 = h_inf
}

PROCEDURE evaluate_fct(v(mV)) {  LOCAL a
	if (s_celsius == 0) {       
        s_celsius = celsius       
        	}
	tadj = Q10KiA^((s_celsius-20)/10) 
	gbtadj =gbQ10KiA^((s_celsius-20)/10)
	tau_m = tauMmax/((exp((v+Vh1m)/Vs1m)+exp(-(v+Vh2m)/Vs2m))+CMM) / tadj           : tauMmax = 1 Vh1m = 35.82 Vs1m = 19.69 Vh2m = 79.69 Vs2m = 12.2 CMM= 0.37
	m_inf1 = 1.0 / (1+exp(-(v+VhMinf)/VsMinf))                                      : VhMinf = 60 VsMinf = 8.5
	
	if (v<-63) {
		tau_h1 = tauHmax/((exp((v+Vh1h)/Vs1h)+exp(-(v+Vh2h)/Vs2h))+CHH) / tadj      : tauHmax = 1 Vh1h = 46.05 Vs1h = 5 Vh2h = 238.4 Vs2h = 37.45 CHH= 0
		}
	else {
		tau_h1 = 19.0/tadj
		}
	h_inf = 1.0/(1+exp((v+VhHinf)/VsHinf))                                          : VhHinf = 78 VsHinf = 6
}
UNITSON