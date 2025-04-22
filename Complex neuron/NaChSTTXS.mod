TITLE NaChSTTXS is a transient ttx-sensitive Na+ current, 
: This current likely consists of NaV1.7 with a little
: NaV1.6 and NaV1.1 mixed in (TMM).
 
COMMENT
        Modification of Na current from Rz and P cells into more diverse and realistic model
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX NaChSTTS
        USEION na READ ena WRITE ina
        RANGE gnabar, amQ10NaTTXR, 	bmQ10NaTTXR, ahQ10NaTTXR, bhQ10NaTTXR, gbarQ10NaTTXR, A_am, B_am, C_am, A_ah, B_ah, C_ah, A_bm, B_bm, C_bm, A_bh, B_bh, C_bh, Vh_m, Vs_m, Vh_h, Vs_h, s_celsius
        :GLOBAL  alpham, alphah, betam, betah 
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
	v (mV)
    celsius    (degC)
    dt (ms)
    gnabar = .21455 (mho/cm2)  : (0.35*0.613) 0.35*0.45
    ena = 60 (mV)
	amQ10NaTTXR = 2.2 (1)
	bmQ10NaTTXR = 14  (1)
	ahQ10NaTTXR = 2.2 (1)
	bhQ10NaTTXR = 2.2 (1)
	gbarQ10NaTTXR = 1.95 (1)
	s_celsius = 0 (degC)

	A_am = 0.03 (/ms)       : 17.235 (/ms) : A for alpha m
	B_am = 28 (mV)          : 7.58 (mV)
	C_am = 15 (mV)          : -11.47 (mV)

	A_ah = 0.045 (/ms)      : 0.23688 (/ms) : A for alpha h
	B_ah = 58 (mV)          : 115 (mV)
	C_ah = 18 (mV)          : 46.33 (mV)

	A_bm = 2.7 (mV)         : 17.235 (/ms) : A for beta m
	B_bm = 53 (mV)          : 66.2 (mV)
	C_bm = 18 (mV)          : 19.8 (mV)

	A_bh = 0.72 (/ms)       : 10.8 (/ms)   : A for beta h
	B_bh = 23 (mV)          : -11.8 (mV)
	C_bh = 14 (mV)          : -11.998 (mV)
	
	Vh_m = 23 (mV)
	Vs_m = 11.142 (mV)
	Vh_h = 61.5 (mV)
	Vs_h = 8.7 (mV)
}
 
STATE {
        m h
}
 
ASSIGNED {
        ina (mA/cm2)
        amq10 
		bmq10
		ahq10 
		bhq10
		gbq10
		m_inf
		h_inf
		m_tau (ms)
		h_tau (ms)
}
 
BREAKPOINT {
		
        SOLVE states METHOD derivimplicit  :cnexp
        ina = gnabar*(m^4)*h*(v - ena)
        }
 
UNITSOFF
 
INITIAL {
	
	M_inf(v)
	H_inf(v)
		
    m = m_inf
    h = h_inf
        
}

DERIVATIVE states {
    M_inf(v)
	H_inf(v)
	M_tau(v)
	H_tau(v)
	
	m' = (m_inf - m)/m_tau
	h' = (h_inf - h)/h_tau
}

:PROCEDURE states() {
:   M_inf(v)
:	H_inf(v)
:	M_tau(v)
:	H_tau(v)

:	m= m + (1-exp(-dt/m_tau))*(m_inf-m)
:	h= h + (1-exp(-dt/h_tau))*(h_inf-h)
:}

UNITSOFF

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
FUNCTION Q10(celsius) {
		if (s_celsius == 0) {       
        s_celsius = celsius       
        	}
	amq10 = amQ10NaTTXR^((s_celsius - 20)/10)
	bmq10 = bmQ10NaTTXR^((s_celsius - 20)/10)
	ahq10 = ahQ10NaTTXR^((s_celsius - 20)/10)
	bhq10 = bhQ10NaTTXR^((s_celsius - 20)/10)
	gbq10 = gbarQ10NaTTXR^((s_celsius - 20)/10)
}
		
:_______________________________________________________________________________________

FUNCTION M_inf(Vm) {
	m_inf = 1/(1+exp(-(Vm+Vh_m)/Vs_m))
}

FUNCTION H_inf(Vm) {
	h_inf = 1/(1+exp((Vm+Vh_h)/Vs_h))
}

FUNCTION M_tau(Vm) {
	Q10(celsius)
	m_tau = 1/((A_am * vtrap(-(v+B_am),C_am))*amq10 + (A_bm * exp(-(v+B_bm)/C_bm))*bmq10)
}

FUNCTION H_tau(Vm) {
	Q10(celsius)
	h_tau = 1/((A_ah * exp(-(v+B_ah)/C_ah))*ahq10 + (A_bh / (exp(-(v+B_bh)/C_bh) + 1))*bhq10)
}
UNITSON

 