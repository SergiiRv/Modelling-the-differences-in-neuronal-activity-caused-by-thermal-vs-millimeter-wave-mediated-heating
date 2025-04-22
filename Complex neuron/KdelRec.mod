TITLE KdelRec.mod   
 
COMMENT
 Classical Delayed Rectifier 
 Initialize this mechanism to steady-state voltage by calling
 rates_gsquid(v) from HOC, then setting m_gsquid=minf_gsquid, etc.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See hh1.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX KdelRec
        USEION k READ ek WRITE ik
        RANGE gkbar, betanbarS, alphanbarS, alphanvh, alphanvs, betanvh, betanvs, Q10KdelRec, gbQ10KdelRec, s_celsius
        :GLOBAL ninf, nexp  
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        :dt (ms)
        gkbar = .006 (mho/cm2)
        ek = -68 (mV)
		Q10KdelRec = 2.3 (1)
        gbQ10KdelRec = 1.44 (1)
		
		
		alphanbarS = 0.024 (1)
		betanbarS = 0.2 (1)
		alphanvh = 17
		alphanvs = 8
		betanvh = 48
		betanvs = 35
		s_celsius = 0 (degC)
		}
 
STATE {
       n
}
 
ASSIGNED {
        ik (mA/cm2)
        ninf 
		nexp
		dbq10
		q10
		tinc
		alpha
		beta
		sum
		
}
 
BREAKPOINT {
		if (s_celsius == 0) {       
        s_celsius = celsius       
        	}
		dbq10 = gbQ10KdelRec^((s_celsius - 20)/10)
        SOLVE states
        ik = dbq10*gkbar*n*n*(v - ek)      
}
 
UNITSOFF
 
INITIAL {
     rates(v)
     n = ninf     
}
PROCEDURE states() {           :Computes state variables n 
		rates(v)               :at the current v and dt.
        n = n + nexp*(ninf-n)
}
 
PROCEDURE rates(v) {
		if (s_celsius == 0) {       
        s_celsius = celsius       
        	}
        q10 = Q10KdelRec^((s_celsius - 20)/10)
        tinc = -dt * q10
           
        alpha = alphanbarS*vtrap(-(v-alphanvh),alphanvs)        :Vh = 17 and Vs=8
        beta = betanbarS*exp(-(v+betanvh)/betanvs)                :Vh = 48 and Vs=35
        sum = alpha + beta
        ninf = alpha/sum
        nexp = 1 - exp(tinc*sum)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

 