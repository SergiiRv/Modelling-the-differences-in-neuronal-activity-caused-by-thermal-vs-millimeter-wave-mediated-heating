TITLE leak channels

COMMENT
    additional current in the total current eqution caused by alterations in membrane capasitance 
	from Hines comment on: hines » Sat Apr 09, 2016 
	d(c*v)/dt + i_ion = i_stim
	or
	c*dv/dt + (dc/dt)*v + i_ion = i_stim
	done by S. Romanenko
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX leak
		NONSPECIFIC_CURRENT il
		RANGE gl, el, KCK, cml
		

}

PARAMETER { 
        v (mV)
        celsius (degC)
        dt (ms)
		gl = 0.0005 (mho/cm2)
        el = -48 (mV)
		KCK =0.001 (1)
		cml
		TK
		LeakTermCoef = 1.24 (1) :1.44
		
}
 
STATE {
    cOld    
}

:LOCAL cOld

ASSIGNED {
 
        il (mA/cm2)
		i1 (mA/cm2)
		i2 (mA/cm2)
		
		:cOld 
		
}

INITIAL {

	cOld = cml
}
 
BREAKPOINT {
		TK = LeakTermCoef^((celsius - 20)/10)
        i1 = gl*(v - el/TK)
		i2 = v*(cml-cOld)/dt*KCK
		il = (i1+i2)*TK
		SOLVE states
}

UNITSOFF
PROCEDURE states() {
	
    cOld = cml
}
UNITSON 



