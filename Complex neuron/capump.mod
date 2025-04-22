NEURON {       
   SUFFIX capump       
   USEION ca READ cai WRITE ica       
   RANGE ica, Ksrink, s_celsius
   }       
       
UNITS {       
   (uM) = (micro/liter)       
   (mM) = (milli/liter)       
   (mA) = (milliamp)       
}       
       
INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }       
       
PARAMETER {       
   final_conc = 0.0001 (mM) 
   tau = 150 (ms) 
   celsius (degC)
   Ksrink = 1 (1)
   s_celsius = 0 (degC)
}       
       
INITIAL { 
cai = 0.0001 
} 
 
ASSIGNED {       
   ica (mA/cm2)       
   cai (mM)       
}       
       
LOCAL Q       

UNITSOFF 
       
BREAKPOINT {       
   if (s_celsius == 0) {       
      s_celsius = celsius       
      Q = 4^((s_celsius - 20)/10)       
		}               
   ica =  Q*(cai-final_conc)/tau*(Ksrink) 
}       
UNITSON  
  