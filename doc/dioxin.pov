#include "colors.inc"
#include "skies.inc"

#declare BAS = true; 
#declare SPF = false;
#declare CST = false;
#declare TRANS = false;

#include "dioxin.inc";

background{White}
sky_sphere{S_Cloud1}

//The molecule
object{
       dioxin
       translate dioxin_center //Translate to (0|0|0) 
      }

camera {
        location <0,0,-12>
        look_at <0,0,0>
       }

light_source {<0,10,0> 
              color White
             } 

