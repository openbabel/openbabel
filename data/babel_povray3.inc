//Include for Babel
//Copyright (c) by Steffen Reith and Oliver Weichold 1996-2002
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation version 2 of the License.
//EMail: streit@streit.cc
//For chemical questions send mail to
//weichold@sas.upenn.edu

//$Id$

#include "colors.inc"
#include "metals.inc"

//Do a ball and sticks model as default (nothing defined)
#ifndef (BAS)
 #ifndef (SPF)
  #ifndef (CST)
   #declare BAS = true;
   #declare SPF = false;
   #declare CST = false;
  #end
 #end
#end

//Default: No Ball and Sticks-model (if undefined)
#ifndef (BAS)
 #declare BAS = false;
#end

//Default: No Spacefill-model (if undefined)
#ifndef (SPF)
 #declare SPF = false;
#end

//Default: No capped sticks-model (if undefined)
#ifndef (CST)
 #declare CST = false;
#end

//Is a model (BAS, SPF or CST) selected ?
#if (!((BAS | SPF) | CST))

 //Terminate parsing
 #error "You have to select BAS (default), SPF or CST" 

#end

//Select only one model
#if ((BAS & SPF) | ((BAS & CST) | (SPF & CST)))

 //Terminate parsing
 #error "You have to select only one model-type" 

#end

//Display copyright
#render "\n\n"
#render "\nbabel31.inc (C) 1996-2002 by\n"
#render "Oliver Weichold (weichold@sas.upenn.edu)\n"
#render "and Steffen Reith (streit@streit.cc)\n"
#render "\n\n"

//TRANS mode off by default (only needed for transparent textures)
#ifndef (TRANS)
#declare TRANS = false;
#end 

#if (TRANS)
#include "glass.inc"
#end

//Finishs for atoms (Nonmetals, Transmetals, Halfmetals and Metals)
#if (TRANS)

 #declare TransFinish = finish {
                                ambient 0.1
                                diffuse 0.1
                                specular 0.8
                                reflection 0.1
                                roughness 0.0003
                                phong 1 phong_size 400
                               }                      

 //Transparent finishs
 #declare MetalFinish      = finish {TransFinish}
 #declare HalfMetalFinish  = finish {TransFinish}
 #declare NonMetalFinish   = finish {TransFinish}
 #declare TransMetalFinish = finish {TransFinish}
 #declare Bond_0_Finish    = finish {TransFinish}
 #declare Bond_1_Finish    = finish {TransFinish}
 #declare Bond_2_Finish    = finish {TransFinish}
 #declare Bond_3_Finish    = finish {TransFinish}
 #declare Bond_4_Finish    = finish {TransFinish}
 #declare Bond_5_Finish    = finish {TransFinish}

 //Filter values for transparent finishs  
 #declare MetalFilter       = 1.00;
 #declare HalfMetalFilter   = 0.95;
 #declare NonMetalFilter    = 1.00;
 #declare TransMetalFilter  = 0.95;
 #declare Bond_0_Filter     = 0.00;
 #declare Bond_1_Filter     = 0.90;
 #declare Bond_2_Filter     = 0.90;
 #declare Bond_3_Filter     = 0.90;
 #declare Bond_4_Filter     = 0.90;
 #declare Bond_5_Filter     = 0.90;

#else

 //Solid finishs
 #declare MetalFinish      = finish {F_MetalA}
 #declare HalfMetalFinish  = finish {F_MetalD}
 #declare NonMetalFinish   = finish {

                            ambient 0.2
                            diffuse 0.8
                            brilliance 1.0 
                            phong 1.0      
                            phong_size 80
                            reflection 0.25

                           }
 #declare TransMetalFinish = finish {F_MetalD}
 #declare Bond_0_Finish    = finish {F_MetalC}
 #declare Bond_1_Finish    = finish {F_MetalC}
 #declare Bond_2_Finish    = finish {F_MetalC}
 #declare Bond_3_Finish    = finish {F_MetalC}
 #declare Bond_4_Finish    = finish {F_MetalC}
 #declare Bond_5_Finish    = finish {F_MetalC}

 //No additional filter (texture is not transparent)
 #declare MetalFilter       = 0.0;
 #declare HalfMetalFilter   = 0.0;
 #declare NonMetalFilter    = 0.0;
 #declare TransMetalFilter  = 0.0;
 #declare Bond_0_Filter     = 0.0;
 #declare Bond_1_Filter     = 0.0;
 #declare Bond_2_Filter     = 0.0;
 #declare Bond_3_Filter     = 0.0;
 #declare Bond_4_Filter     = 0.0;
 #declare Bond_5_Filter     = 0.0;

#end

//Thickness of Bonds (see Makro Bond_order)
#if (CST)

 // For capped-stick models you must use _one_ Bond_Thick ! 
 #declare Bond_Thick_0 = 0.1;
 #declare Bond_Thick_1 = Bond_Thick_0; 
 #declare Bond_Thick_2 = Bond_Thick_0;
 #declare Bond_Thick_3 = Bond_Thick_0;
 #declare Bond_Thick_4 = Bond_Thick_0;
 #declare Bond_Thick_5 = Bond_Thick_0;

#else

 #declare Bond_Thick_0 = 0.1;
 #declare Bond_Thick_1 = 0.1;
 #declare Bond_Thick_2 = 0.1;
 #declare Bond_Thick_3 = 0.1;
 #declare Bond_Thick_4 = 0.1;
 #declare Bond_Thick_5 = 0.1;

#end

//Colors for Bonds (see Makro Bond_order)
#declare Bond_0_Color = rgbf <1.00,0.725,0.275,Bond_0_Filter>; 
#declare Bond_1_Color = rgbf <1.00,0.725,0.275,Bond_1_Filter>;
#declare Bond_2_Color = rgbf <1.00,0.725,0.275,Bond_2_Filter>;
#declare Bond_3_Color = rgbf <1.00,0.725,0.275,Bond_3_Filter>;
#declare Bond_4_Color = rgbf <1.00,0.725,0.275,Bond_4_Filter>;
#declare Bond_5_Color = rgbf <1.00,0.725,0.275,Bond_5_Filter>;

#declare Bond_0_Color = rgbf <0.95,0.95,0.1,Bond_0_Filter>;
#declare Bond_1_Color = rgbf <0.95,0.95,0.1,Bond_1_Filter>;
#declare Bond_2_Color = rgbf <0.95,0.95,0.1,Bond_2_Filter>;
#declare Bond_3_Color = rgbf <0.95,0.95,0.1,Bond_3_Filter>;
#declare Bond_4_Color = rgbf <0.95,0.95,0.1,Bond_4_Filter>;
#declare Bond_5_Color = rgbf <0.95,0.95,0.1,Bond_5_Filter>;

//The bond declarations
#declare bond_0 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_0
                            texture {

                             pigment {color Bond_0_Color}
                             finish {Bond_0_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }
#declare bond_1 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_1
                            texture {

                             pigment {color Bond_1_Color}
                             finish {Bond_1_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }
#declare bond_2 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_2
                            texture {

                             pigment {color Bond_2_Color}
                             finish {Bond_2_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }
#declare bond_3 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_3
                            texture {

                             pigment {color Bond_3_Color}
                             finish {Bond_3_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }
#declare bond_4 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_4
                            texture {

                             pigment {color Bond_4_Color}
                             finish {Bond_4_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }
#declare bond_5 = cylinder {
                            <0,0,0>
                            <1,0,0>
                            Bond_Thick_5
                            texture {

                             pigment {color Bond_5_Color}
                             finish {Bond_5_Finish}

                            }
#if (TRANS)
			    interior{ I_Glass}
#end
                           }

//Generic definition of an atom
#if (BAS)
#render "\n\nCalculate a ball and sticks model of :"
#declare atom = sphere {<0,0,0>,0.6}
#end
#if (SPF)  
#render "\n\nCalculate a spacefill model of :"
#declare atom = sphere {<0,0,0>,1.66}
#end
#if (CST)  
#render "\n\nCalculate a capped sticks model of :"
#declare atom = sphere {<0,0,0>,1.0}
#end

//For Povray 3.1 only
#version 3.1;

//For capped sticks models we use Bond_Thick_0 for all atom radii
#if (CST)

 #declare Radius_C1 = Bond_Thick_0;
 #declare Radius_C2 = Bond_Thick_0;
 #declare Radius_C3 = Bond_Thick_0;
 #declare Radius_H  = Bond_Thick_0;
 #declare Radius_D =  Bond_Thick_0;
 #declare Radius_O  = Bond_Thick_0;
 #declare Radius_N  = Bond_Thick_0;
 #declare Radius_S  = Bond_Thick_0;
 #declare Radius_P  = Bond_Thick_0;
 #declare Radius_F  = Bond_Thick_0;
 #declare Radius_Cl = Bond_Thick_0;
 #declare Radius_Br = Bond_Thick_0;
 #declare Radius_I  = Bond_Thick_0;
		      
 #declare Radius_Si = Bond_Thick_0;
 #declare Radius_B  = Bond_Thick_0;
 #declare Radius_Ge = Bond_Thick_0;
 #declare Radius_Se = Bond_Thick_0;

 #declare Radius_Ca  = Bond_Thick_0;
 #declare Radius_Ca2 = Bond_Thick_0;
 #declare Radius_Sn  = Bond_Thick_0;
 #declare Radius_Sn2 = Bond_Thick_0;
 #declare Radius_Sn4 = Bond_Thick_0;
 #declare Radius_Pb  = Bond_Thick_0;
 #declare Radius_Pb2 = Bond_Thick_0;
 #declare Radius_Pb4 = Bond_Thick_0;
 #declare Radius_Te  = Bond_Thick_0;
 #declare Radius_Te4 = Bond_Thick_0;
 #declare Radius_Te6 = Bond_Thick_0;
 #declare Radius_Na  = Bond_Thick_0;
 #declare Radius_Na1 = Bond_Thick_0;
 #declare Radius_K   = Bond_Thick_0;
 #declare Radius_K1  = Bond_Thick_0;
 #declare Radius_Li  = Bond_Thick_0;
 #declare Radius_Li1 = Bond_Thick_0;
 #declare Radius_Al  = Bond_Thick_0;
 #declare Radius_Al3 = Bond_Thick_0;

 #declare Radius_Mn  = Bond_Thick_0;
 #declare Radius_Mn2 = Bond_Thick_0;
 #declare Radius_Mn3 = Bond_Thick_0;
 #declare Radius_Mn4 = Bond_Thick_0;
 #declare Radius_W   = Bond_Thick_0;
 #declare Radius_W4  = Bond_Thick_0;
 #declare Radius_W6  = Bond_Thick_0;
 #declare Radius_Fe  = Bond_Thick_0;
 #declare Radius_Fe2 = Bond_Thick_0;
 #declare Radius_Fe3 = Bond_Thick_0;
 #declare Radius_Pd  = Bond_Thick_0;
 #declare Radius_Pd2 = Bond_Thick_0;
 #declare Radius_Pd4 = Bond_Thick_0;
 #declare Radius_Ni  = Bond_Thick_0;
 #declare Radius_Ni2 = Bond_Thick_0;
 #declare Radius_Ni3 = Bond_Thick_0;
 #declare Radius_Cu  = Bond_Thick_0;
 #declare Radius_Cu1 = Bond_Thick_0;
 #declare Radius_Cu2 = Bond_Thick_0;
 #declare Radius_V   = Bond_Thick_0;
 #declare Radius_V2  = Bond_Thick_0;
 #declare Radius_V3  = Bond_Thick_0;
 #declare Radius_V4  = Bond_Thick_0;
 #declare Radius_V5  = Bond_Thick_0;
 #declare Radius_Zn  = Bond_Thick_0;
 #declare Radius_Zn2 = Bond_Thick_0;

#else

 //Atom radii in A as published in J. Emsley, The Elements, 
 //Clarender Press, Oxford 1991

 //All uncharged radii are kovalent. The values for Na and V are calculated 
 //from the atomic radius (10% less).
 //We use the kovalent radii for the description of metal atoms. If you want to
 //use ion radii, just replace the radius in the atom description. 
 //E.g., replace
 //
 //scale Radius_Zn
 //
 //by
 //
 //scale Radius_Zn2
 //
 //for creating a sphere with the ion radius of Zn2+.

 #declare Radius_C1 = 0.770;
 #declare Radius_C2 = 0.670;
 #declare Radius_C3 = 0.600;

 #declare Radius_H  = 0.370;
 #declare Radius_D  = 0.370;
 #declare Radius_O  = 0.660;
 #declare Radius_N  = 0.700;
 #declare Radius_S  = 1.040;
 #declare Radius_P  = 1.100;
 #declare Radius_F  = 0.580;
 #declare Radius_Cl = 0.990; 
 #declare Radius_Br = 1.142;
 #declare Radius_I  = 1.330;

 #declare Radius_Si = 1.170;
 #declare Radius_B  = 0.880;
 #declare Radius_Ge = 1.220;
 #declare Radius_Se = 1.170;

 #declare Radius_Ca  = 1.740;
 #declare Radius_Ca2 = 1.060;
 #declare Radius_Sn  = 1.400;
 #declare Radius_Sn2 = 0.930;
 #declare Radius_Sn4 = 0.740;
 #declare Radius_Pb  = 1.540;
 #declare Radius_Pb2 = 1.320;
 #declare Radius_Pb4 = 0.890;
 #declare Radius_Te  = 1.370;
 #declare Radius_Te4 = 0.970;
 #declare Radius_Te6 = 0.560;
 #declare Radius_Na  = 1.380;
 #declare Radius_Na1 = 0.980;
 #declare Radius_K   = 2.030;
 #declare Radius_K1  = 1.330;
 #declare Radius_Li  = 1.230;
 #declare Radius_Li1 = 0.780;
 #declare Radius_Al  = 1.250;
 #declare Radius_Al3 = 0.570;

 #declare Radius_Mn  = 1.170;
 #declare Radius_Mn2 = 0.910;
 #declare Radius_Mn3 = 0.700;
 #declare Radius_Mn4 = 0.520;
 #declare Radius_W   = 1.300;
 #declare Radius_W4  = 0.680;
 #declare Radius_W6  = 0.620;
 #declare Radius_Fe  = 1.165;
 #declare Radius_Fe2 = 0.820;
 #declare Radius_Fe3 = 0.670;
 #declare Radius_Pd  = 1.280;
 #declare Radius_Pd2 = 0.860;
 #declare Radius_Pd4 = 0.640;
 #declare Radius_Ni  = 1.150;
 #declare Radius_Ni2 = 0.780;
 #declare Radius_Ni3 = 0.620;
 #declare Radius_Cu  = 1.170;
 #declare Radius_Cu1 = 0.960;
 #declare Radius_Cu2 = 0.720;
 #declare Radius_V   = 1.188;
 #declare Radius_V2  = 0.720;
 #declare Radius_V3  = 0.650;
 #declare Radius_V4  = 0.610;
 #declare Radius_V5  = 0.590;
 #declare Radius_Zn  = 1.250;
 #declare Radius_Zn2 = 0.830;

#end

//Atom descriptions

//Section: Nonmetals
#declare Color_C1 = rgbf <0.00,0.00,0.00,NonMetalFilter>;
#declare Atom_C1 = object { 
                           
                    atom
                    scale Radius_C1
                    texture {
 
                     pigment {color Color_C1}
                     finish {NonMetalFinish}
 
                    }                           
#if (TRANS)
	            interior{ I_Glass}
#end

                   }
#declare Atom_C	   = Atom_C1
#declare Atom_Cac  = Atom_C1
#declare Atom_Cbl  = Atom_C1
#declare Atom_Cmi  = Atom_C1
#declare Atom_Crad = Atom_C1
#declare Atom_Cx   = Atom_C1
#declare Atom_Nr   = Atom_C1

#declare Radius_Cac  = Radius_C1;
#declare Radius_Cbl  = Radius_C1;
#declare Radius_Cmi  = Radius_C1;
#declare Radius_Crad = Radius_C1;
#declare Radius_Cx   = Radius_C1;
#declare Radius_Nr   = Radius_C1;

#declare Color_C   = Color_C1;
#declare Color_Cac  = Color_C1;
#declare Color_Cbl  = Color_C1;
#declare Color_Cmi  = Color_C1;
#declare Color_Crad = Color_C1;
#declare Color_Cx   = Color_C1;
#declare Color_Nr   = Color_C1;

#declare Color_C2 = rgbf <0.00,0.00,0.00,NonMetalFilter>;
#declare Atom_C2  = object { 

                     atom
                     scale Radius_C2
                     texture {
 
                      pigment {color Color_C2}
                      finish {NonMetalFinish}
 
                     }                           
#if (TRANS)
	             interior{ I_Glass}
#end
                    }

#declare Color_C3 = rgbf <0.00,0.00,0.00,NonMetalFilter>;
#declare Atom_C3 = object { 

                    atom
                    scale Radius_C3
                    texture {
 
                     pigment {color Color_C3}
                     finish {NonMetalFinish}

                    }                           
#if (TRANS)
	            interior{ I_Glass}
#end
                   }
#declare Atom_Car = Atom_C3
#declare Atom_Cpl = Atom_C3
#declare Atom_Cx  = Atom_C3

#declare Radius_Car = Radius_C3;
#declare Radius_Cpl = Radius_C3;
#declare Radius_Cx  = Radius_C3;

#declare Color_Car = Color_C3;
#declare Color_Cpl = Color_C3;
#declare Color_Cx  = Color_C3;

#declare Color_H = rgbf <1.00,1.00,1.00,NonMetalFilter>;
#declare Atom_H = object { 
                   
                   atom
                   scale Radius_H
                   texture {
 
                    pigment {color Color_H}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_HC  = Atom_H
#declare Atom_HO  = Atom_H
#declare Atom_HN  = Atom_H
#declare Atom_Hpl = Atom_H

#declare Radius_HC  = Radius_H;
#declare Radius_HO  = Radius_H;
#declare Radius_HN  = Radius_H;
#declare Radius_Hpl = Radius_H;

#declare Color_HC  = Color_H;
#declare Color_HO  = Color_H;
#declare Color_HN  = Color_H;
#declare Color_Hpl = Color_H;

#declare Color_D = rgbf <0.85,0.85,0.95,NonMetalFilter>;
#declare Atom_D = object { 
                          
                   atom
                   scale Radius_D
                   texture {
 
                    pigment {color Color_D}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_DC = Atom_D

#declare Radius_DC = Radius_D;

#declare Color_DC = Color_D;

#declare Color_O = rgbf <1.00,0.00,0.00,NonMetalFilter>;
#declare Atom_O = object { 

                   atom
                   scale Radius_O
                   texture {

                    pigment {color Color_O}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_O3   = Atom_O
#declare Atom_O2   = Atom_O
#declare Atom_Omi  = Atom_O
#declare Atom_OCO2 = Atom_O
#declare Atom_Oany = Atom_O
#declare Atom_Oco2 = Atom_O

#declare Radius_O3   = Radius_O;
#declare Radius_O2   = Radius_O;
#declare Radius_Omi  = Radius_O;
#declare Radius_OCO2 = Radius_O;
#declare Radius_Oany = Radius_O;
#declare Radius_Oco2 = Radius_O;

#declare Color_O3   = Color_O;
#declare Color_O2   = Color_O;
#declare Color_Omi  = Color_O;
#declare Color_OCO2 = Color_O;
#declare Color_Oany = Color_O;
#declare Color_Oco2 = Color_O;

#declare Color_N = rgbf <0.14,0.42,0.56,NonMetalFilter>;
#declare Atom_N = object { 
     
                   atom
                   scale Radius_N
                   texture {
                         
                    pigment {color Color_N}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_NPL3 = Atom_N
#declare Atom_N2   = Atom_N
#declare Atom_N4   = Atom_N
#declare Atom_N3   = Atom_N
#declare Atom_Npl  = Atom_N
#declare Atom_N2pl = Atom_N
#declare Atom_Naz  = Atom_N
#declare Atom_N1   = Atom_N
#declare Atom_Nox  = Atom_N
#declare Atom_Ntr  = Atom_N
#declare Atom_NC   = Atom_N
#declare Atom_Nar  = Atom_N
#declare Atom_Nam  = Atom_N
#declare Atom_Ngpl = Atom_N
#declare Atom_N3pl = Atom_N
#declare Atom_Nany = Atom_N

#declare Radius_NPL3 = Radius_N;
#declare Radius_N2   = Radius_N;
#declare Radius_N4   = Radius_N;
#declare Radius_N3   = Radius_N;
#declare Radius_Npl  = Radius_N;
#declare Radius_N2pl = Radius_N;
#declare Radius_Naz  = Radius_N;
#declare Radius_N1   = Radius_N;
#declare Radius_Nox  = Radius_N;
#declare Radius_Ntr  = Radius_N;
#declare Radius_NC   = Radius_N;
#declare Radius_Nar  = Radius_N;
#declare Radius_Nam  = Radius_N;
#declare Radius_Ngpl = Radius_N;
#declare Radius_N3pl = Radius_N;
#declare Radius_Nany = Radius_N;

#declare Color_NPL3 = Color_N;
#declare Color_N2   = Color_N;
#declare Color_N4   = Color_N;
#declare Color_N3   = Color_N;
#declare Color_Npl  = Color_N;
#declare Color_N2pl = Color_N;
#declare Color_Naz  = Color_N;
#declare Color_N1   = Color_N;
#declare Color_Nox  = Color_N;
#declare Color_Ntr  = Color_N;
#declare Color_NC   = Color_N;
#declare Color_Nar  = Color_N;
#declare Color_Nam  = Color_N;
#declare Color_Ngpl = Color_N;
#declare Color_N3pl = Color_N;
#declare Color_Nany = Color_N;

#declare Color_F = rgbf<1.00,1.00,0.88,NonMetalFilter>;
#declare Atom_F = object {

                   atom
                   scale Radius_F
                   texture {
                   
                    pigment {color Color_F}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }

#declare Color_Cl = rgbf <0.20,0.80,0.20,NonMetalFilter>;
#declare Atom_Cl = object {

                    atom
                    scale Radius_Cl
                    texture {
                
                     pigment {color Color_Cl}
                     finish {NonMetalFinish}
 
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Br = rgbf <0.65,0.16,0.16,NonMetalFilter>;
#declare Atom_Br = object {

                    atom
                    scale Radius_Br
                    texture {
                    
                     pigment {color Color_Br}
                     finish {NonMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_I = rgbf<0.58,0.44,0.86,NonMetalFilter>;
#declare Atom_I = object {

                   atom
                   scale Radius_I
                   texture {
                  
                    pigment {color Color_I}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }

#declare Color_P = rgbf <1.00,0.50,0.00,NonMetalFilter>;
#declare Atom_P = object {

                   atom
                   scale Radius_P
                   texture {
                   
                    pigment {color Color_P}
                    finish{NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_Pac  = Atom_P
#declare Atom_P3   = Atom_P
#declare Atom_Pox  = Atom_P
#declare Atom_P3pl = Atom_P

#declare Radius_Pac  = Radius_P;
#declare Radius_P3   = Radius_P;
#declare Radius_Pox  = Radius_P;
#declare Radius_P3pl = Radius_P;

#declare Color_Pac  = Color_P;
#declare Color_P3   = Color_P;
#declare Color_Pox  = Color_P;
#declare Color_P3pl = Color_P;

#declare Color_S = rgbf <1.00,1.00,0.00,NonMetalFilter>;
#declare Atom_S = object {

                   atom
                   scale Radius_S
                   texture {
                         
                    pigment {color Color_S}
                    finish {NonMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }
#declare Atom_S3pl = Atom_S
#declare Atom_S3   = Atom_S
#declare Atom_S2   = Atom_S
#declare Atom_Sac  = Atom_S
#declare Atom_Sox  = Atom_S
#declare Atom_Sany = Atom_S
#declare Atom_So   = Atom_S

#declare Radius_S3pl = Radius_S;
#declare Radius_S3   = Radius_S;
#declare Radius_S2   = Radius_S;
#declare Radius_Sac  = Radius_S;
#declare Radius_Sox  = Radius_S;
#declare Radius_Sany = Radius_S;
#declare Radius_So   = Radius_S;

#declare Color_S3pl = Color_S;
#declare Color_S3   = Color_S;
#declare Color_S2   = Color_S;
#declare Color_Sac  = Color_S;
#declare Color_Sox  = Color_S;
#declare Color_Sany = Color_S;
#declare Color_So   = Color_S;

//A special color for halfmetals
#declare HalfMetalColor = rgbf <0.6,0.6,0.85,HalfMetalFilter>; 

//Section: HalfMetals
#declare Color_B = rgbf <0,1.0,0,HalfMetalFilter>;
#declare Atom_B = object {

                   atom
                   scale Radius_B
                   texture {
                   
                    pigment {color Color_B}
                    finish {HalfMetalFinish}
                   
                   }
#if (TRANS)
                   interior{ I_Glass}
#end
                  }
#declare Atom_Bac = Atom_B
#declare Atom_Box = Atom_B

#declare Radius_Bac = Radius_B;
#declare Radius_Box = Radius_B;

#declare Color_Bac = Color_B;
#declare Color_Box = Color_B;

#declare Color_Si = HalfMetalColor;
#declare Atom_Si = object {
 
                    atom
                    scale Radius_Si
                    texture {

                     pigment {color Color_Si}
                     finish {HalfMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }


#declare Color_Ge = HalfMetalColor;
#declare Atom_Ge = object { 

                    atom
                    scale Radius_Ge
                    texture {

                     pigment {color Color_Ge}
                     finish {HalfMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }
#declare Color_Se = HalfMetalColor;
#declare Atom_Se = object { 

                    atom
                    scale Radius_Se
                    texture {

                     pigment {color Color_Se}
                     finish {HalfMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

//Metals
#declare Color_Metal = rgbf <0.94,0.93,0.90,MetalFilter>;

#declare Color_Te = rgbf <0.6,0.6,0.6,MetalFilter>;
#declare Atom_Te = object { 

                    atom
                    scale Radius_Te
                    texture {

                     pigment {color Color_Te}
                     finish {MetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Li = Color_Metal;
#declare Atom_Li = object { 

                    atom
                    scale Radius_Li                          
                    texture {
                     
                     pigment {color Color_Li}
                     finish {MetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Na = Color_Metal;
#declare Atom_Na = object { 

                    atom
                    scale Radius_Na
                    texture {
                      
                     pigment {color Color_Na}
                     finish {MetalFinish}
   
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_K = Color_Metal;
#declare Atom_K = object { 

                   atom
                   scale Radius_K
                   texture {
                           
                    pigment {color Color_K}
                    finish {MetalFinish}

                   }
#if (TRANS)
	            interior{ I_Glass}
#end
                  }

#declare Color_Ca = Color_Metal;
#declare Atom_Ca = object {
 
                    atom
                    scale Radius_Ca
                    texture {
                     
                     pigment {color Color_Ca}
                     finish {MetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Al = Color_Metal;
#declare Atom_Al = object { 

                    atom
                    scale Radius_Al
                    texture {
                    
                     pigment {color Color_Al}
                     finish {MetalFinish}
 
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Sn = Color_Metal;
#declare Atom_Sn = object { 

                    atom
                    scale Radius_Sn
                    texture {
                    
                     pigment {color Color_Sn}
                     finish {MetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Pb = Color_Metal;
#declare Atom_Pb = object { 

                    atom
                    scale Radius_Pb
                    texture {
                     
                     pigment {color Color_Pb}
                     finish {MetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

//TransMetals
#declare Color_TransMetal = rgbf <0.94,0.93,0.90,TransMetalFilter>;

#declare Color_W = Color_TransMetal;
#declare Atom_W = object { 

                   atom
                   scale Radius_W
                   texture {

                    pigment {color Color_W}
                    finish {TransMetalFinish}

                   }
#if (TRANS)
	           interior{ I_Glass}
#end
                  }

#declare Color_Pd = Color_TransMetal;
#declare Atom_Pd = object { 

                    atom
                    scale Radius_Pd
                    texture {
  
                     pigment {color Color_Pd}
                     finish {TransMetalFinish}
 
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Ni = Color_TransMetal;
#declare Atom_Ni = object { 

                    atom
                    scale Radius_Ni
                    texture {

                     pigment {color Color_Ni}
                     finish {TransMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Cu = rgbf <0.50,0.25,0.15,TransMetalFilter>;
#declare Atom_Cu = object { 

                    atom
                    scale Radius_Cu

                    texture {

                     pigment {color Color_Cu}
                     finish {TransMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Mn = Color_TransMetal;
#declare Atom_Mn = object { 

                    atom
                    scale Radius_Mn
                    texture {

                     pigment {color Color_Mn}
                     finish {TransMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }

#declare Color_Fe = Color_TransMetal;
#declare Atom_Fe = object { 

                    atom
                    scale Radius_Fe

                    texture {

                     pigment {color Color_Fe}
                     finish {TransMetalFinish}

                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   }


#declare Color_V = Color_TransMetal;
#declare Atom_V = object {
                    
                   atom
                   scale Radius_V
                   texture {
                     
                    pigment {color Color_V}
                    finish {TransMetalFinish}
                     
                   } 
#if (TRANS)
	           interior{ I_Glass}
#end
                  }

#declare Atom_V2 = object {
                    
                   atom
                   scale Radius_V2
                   texture {
                     
                    pigment {color Color_V}
                    finish {TransMetalFinish}
                     
                   } 
#if (TRANS)
	           interior{ I_Glass}
#end
                  }  

#declare Atom_V3 = object {
                    
                   atom
                   scale Radius_V3
                   texture {
                     
                    pigment {color Color_V}
                    finish {TransMetalFinish}
                     
                   } 
#if (TRANS)
	           interior{ I_Glass}
#end
                  }  

#declare Atom_V4 = object {
                    
                   atom
                   scale Radius_V4
                   texture {
                     
                    pigment {color Color_V}
                    finish {TransMetalFinish}
                     
                   } 
#if (TRANS)
	           interior{ I_Glass}
#end
                  } 

#declare Atom_V5 = object {
                    
                   atom
                   scale Radius_V5
                   texture {
                     
                    pigment {color Color_V}
                    finish {TransMetalFinish}
                     
                   } 
#if (TRANS)
	           interior{ I_Glass}
#end
                  }  


#declare Color_Zn = Color_TransMetal;
#declare Atom_Zn = object {
                    
                    atom
                    scale Radius_Zn
                    texture {
                     
                     pigment {color Color_Zn} 
                     finish {TransMetalFinish}
                     
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   } 

#declare Atom_Zn2 = object {
                    
                    atom
                    scale Radius_Zn2
                    texture {
                     
                     pigment {color Color_Zn} 
                     finish {TransMetalFinish}
                     
                    }
#if (TRANS)
	            interior{ I_Glass}
#end
                   } 
