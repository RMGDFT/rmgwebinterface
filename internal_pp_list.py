species_list_sg15 =[ 
"AG" ,	
"AL" ,	
"AR" ,	
"AS" ,	
"AU" ,	
"BA" ,	
"BE" ,	
"BI" ,	
"BR" ,	
"B" ,	
"CA" ,	
"CD" ,	
"CL" ,	
"CO" ,	
"CR" ,	
"C" ,	
"CS" ,	
"CU" ,	
"FE" ,	
"F" ,	
"GA" ,	
"GE" ,	
"HE" ,	
"HF" ,	
"HG" ,	
"H" ,	
"IN" ,	
"IR" ,	
"I" ,	
"KR" ,	
"K" ,	
"LA" ,	
"LI" ,	
"MG" ,	
"MN" ,	
"MO" ,	
"NA" ,	
"NB" ,	
"NE" ,	
"NI" ,	
"N" ,	
"O" ,	
"OS" ,	
"PB" ,	
"PD" ,	
"PO" ,	
"P" ,	
"PT" ,	
"RB" ,	
"RE" ,	
"RH" ,	
"RU" ,	
"SB" ,	
"SC" ,	
"SE" ,	
"SI" ,	
"SN" ,	
"SR" ,	
"S" ,	
"TA" ,	
"TC" ,	
"TE" ,	
"TI" ,	
"TL" ,	
"V" ,	
"W" ,	
"XE" ,	
"Y" ,	
"ZN" ,	
"ZR" ]	

species_list_ncpp = ["AG"   ,
"AL"   ,
"AR"   ,
"AS"   ,
"AU"   ,
"BA"   ,
"BE"   ,
"BI"   ,
"BR"   ,
"B"   , 
"CA"   ,
"CD"   ,
"CL"   ,
"CO"   ,
"CR"   ,
"CS"   ,
"C"   , 
"CU"   ,
"FE"   ,
"F"   , 
"GA"   ,
"GE"   ,
"HE"   ,
"HF"   ,
"HG"   ,
"H"   , 
"IN"   ,
"IR"   ,
"I"   , 
"KR"   ,
"K"   , 
"LA"   ,
"LI"   ,
"LU"   ,
"MG"   ,
"MN"   ,
"MO"   ,
"NA"   ,
"NB"   ,
"NE"   ,
"NI"   ,
"N"   , 
"OS"   ,
"O"   , 
"PB"   ,
"PD"   ,
"PO"   ,
"P"   , 
"PT"   ,
"RB"   ,
"RE"   ,
"RH"   ,
"RN"   ,
"RU"   ,
"SB"   ,
"SC"   ,
"SE"   ,
"SI"   ,
"SN"   ,
"SR"   ,
"S"   , 
"TA"   ,
"TC"   ,
"TE"   ,
"TI"   ,
"TL"   ,
"V"   , 
"U"   , 
"W"   , 
"XE"   ,
"Y"   , 
"ZN"   ,
"ZR"   ]

species_list_uspp = [
"AG"   ,
"AL"   ,
"AS"   ,
"AU"   ,
"BA"   ,
"BE"   ,
"B"    ,
"BR"   ,
"CA"   ,
"CD"   ,
"CL"   ,
"CO"   ,
"C"    ,
"CR"   ,
"CS"   ,
"CU"   ,
"FE"   ,
"F"    ,
"GA"   ,
"GE"   ,
"HF"   ,
"HG"   ,
"H"    ,
"IN"   ,
"I"    ,
"IR"   ,
"K"    ,
"LA"   ,
"LI"   ,
"MG"   ,
"MN"   ,
"MO"   ,
"NA"   ,
"NB"   ,
"NI"   ,
"N"    ,
"O"    ,
"OS"   ,
"PB"   ,
"PD"   ,
"P"    ,
"PT"   ,
"RB"   ,
"RE"   ,
"RH"   ,
"RU"   ,
"SB"   ,
"SC"   ,
"SE"   ,
"SI"   ,
"SN"   ,
"S"    ,
"SR"   ,
"TA"   ,
"TC"   ,
"TE"   ,
"TI"   ,
"TL"   ,
"V"    ,
"W"    ,
"Y"    ,
"ZN"   ,
"ZR"   ]

for i in range(len(species_list_sg15)):
    species_list_sg15[i] = species_list_sg15[i].lower()
for i in range(len(species_list_ncpp)):
    species_list_ncpp[i] = species_list_ncpp[i].lower()
for i in range(len(species_list_uspp)):
    species_list_uspp[i] = species_list_uspp[i].lower()

num_orbitals_dict = {
"Ag" :13,
"Al" :4,
"Ar" :9,
"As" :9,
"Au" :9,
"Ba" :8,
"Be" :4,
"Bi" :13,    
"B" :4,
"Br" :4,
"Ca" :8,
"Cd" :3,
"Cl" :4,
"C" :4,
"Co" :13,
"Cr" :13,
"Cs" :8,
"Cu" :13,
"Fe" :13,
"F" :4,
"Ga" :9,
"Ge" :9,
"He" :2,
"Hf" :13,
"Hg" :13,
"H" :2,
"In" :9,
"I" :4,
"Ir" :13,
"K" :8,
"Kr" :9,
"Li" :4,
"Mg" :9,
"Mn" :13,
"Mo" :13,
"Na" :8,
"Nb" :9,
"Ne" :9,
"Ni" :13,
"N" :4,
"O" :4,
"Os" :13,
"Pb" :13,
"Pd" :13,
"P" :4,
"Po" :13,
"Pt" :13,
"Rb" :8,
"Re" :13,
"Rh" :13,
"Rn" :13,
"Ru" :13,
"Sb" :13,
"Sc" :9,
"Se" :13,
"Si" :4,
"Sn" :9,
"S" :4,
"Sr" :9,
"Ta" :9,
"Tc" :13,
"Te" :13,
"Ti" :13,
"Tl" :13,
"V" :13,
"W" :13,
"Xe" :9,
"Y" :9,
"Zn" :13,
"Zr" :9
}
