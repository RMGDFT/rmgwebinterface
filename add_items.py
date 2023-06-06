
import streamlit as st
from internal_pp_list import *
import math
def add_pseudo(species):  
    expand_ = st.expander("PSEUDOPOTENTIAL")
    with expand_:
        st.subheader('Pseudopotentials')
        cstart, col1 = st.columns([0.1,1])
        pseudo_dict={}
        internal_pp = True 
        pp_type = st.radio("Select pseudopotentials", ["SG15(NC)", "Pseudo Dojo(NC)", "GBRV-1.5(US)", "SelectFiles"], 
            help = """
            ONCV: build-in normal conserving Haman pseudopotentials
            GBRV: build-in ultra-soft pseudopotentials
            SelectFiles: need to be in upf2 or xml format""")
        pseudolines = ""  
        if pp_type == "SG15(NC)":
            print(species_list_sg15)
            for sp in species:
                if sp.lower() not in species_list_sg15:
                    internal_pp = False
                    st.markdown("no internal SG15 pseudopotential available for specie " + sp )
                    st.markdown("select pseudopotentials by yourself ")
                    break;

            pseudolines = '#******* Pseudopotentials *******   \n'
            pseudolines += 'internal_pseudo_type = "sg15"  \n'
            pseudolines += '#use Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotenitals  \n'
            pseudolines += '#those pseudopotentials are built in with RMG  \n'
        if pp_type == "GBRV-1.5(US)":
            print(species_list_uspp)
            for sp in species:
                if sp.lower() not in species_list_uspp:
                    internal_pp = False
                    st.markdown("no internal GBRV pseudopotential available for specie " + sp )
                    st.markdown("select pseudopotentials by yourself ")
                    break;
            pseudolines = '#******* Pseudopotentials *******   \n'
            pseudolines += 'internal_pseudo_type = "ultrasoft"  \n'
            pseudolines += '#use Vanderbilt ultrasoft (GBRV) pseudopotenitals  \n'
            pseudolines += '#those pseudopotentials are built in with RMG   \n'
        if pp_type == "Pseudo Dojo(NC)":
            for sp in species:
                if sp.lower() not in species_list_ncpp:
                    internal_pp = False
                    st.markdown("no internal Pseudo Dojo pseudopotential available for specie " + sp )
                    st.markdown("select pseudopotentials by yourself ")
                    break;

            pseudolines = '#******* Pseudopotentials *******   \n'
            pseudolines += 'internal_pseudo_type = "nc_accuracy"  \n'
            pseudolines += '#use Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotenitals  \n'
            pseudolines += '#those pseudopotentials are built in with RMG  \n'

        cstart, col1 = st.columns([0.2,1])
        with col1:
            if pp_type == "SelectFiles" or not internal_pp:
               pseudo_dir = st.text_input("pseudopotential file directory", value="./", on_change=None)
               pseudolines = '#******* Pseudopotentials *******   \n'
               pseudolines += 'pseudo_dir = "' + pseudo_dir + '"  \n'
               pseudolines += 'pseudopotential = "  \n' 
               for sp in species:
                   pseudo_dict[sp] = st.text_input(sp+":", value=sp+".UPF", on_change=None)
                   pseudolines += sp + '   ' + pseudo_dict[sp] +'  \n'

               pseudolines += '"  \n' 

        write_pseudopotential_plots = st.checkbox("flag to write pseudopotential plots", False)
        pseudolines += 'write_pseudopotential_plots ="%s"  \n'%str(write_pseudopotential_plots)

        tune_pp = st.checkbox("tune the localization of pseudopotenitals ", False)
        if tune_pp:
            cstart, col1, col2 = st.columns([0.1,1,1])
            localize_projectors = col1.checkbox("non-local projector localization", value = True, 
                    help = "false: non-local projectors will spread in whole space, similar to plane wave codes")
            localize_localpp = col2.checkbox("local potential localization", value = True, 
                    help = "fasle: pseudopotential's local part spreads in the whole space")
            max_nlradius = col1.number_input("max radius of non-local projector", 100.0)
            min_nlradius = col2.number_input("max radius of non-local projector", 2.0)
            max_qradius = col1.number_input("max radius of q functions in Ultrasoft PP", 100.0)
            min_qradius = col2.number_input("min radius of q functions in Ultrasoft PP", 2.0)

            pseudolines += 'localize_localpp ="%s"  \n'%str(localize_localpp)
            pseudolines += 'localize_projectors ="%s"  \n'%str(localize_projectors)
            pseudolines += 'max_nlradius ="%f"  \n'%max_nlradius
            pseudolines += 'min_nlradius ="%f"  \n'%min_nlradius
            pseudolines += 'max_qradius ="%f"  \n'%max_qradius
            pseudolines += 'min_qradius ="%f"  \n'%min_qradius

    pseudolines += '  \n'
    return pseudolines    
def add_kpoint_mesh(cell):
    cs, col1, col2, col3 = st.columns([0.2,1,1,1])
    with col1:
        k_delta = st.number_input("kdelta(2PI/bohr)", value=0.2, help ="use kdelta to estimate kmesh")
    recip_lat= cell.reciprocal_latticevectors()
    for i in range(3):
        for j in range(3):
            recip_lat[i][j] = recip_lat[i][j]/cell.lengthscale
    kmesh_init = [max(1, int(round(b.length()/k_delta))) for b in recip_lat]

    kmesh_init_str =""
    for i in range(3):
        kmesh_init_str += str(kmesh_init[i]) + "  "
    with col2:
        kmesh_str = st.text_input("kpoint mesh", value=kmesh_init_str)
    with col3:
        kshift_str = st.text_input("kpoint shift", value="0 0 0", help="0 0 0 including Gamma point")
    kpointlines = 'kpoint_mesh="' + kmesh_str +'"  \n'    
    kpointlines += 'kpoint_is_shift="' + kshift_str +'"   \n'    
    with col2:
        kdist = st.text_input("kpoints distribution", value="1", 
                help =" control the parallel over kpoints")
    kpointlines +='kpoint_distribution = "' + kdist +'"   \n'    
    return kpointlines

          
def add_kpoint_text():
    cs, col1 = st.columns([0.2,1])
    kp_list_str=col1.text_area("K point list in unit of reciprocal lattice vectors and its weight", "0.0  0.0  0.0  1.0")
    kp_list = kp_list_str.split("\n")
    kpoints = ""
    kpointlines = 'kpoint_mesh = "-1 1 1"  \n'
    num_kpt = 0
    for kp in kp_list:
        if(len(kp.split()) ==4):
            num_kpt+=1
            kpoints += kp + '  \n'
    kpointlines += 'kpoints = "  \n'
    kpointlines += kpoints
    kpointlines += '"  \n'
    col1.markdown(kpointlines)
    if num_kpt == 0:
        st.markdown("kpoint list need to be kx, ky, kz, weight, 4 numbers in a row")
    return kpointlines

def add_kbandstr_lines():
    cs, col1 = st.columns([0.2,1])
    kp_list_str=col1.text_area("special lines for band structure calculation",
       '''   0.0   0.0   0.0   0  G
   0.5   0.0   0.0   20 X''', help ="kx, ky, kz, num, symbol, in unit of reciprocal lattice vector, num: number of kpoinks to previous special kpoint. symbol for plot")
    kpointlines = 'kpoints_bandstructure = "  \n'
    kp_list = kp_list_str.split("\n")
    for kp in kp_list:
        if(len(kp.split()) == 5):
            kpointlines += kp + "  \n"
        else:
            st.markdown("format is wrong for kpoint lines for band structure")

    kpointlines += '"  \n'
    col1.markdown(kpointlines)
    return kpointlines


def add_kpoints(cell):
    expand_ = st.expander("K POINTS")
    kpointlines = '#********* K POINT SETUP *********  \n'   
    with expand_:
        kp_method = st.radio("use gamma point, a mesh or a list", ["gamma", "use mesh", "use list"])
        if kp_method == "gamma":
            kpointlines += 'kpoint_mesh = "1 1 1"  \n'
            kpointlines += 'kpoint_is_shift = "0 0 0"  \n'
        elif kp_method == "use mesh":
            kpointlines += add_kpoint_mesh(cell)
        if kp_method == "use list":
            kpointlines += add_kpoint_text()
        kp_bandstr = st.radio("kpoints for band structure", ["None", "use special lines", "use list"])
        if kp_method != "use list" and kp_bandstr == "use list":
            kpointlines += add_kpoint_text()
        if kp_bandstr == "use special lines":    
            kpointlines += add_kbandstr_lines()
    kpointlines += '  \n'
    return kpointlines
def add_control():
    expand_ = st.expander("CONTROL OPTIONS")
    extra_lines =""
    ctrl_lines =""
    with expand_:
        start_mode = st.radio("start mode", 
                ["LCAO Start", "Restart From File", "Random Start",
                  "FIREBALL Start", "Gaussian Start",
                    "Modified LCAO Start"])
        calculation_mode= st.radio("calculation mode", 
                ["Quench Electrons  ", 
                    "Relax Structure  ", 
                    "Constant Volume And Energy  ",
                    "Constant Temperature And Energy   ",
                    "Constant Pressure And Energy  ", 
                    "Plot  ", 
                    "Psi Plot  ", 
                    "Band Structure Only  ", 
                    "NEB Relax  ",
                    "Dimer Relax  ",
                    "TDDFT ",
                    "STM",
                    "NSCF"
                    ])
        if calculation_mode == "TDDFT ":
            cs, col1,col2 = st.columns([0.1,1,1,1])

            restart_tddft = col1.checkbox("restart TDDFT?", False)
            tddft_mode = col2.radio("TDDFT mode", ["electric field", "point charge"])
            tddft_steps = col1.number_input("number tddft steps", 2000)
            tddft_time_step = col2.number_input("tddft time step in atomic unit", 0.2) 

            extra_lines += 'restart_tddft = "%s"  \n'%str(restart_tddft) 
            extra_lines += 'tddft_mode = "%s"  \n'%tddft_mode
            extra_lines += 'tddft_steps = "%d"  \n'%tddft_steps
            extra_lines += 'tddft_time_step = "%f"  \n'%tddft_time_step

            if tddft_mode == "electric field":
                electric_field_magnitude = col1.number_input("E field value", 0.001)
                electric_field_vector = col2.text_input("E field direction", "1  0  0")
                extra_lines += 'electric_field_magnitude = "%f"  \n'%electric_field_magnitude
                extra_lines += 'electric_field_vector = "%s"  \n'%electric_field_vector 
            else:
                tddft_qpos = col1.text_input("point charge position", "0.0  0.0  0.0")
                tddft_qgau = col2.number_input("Gaussian for point charge", 0.1)
                extra_lines += 'tddft_qpos = "%s"  \n'%tddft_qpos
                extra_lines += 'tddft_qgau = "%f"  \n'%tddft_qgau

        subdiag_driver = st.radio("diagonalizatoin libs",
                ["auto", "lapack", "scalapack", "magma", 
                 "cusolver", "elpa", "rocsolver"])
        if subdiag_driver == "scalapack":
            blk_size = st.number_input("block dim for scalapack", 64)
            extgra_lines += 'scalapack_block_factor = "%d"  \n'%blk_size
        kohn_sham_solver=st.radio("kohn_sham_solver", ["davidson", "multigrid"],
               help="Davidson is prefered for a small system and multigrid for a large system")


        more_ctrl = st.checkbox("check the box for more control options", False)
        if more_ctrl:
            if kohn_sham_solver == "davidson":
                cs, col1,col2,col3 = st.columns([0.1,1,1,1])
                davidson_multiplier = col1.number_input("davidson_multiplier",0)
                davidson_max_steps  = col2.number_input("davidson_max_steps", 8)
                davidson_premg      = col3.number_input("davidson_premg", 4, help = "number of multigrid steps before davidson") 
                extra_lines += 'davidson_multiplier = "%d"  \n'%davidson_multiplier
                extra_lines += 'davidson_max_steps  = "%d"  \n'%davidson_max_steps
                extra_lines += 'davidson_premg      = "%d"  \n'%davidson_premg
            else:
                cs, col1,col2,col3= st.columns([0.1,1,1,1])
                kohn_sham_mg_levels = col1.number_input("kohn_sham_mg_levels", -1, 
                    help = "negative: code determines by automatically")
                kohn_sham_pre_smoothing = col2.number_input("kohn_sham_pre_smoothing", 2)
                kohn_sham_post_smoothing = col3.number_input("kohn_sham_post_smoothing", 2)
                kohn_sham_mucycles = col1.number_input("kohn_sham_mucycles", 2)
                kohn_sham_coarse_time_step = col2.number_input("kohn_sham_coarse_time_step", 1.0)
                kohn_sham_time_step = col3.number_input("kohn_sham_time_step", 0.66)
                kohn_sham_mg_timestep = col1.number_input("kohn_sham_mg_timestep", 0.66)

                extra_lines += 'kohn_sham_mg_levels = "%d"  \n'%kohn_sham_mg_levels
                extra_lines += 'kohn_sham_pre_smoothing = "%d"  \n'%kohn_sham_pre_smoothing
                extra_lines += 'kohn_sham_post_smoothing = "%d"  \n'%kohn_sham_post_smoothing
                extra_lines += 'kohn_sham_mucycles = "%d"  \n'%kohn_sham_mucycles
                extra_lines += 'kohn_sham_coarse_time_step = "%f"  \n'%kohn_sham_coarse_time_step
                extra_lines += 'kohn_sham_time_step = "%f"  \n'%kohn_sham_time_step
                extra_lines += 'kohn_sham_mg_timestep = "%f"  \n'%kohn_sham_mg_timestep
            poisson_solver = st.radio("Poisson Solver",
                    ["pfft", "multigrid"])
            if poisson_solver == "multigrid":
                cs, col1,col2,col3 = st.columns([0.1,1,1,1])
                poisson_mg_levels = col1.number_input("poisson_mg_levels", -1)
                poisson_pre_smoothing = col2.number_input(" poisson_pre_smoothing", 2)
                poisson_post_smoothing = col3.number_input(" poisson_post_smoothing", 1)
                poisson_mucycles = col1.number_input(" poisson_mucycles", 3)
                poisson_finest_time_step = col2.number_input(" poisson_finest_time_step", 1.0)
                poisson_coarse_time_step = col3.number_input(" poisson_coarse_time_step", 0.8)
                poisson_coarsest_steps = col1.number_input(" poisson_coarsest_steps", 25)
                hartree_max_sweeps = col2.number_input(" hartree_max_sweeps", 10)
                hartree_min_sweeps = col3.number_input(" hartree_min_sweeps", 5)
                extra_lines += 'poisson_mg_levels = "%d"  \n'%poisson_mg_levels
                extra_lines += 'poisson_pre_smoothing = "%d"  \n'% poisson_pre_smoothing
                extra_lines += 'poisson_post_smoothing = "%d"  \n'% poisson_post_smoothing
                extra_lines += 'poisson_mucycles = "%d"  \n'% poisson_mucycles
                extra_lines += 'poisson_finest_time_step = "%f"  \n'% poisson_finest_time_step
                extra_lines += 'poisson_coarse_time_step = "%f"  \n'% poisson_coarse_time_step
                extra_lines += 'poisson_coarsest_steps = "%d"  \n'% poisson_coarsest_steps
                extra_lines += 'hartree_max_sweeps = "%d"  \n'% hartree_max_sweeps
                extra_lines += 'hartree_min_sweeps = "%d"  \n'% hartree_min_sweeps


            relax_mass = st.radio("mass for atoms", ["Atomic", "Equal"], 
                    help="equal mas for fast relax may help in some cases")
            dos_method = st.radio("density of state calc", 
                  ["tetrahedra", "Gaussian"])
            if dos_method == "Gaussian":
                dos_broading = st.number_input("Gaissian broading in eV", 0.1)
                extgra_lines += 'dos_broading = "%f"  \n'%dos_broading

            occupations_type = st.radio("occupation type",
                    ["Fermi Dirac", "Fixed", "Cold Smearing", "MethfesselPaxton"])
            if occupations_type != "Fixed":
                cs, col1,col2 = st.columns([0.1,1,1])
                occ_smear = col1.number_input("occupation smear in eV", value =0.04)
                MP_order = col2.number_input("Order of Methefessel Paxton Occupation", value=2)


            md_tem_ctrl = st.radio("MD temperature control", 
                    ["Nose Hoover Chains","Anderson Rescaling"])
            md_integration_order = st.radio("MD Integration order",
                    ["2nd Velocity Verlet",
                     "3rd Beeman-Velocity Verlet",
                     "5th Beeman-Velocity Verlet"])
            md_number_of_nose_thermostats = st.number_input("Number of Nosethermostats", 5)
            ctrl_lines += 'relax_mass          ="' +relax_mass +'"  \n'
            ctrl_lines += 'dos_method          ="' +dos_method +'"  \n'
            ctrl_lines += 'occupations_type    ="' +occupations_type +'"  \n'
            ctrl_lines += 'occupation_electron_temperature_eV="%f"  \n'%occ_smear
            ctrl_lines += 'MP_order="%d"  \n'%MP_order
            ctrl_lines += 'poisson_solver      ="' +poisson_solver +'"  \n'
            ctrl_lines += 'md_temperature_control    ="' +md_tem_ctrl +'"  \n'
            ctrl_lines += 'md_integration_order="' +md_integration_order +'"  \n'
            ctrl_lines += 'md_number_of_nose_thermostats ="%d"  \n'%md_number_of_nose_thermostats

        ctrl_lines = "#******* CONTROL OPTIONS *******  \n"
        ctrl_lines += 'start_mode          ="' +start_mode +'"  \n'
        ctrl_lines += 'calculation_mode    ="' +calculation_mode +'"  \n'
        ctrl_lines += 'kohn_sham_solver    ="' +kohn_sham_solver +'"  \n'
        ctrl_lines += 'subdiag_driver      ="' +subdiag_driver +'"  \n'
        ctrl_lines += '#auto: if cuda available, use cusolver, otherwise use lapack for n<128 and scaplack for large system  \n'
        ctrl_lines += extra_lines
        ctrl_lines += '  \n'
        return ctrl_lines

def add_grid(cell):
    expand_ = st.expander("REAL SPACE GRID")
    with expand_:
        cs, col1, col2 = st.columns([0.1,2,2])
        cutoff = col1.number_input("equivalent cutoff energy(Ry) in plane wave", value=70.0, step=5.0,
                    help ="approximate equivalent cutoff energy in plane wave code, it may be different for different lattice types")
        grid_spacing_0 = 3.1415926/math.sqrt(cutoff)
        grid_spacing = col2.number_input("grid spacing(bohr)", value=grid_spacing_0,
                    help ="use grid spacing to determine the real space grid")
        if cell.unit == "angstrom" :
            grid_spacing = grid_spacing * 0.529177

            
        nx = int(round(cell.a/grid_spacing))
        ny = int(round(cell.b/grid_spacing))
        nz = int(round(cell.c/grid_spacing))
        i2 = 1
        for i in range(4):
            i2 *= 2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2
            h_max = max(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            h_min = min(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            anisotropy = h_max/h_min
            if(anisotropy > 1.1): break
        if i2 == 2:
            st.markdown("reduce grid spacing, anisotropy too large %f"%anisotropy) 
        else:
            i2 = i2//2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2
        grids_str = col1.text_input("number of grid Nx, Ny, Nz", value="%d %d %d"%(nx1, ny1, nz1))

        hx = cell.a/int(grids_str.split()[0])
        hy = cell.b/int(grids_str.split()[1])
        hz = cell.c/int(grids_str.split()[2])
        st.markdown("final grid spacing: hx =%f hy=%f hz=%f "%(hx,hy,hz) + cell.unit)
        anisotropy = max(hx,hy,hz)/min(hx,hy,hz)
        st.markdown("grid anisotropy =%f"%anisotropy)
        if(anisotropy >=1.1):
            st.markdown('<p style="color:red;">WARNGING: too big grid anisotropy, need to be <1.1 rmg wont run</p>', unsafe_allow_html=True)
        pot_grid= col2.number_input("rho pot grid refinement", value=2)
        grid_lines ='#******** REAL SPACE GRID ********   \n'
        grid_lines += 'wavefunction_grid="'+grids_str+'"  \n'
        grid_lines += 'potential_grid_refinement="%d"  \n'%pot_grid
        grid_lines += '  \n'

    return grid_lines
def add_scf():
    expand_ = st.expander("SCF & CONVERGENCE CONTROL")
    with expand_:
        cs, col1, col2, col3 = st.columns([0.1,1,1,1])
        max_scf_steps= col1.text_input("max scf steps", value="40")
        max_md_steps= col2.text_input("max md or relax steps", value="10")
        max_exx_steps= col3.text_input("max Exx steps for hybrid or HF", value="20")
        e_err= col1.text_input("energy convergence criterion", value="1.0e-9")
        rms_err= col2.text_input("rms convergence criterion", value="1.0e-7")
        precon_thres= col3.text_input("preconditioner threshold", value="0.0001")
        exx_convergence_criterion = col1.text_input("Exx convergence criterion for hybrid functional", value="1.0e-9") 
        vexx_fft_threshold = col2.text_input("Exx threshold for switch singlt to doulbe precision", value="1.0e-14")
        scf_lines = 'max_scf_steps = "'+max_scf_steps + '"  \n'
        scf_lines += 'max_md_steps = "'+max_md_steps + '"  \n'
        scf_lines += 'max_exx_steps = "'+max_md_steps + '"  \n'
        scf_lines += 'energy_convergence_criterion="' + e_err + '"  \n'
        scf_lines += 'rms_convergence_criterion = "' + rms_err +'"  \n'
        scf_lines += 'preconditioner_threshold = "' + precon_thres + '"  \n'
        scf_lines += 'exx_convergence_criterion = "' + exx_convergence_criterion + '"  \n'
        scf_lines += 'vexx_fft_threshold = "' + vexx_fft_threshold + '"  \n'


    return scf_lines


def add_mixing():
    expand_ = st.expander("MIXING OPTIONS")
    with expand_:
        charge_mixing_type = st.radio("charge density mixing type", 
                ["Broyden", "Pulay", "linear"])
        cs, col1, col2, col3 = st.columns([0.1,1,1,1])
        mix = col1.text_input("charge density mixing parameter",value="0.5") 
        mix_scale = col2.text_input("charge density mixing scale",value="0.5") 
        mix_order = col1.text_input("Broyden or Pulay order", value="5")
        refresh_step = col2.text_input("Broyden or Pulay refresh step", value="100")
        mixing_lines  = 'charge_mixing_type = "'+ charge_mixing_type +'"  \n'
        mixing_lines += 'charge_density_mixing ="' + mix +'"  \n'
        if charge_mixing_type == "Broyden":
            mixing_lines += 'charge_broyden_order = "' + mix_order + '"  \n'
            mixing_lines += 'charge_broyden_scale = "' +mix_scale + '"  \n'
            mixing_lines += 'charge_broyden_refresh = "' +refresh_step + '"  \n'
        elif charge_mixing_type == "Pulay":
            mixing_lines += 'charge_pulay_order = "' + mix_order + '"  \n'
            mixing_lines += 'charge_pulay_scale = "' +mix_scale + '"  \n'
            mixing_lines += 'charge_pulay_refresh = "' +refresh_step + '"  \n'
            pulay_gspace = col1.checkbox("Pulay mixing in G space", False)
            drho_precond = col2.checkbox("scale q^2/(q^2+q0^2)", False)
            drho_precond_q0= col1.number_input("q0 value", 0.5)
            mixing_lines += 'charge_pulay_Gspace = "' + str(pulay_gspace)+ '"  \n'
            mixing_lines += 'drho_precond = "' +str(drho_precond) + '"  \n'
            mixing_lines += 'drho_precond_q0="%f"  \n'%drho_precond_q0

    return mixing_lines

def add_xc(species):
    expand_ = st.expander("EXCHANGE CORRELATION POTENTIAL")
    xc_lines = '#*****Exchange Correlation ******  \n'
    with expand_:
        xc_type = st.radio("exchange correlation type", 
                ["AUTO_XC", "LDA", "GGA XB CP", "PW91", "GGA BLYP", "GGA PBE",
                  "REVPBE", "PW86PBE", "PBESOL", "PBE0", "HSE", "B3LYP", "gaupbe", 
                  "vdw-df", "VDW-DF", "hartree-fock"], 
                help = "AUTO_XC: XC will be determined from pseudopotential")
        xc_lines += 'exchange_correlation_type="'+xc_type +'"  \n'
        xc_lines += '#AUTO_XC: XC will be determined from pseudopotential  \n'
        vdw_corr = st.radio("empirical van der Waals correction", 
                ["None", "DFT-D2", "Grimme-D2","DFT-D3"])
        cs, col1, col2 = st.columns([0.1,1,1])
        if xc_type in ["PBESOL", "PBE0", "HSE", "B3LYP", "gaupbe"]:
            exx_mode = col1.radio("Exx mode", ["Local fft", "Distributed fft"])

            exxdiv_treatment = col2.radio("Exx divergence treatment", 
                    ["gygi-baldereschi", "none"])
            x_gamma_extrapolation = col1.checkbox("x_gamma_extrapolation", True)
            exx_fracton = col2.text_input("the fraction of Exx for hybrid functional", value="-1.0", 
                help="negative value: the fraction determined by code for different hybrid functionals")
            xc_lines += 'exx_mode = "' + exx_mode + '"  \n'
            xc_lines += 'exxdiv_treatment = "' + exxdiv_treatment +'"  \n'
            xc_lines += 'x_gamma_extrapolation ="' + str(x_gamma_extrapolation) +'"  \n'
            xc_lines += 'exx_fracton = "' + exx_fracton +'"  \n'
        if xc_type in ["vdw-df", "VDW-DF"]:
            vdwdf_grid = col2.radio("grid for vdw corr",
                    ["Coarse", "Fine"])
            vdwdf_kernel_filepath = col1.text_input("van der Waals Kernel file", "vdW_kernel_table")

            xc_lines += 'vdwdf_grid_type     ="' +vdwdf_grid +'"  \n'
            xc_lines += 'vdwdf_kernel_filepath ="%s"  \n'%vdwdf_kernel_filepath 

        if vdw_corr != "None":
            xc_lines += 'vdw_corr            ="' +vdw_corr +'"  \n'

        ldaU_mode = st.radio("LDA+U type",["None","Simple"])
        if(ldaU_mode == "Simple"):
            xc_lines += 'ldaU_mode = "%s"  \n'%ldaU_mode
            cs, col1, col2 = st.columns([0.1,1,1])
            Hubbard_U = col1.text_area("HUbbard U for species", "", 
                help= "Ni 6.5 3d 0.0 0.0 0.0 for each specie ")
            xc_lines += 'Hubbard_U ="  \n' + Hubbard_U + '  \n"  \n'
            ldau_mixing_type = col1.radio("mixing type for ldau occupations", ["Linear", "Pulay"])
            xc_lines += 'ldau_mixing_type = "%s"  \n'%ldau_mixing_type

            ldaU_radius = col2.number_input("orbital radius for LDA+U projection", 9.0)
            xc_lines += 'ldaU_radius = "%f"  \n'%ldaU_radius

            ldau_mixing = col1.number_input("mixing fractions", 1.0)
            xc_lines += 'ldau_mixing = "%f"  \n'%ldau_mixing
            if ldau_mixing_type == "Pulay":
                ldau_pulay_order = col2.number_input("Pulay order for lda+u mixing", 5)
                ldau_pulay_scale = col1.number_input("Pulay scale for lda+u mixing", 0.8)
                ldau_pulay_refresh = col2.number_input("Pulay refresh steps for lda+u mixing", 100)
                xc_lines += 'ldau_pulay_order = "%d"  \n'%ldau_pulay_order
                xc_lines += 'ldau_pulay_scale = "%f"  \n'% ldau_pulay_scale
                xc_lines += 'ldau_pulay_refresh = "%d"  \n'%ldau_pulay_refresh

    xc_lines += '  \n'
    return xc_lines

def add_qmcpack():
    expand_ = st.expander("QMCPACK INTERFACE")
    with expand_:
        cs, col1, col2 = st.columns([0.1,1,1])
        qmcpack = col1.checkbox("Write out file for QMCPACK")
        qmcpack_lines = 'write_qmcpack_restart = "' + str(qmcpack) + '"  \n'
        cs, col1, col2 = st.columns([0.1,1,1])
        if qmcpack:
            exx_integrals_filepath = col1.text_input("file name for afqmc", value="afqmc_rmg")
            ExxIntCholosky = col1.checkbox("Cholesky factorization for Vexx", True)
            ExxCholMax = col2.text_input("maximum Cholesky vectors", value="8")
            exx_int_flag = col2.checkbox("Calculate Exack exchange integrals", True)
            qmc_nband = col1.text_input("number of bands for qmcpack", value="0", 
                    help="default value 0: use the number of states")

            qmcpack_lines +='exx_integrals_filepath = "' + exx_integrals_filepath +'"  \n'
            qmcpack_lines +='ExxIntCholosky = "' + str(ExxIntCholosky) +'"  \n'
            qmcpack_lines +='ExxCholMax = "' +  ExxCholMax + '"  \n'
            qmcpack_lines +='exx_int_flag = "' +  str(exx_int_flag) +'"  \n'
            qmcpack_lines +='qmc_nband = "' + qmc_nband +'"  \n'
    return qmcpack_lines


def add_lattice(bounding_box):
    expand_ = st.expander("LATTICE INFO in unit of Anstrom")
    lattvec = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    #estimate the a, b, c = bounding box + 5 Angstrom
    a = bounding_box[1] - bounding_box[0] +5.0
    b = bounding_box[3] - bounding_box[2] +5.0
    c = bounding_box[5] - bounding_box[4] +5.0
    lattvec = [[a,0.0,0.0],[0.0,b,0.0],[0.0,0.0,c]]

    ibrav = 0
    with expand_:
        st.markdown("min_x = %f, max_x = %f"%(bounding_box[0], bounding_box[1]))
        st.markdown("min_y = %f, max_x = %f"%(bounding_box[2], bounding_box[3]))
        st.markdown("min_z = %f, max_z = %f"%(bounding_box[4], bounding_box[5]))
        cs, col1 = st.columns([0.1,1])
        ibrav_str = st.radio("Bravais lattice type", 
                ["Orthorhombic", "Simple Cubic", "FCC", "BCC", "Hexagonal", "do not know"],
                help = "choose do not know for others")
        cs, col1,col2, col3 = st.columns([0.1,1,1,1])
        if ibrav_str == "do not know":
            ibrav = 0
            lattvec_str = col1.text_area("lattice vector in Angstrom", 
                    help = " must be 3x3 numbers")
            mat = lattvec_str.split("\n")
            if len(mat) == 3:
                for i in range(3):
                    vec = mat[i].split()
                    for j in range(3):
                        lattvec[i][j] = float(vec[j])
            a = sqrt(lattvec[0][0] *lattvec[0][0] +lattvec[0][1] *lattvec[0][1] +lattvec[0][2] *lattvec[0][2] )
            b = sqrt(lattvec[1][0] *lattvec[1][0] +lattvec[1][1] *lattvec[1][1] +lattvec[1][2] *lattvec[1][2] )
            c = sqrt(lattvec[2][0] *lattvec[2][0] +lattvec[2][1] *lattvec[2][1] +lattvec[2][2] *lattvec[2][2] )
        elif ibrav_str == "Simple Cubic":
            ibrav = 1
            a = col1.number_input("length a", value=a)
            b = a
            c = a
            lattvec = [[a,0.0,0.0],[0.0,b,0.0],[0.0,0.0,c]]
        elif ibrav_str == "FCC":
            ibrav = 2
            a = col1.number_input("length a", value=a)
            b = a
            c = a
            lattvec = [[0.5 * a,0.5*a,0.0],[0.5*a,0.0,0.5*a],[0.0,0.5*a,0.5*a]]
        elif ibrav_str =="BCC":       
            ibrav = 3
            a = col1.number_input("length a", value=a)
            b = a
            c = a
            lattvec = [[-0.5 * a,0.5*a, 0.5*a],[0.5*a,-0.5*a,0.5*a],[0.5*a,0.5*a,-0.5*a]]
        elif ibrav_str == "Orthorhombic": 
            ibrav = 8
            a = col1.number_input("length a", value=a)
            b = col2.number_input("length b", value=b)
            c = col3.number_input("length c", value=c)
            lattvec = [[a,0.0,0.0],[0.0,b,0.0],[0.0,0.0,c]]
        elif ibrav_str == "Hexagonal":
            ibrav = 4
            a = col1.number_input("length a", value=a)
            b = a
            c = col3.number_input("length c", value=c)
            lattvec = [[a,0.0,0.0],[-0.5*a,sqrt(3.0)/2 * a,0.0],[0.0,0.0,c]]
    return (ibrav, a,b,c, lattvec)            
def add_IOctrl():
    expand_ = st.expander("IO: files and paths")
    with expand_:
        cs, col1, col2 = st.columns([0.2,1,1])
        verbose = col1.checkbox("print out more in log file if True", False)
        cs, col1, col2 = st.columns([0.2,1,1])
        input_wave_function_file = col1.text_input("input wave function file", "Waves/wave.out")
        output_wave_function_file = col2.text_input("output wave function file", "Waves/wave.out")
        write_serial_restart = col1.checkbox("write a serial file for restart", False)
        read_serial_restart = col2.checkbox("restart from a serial file", False)
        compressed_infile = col1.checkbox("read the compressed file for restart", True)
        compressed_outfile = col2.checkbox("compress the out wave file", True)
        input_tddft_file = col1.text_input("input tddft file", "Waves/wave_tddft.out")
        output_tddft_file = col2.text_input("output TDDFT file", "Waves/wave_tddft.out")
        nvme_weights = col1.checkbox("map nonlocal projectors to disk", False)
        nvme_work = col2.checkbox("map work arrays to disk", False)
        nvme_orbitals = col1.checkbox("map orbitals to disk", False)
        nvme_qfunctons = col2.checkbox("map qfunctions to disk", False)
        nvme_weights_filepath = col1.text_input("nvme directory for non-local projectors", "Weights/")
        nvme_work_filepath = col2.text_input("nvme directory for work arrays", "Work/")
        nvme_orbitals_filepath = col1.text_input(" nvme directory for orbitals", "Orbitals/")
        qfunction_filepath = col2.text_input("nvme directory for Qfunction", "Qfunctions/")
        cube_rho = col1.checkbox("output rho in cube format", True)
        output_rho_xsf = col2.checkbox("output rho in xsf format", False)
        cube_vh = col1.checkbox("output vh in cube format",  False)
        cube_pot = col2.checkbox("output pot in cube format", False)
        write_data_period = col1.number_input("steps to write the restart file", 5)
        write_eigvals_period = col2.number_input("steps to write eigenvalues",5)    
        cube_states_list = col1.text_input("list of states to plot in cube format", "", 
                help="0,1-3,6,9 will print states 0, 1 to 3, 6 and 9")

    IO_lines = ""
    IO_lines += 'verbose = "%s"  \n'%str(verbose)
    IO_lines += 'input_wave_function_file = "%s"  \n'%input_wave_function_file  
    IO_lines += 'output_wave_function_file = "%s"  \n'%output_wave_function_file 
    IO_lines += 'write_serial_restart = "%s"  \n'%str(write_serial_restart) 
    IO_lines += 'read_serial_restart = "%s"  \n'%str(read_serial_restart) 
    IO_lines += 'compressed_infile = "%s"  \n'%str(compressed_infile) 
    IO_lines += 'compressed_outfile = "%s"  \n'%str(compressed_outfile) 
    IO_lines += 'input_tddft_file = "%s"  \n'%input_tddft_file
    IO_lines += 'output_tddft_file = "%s"  \n'%output_tddft_file
    IO_lines += 'nvme_weights = "%s"  \n'%str(nvme_weights) 
    IO_lines += 'nvme_work = "%s"  \n'%str(nvme_work) 
    IO_lines += 'nvme_orbitals = "%s"  \n'%str(nvme_orbitals) 
    IO_lines += 'nvme_qfunctons = "%s"  \n'%str(nvme_qfunctons) 
    IO_lines += 'nvme_weights_filepath = "%s"  \n'%nvme_weights_filepath
    IO_lines += 'nvme_work_filepath = "%s"  \n'%nvme_work_filepath
    IO_lines += 'nvme_orbitals_filepath = "%s"  \n'%nvme_orbitals_filepath
    IO_lines += 'qfunction_filepath = "%s"  \n'%qfunction_filepath
    IO_lines += 'cube_rho = "%s"  \n'%str(cube_rho) 
    IO_lines += 'output_rho_xsf = "%s"  \n'%str(output_rho_xsf) 
    IO_lines += 'cube_vh = "%s"  \n'%str(cube_vh) 
    IO_lines += 'cube_pot = "%s"  \n'%str(cube_pot) 
    if cube_states_list != "":
        IO_lines += 'cube_states_list = "%s"  \n'%cube_states_list
    IO_lines += 'write_data_period = "%d"  \n'%write_data_period 
    IO_lines += 'write_eigvals_period = "%d"  \n'%write_eigvals_period
    return IO_lines
def add_spin(species, atoms):
    expand_ = st.expander("SPIN and MAGNETIZATION")
    dict_mag_species = {}
    dict_mag_species = {}
    angle1_species = {}
    angle2_species = {}

    mag = []
    for atom in atoms:
        mag.append([0.0, 0.0, 0.0])
    spin_lines = ""
    with expand_:
        nspin_str = st.radio("spin setup", ["None", "spin polarization", "spin orbit coupling"])

        if(nspin_str == "spin polarization"):
            spin_lines += 'spin_polarization = "True"  \n'
            AFM = st.checkbox("Anti-Ferromagnetic?", False)
            spin_lines += 'AFM = "%s"  \n'%str(AFM) 

            s_or_a = st.radio("Init Magnetization", ["by species", "by atoms"])
            if s_or_a == "by species":
                for sp in species:
                    dict_mag_species[sp] = st.number_input(sp, min_value = -1.0, max_value = 1.0, value =0.0, 
                            help="spin up and down density: (0.5+x, 0.5-x) of total atomic charge density")
                for i in range(len(atoms)):
                    mag[i][0] = dict_mag_species[atoms[i][0]]
            else:
                for i in range(len(atoms)):
                    tem_str = "atom " + str(i) +": up down spin difference"
                    mag[i][0] = st.number_input(tem_str, 0.0)
        
        if(nspin_str == "spin orbit coupling"):
            spin_lines += 'spinorbit = "True"  \n'
            spin_lines += 'noncollinear = "True"  \n'
            s_or_a = st.radio("Init Magnetization", ["by species", "by atoms"])
            cs, col1, col2, col3 = st.columns([0.2,1,1,1])
            if s_or_a == "by species":
                for sp in species:
                    dict_mag_species[sp] = col1.number_input(sp + " mag", 0.0, 
                            help="spin up and down density: (0.5+x, 0.5-x) of total atomic charge density")
                    angle1_species[sp] = col2.number_input(sp + " angle1", 0, help = "180 indicate the -z direction") 
                    angle2_species[sp] = col3.number_input(sp + " angle2", 0, help = "directiopn in xy plane") 
                for i in range(len(atoms)):
                    mag[i][0] = dict_mag_species[atoms[i][0]]
                    mag[i][1] = angle1_species[atoms[i][0]]
                    mag[i][2] = angle2_species[atoms[i][0]]
            else:
                for i in range(len(atoms)):
                    tem_str = "atom " + str(i) +": "
                    mag[i][0] = col1.number_input(tem_str + " mag", 0.0)
                    mag[i][1] = col2.number_input(tem_str + "angle1", 0)
                    mag[i][2] = col3.number_input(tem_str + "angle2", 0)


    return spin_lines, mag


def add_misc():
    expand_ = st.expander("MISC OPTIONS")
    misc_lines = ""
    with expand_:
        col0, col1, col2 = st.columns([1,1,1])

        dftd3_version= col0.number_input("dftd3_version", 3)
        misc_lines += 'dftd3_version = "%d"  \n'%dftd3_version
        charge_analysis_period= col1.number_input("charge_analysis_period", 0)
        misc_lines += 'charge_analysis_period = "%d"  \n'%charge_analysis_period
        omp_threads_per_node= col2.number_input("omp_threads_per_node", 0)
        misc_lines += 'omp_threads_per_node = "%d"  \n'%omp_threads_per_node
        fd_allocation_limit= col0.number_input("fd_allocation_limit", 65536)
        misc_lines += 'fd_allocation_limit = "%d"  \n'%fd_allocation_limit
        rmg_threads_per_node= col1.number_input("rmg_threads_per_node", 0)
        misc_lines += 'rmg_threads_per_node = "%d"  \n'%rmg_threads_per_node
        unoccupied_states_per_kpoint= col2.number_input("unoccupied_states_per_kpoint", 10)
        misc_lines += 'unoccupied_states_per_kpoint = "%d"  \n'%unoccupied_states_per_kpoint
        state_block_size= col0.number_input("state_block_size", 64)
        misc_lines += 'state_block_size = "%d"  \n'%state_block_size
        extra_random_lcao_states= col1.number_input("extra_random_lcao_states", 0)
        misc_lines += 'extra_random_lcao_states = "%d"  \n'%extra_random_lcao_states
        kohn_sham_fd_order= col2.number_input("kohn_sham_fd_order", 8)
        misc_lines += 'kohn_sham_fd_order = "%d"  \n'%kohn_sham_fd_order
        force_grad_order= col0.number_input("force_grad_order", 8)
        misc_lines += 'force_grad_order = "%d"  \n'%force_grad_order
        non_local_block_size= col1.number_input("non_local_block_size", 512)
        misc_lines += 'non_local_block_size = "%d"  \n'%non_local_block_size
        dynamic_time_delay= col2.number_input("dynamic_time_delay", 5)
        misc_lines += 'dynamic_time_delay = "%d"  \n'%dynamic_time_delay
        dynamic_time_counter= col0.number_input("dynamic_time_counter", 0)
        misc_lines += 'dynamic_time_counter = "%d"  \n'%dynamic_time_counter
        scf_steps_offset= col1.number_input("scf_steps_offset", 0)
        misc_lines += 'scf_steps_offset = "%d"  \n'%scf_steps_offset
        total_scf_steps_offset= col2.number_input("total_scf_steps_offset", 0)
        misc_lines += 'total_scf_steps_offset = "%d"  \n'%total_scf_steps_offset
        md_steps_offset= col0.number_input("md_steps_offset", 0)
        misc_lines += 'md_steps_offset = "%d"  \n'%md_steps_offset
        coalesce_factor= col1.number_input("coalesce_factor", 4)
        misc_lines += 'coalesce_factor = "%d"  \n'%coalesce_factor
        folded_spectrum_iterations= col2.number_input("folded_spectrum_iterations", 2)
        misc_lines += 'folded_spectrum_iterations = "%d"  \n'%folded_spectrum_iterations
        vxc_diag_nmin= col0.number_input("vxc_diag_nmin", 1)
        misc_lines += 'vxc_diag_nmin = "%d"  \n'%vxc_diag_nmin
        vxc_diag_nmax= col1.number_input("vxc_diag_nmax", 1)
        misc_lines += 'vxc_diag_nmax = "%d"  \n'%vxc_diag_nmax
        num_wanniers= col2.number_input("num_wanniers", 0)
        misc_lines += 'num_wanniers = "%d"  \n'%num_wanniers
        wannier90_scdm= col0.number_input("wannier90_scdm", 0)
        misc_lines += 'wannier90_scdm = "%d"  \n'%wannier90_scdm
        md_temperature= col1.number_input("md_temperature", 300.0)
        misc_lines += 'md_temperature = "%f"  \n'%md_temperature
        md_nose_oscillation_frequency_THz= col2.number_input("md_nose_oscillation_frequency_THz", 15.59)
        misc_lines += 'md_nose_oscillation_frequency_THz = "%f"  \n'%md_nose_oscillation_frequency_THz
        filter_factor= col0.number_input("filter_factor", 0.25)
        misc_lines += 'filter_factor = "%f"  \n'%filter_factor
        potential_acceleration_constant_step= col1.number_input("potential_acceleration_constant_step", 0.0)
        misc_lines += 'potential_acceleration_constant_step = "%f"  \n'%potential_acceleration_constant_step
        ionic_time_step= col2.number_input("ionic_time_step", 50.0)
        misc_lines += 'ionic_time_step = "%f"  \n'%ionic_time_step
        ionic_time_step_increase= col0.number_input("ionic_time_step_increase", 1.1)
        misc_lines += 'ionic_time_step_increase = "%f"  \n'%ionic_time_step_increase
        ionic_time_step_decrease= col1.number_input("ionic_time_step_decrease", 0.5)
        misc_lines += 'ionic_time_step_decrease = "%f"  \n'%ionic_time_step_decrease
        max_ionic_time_step= col2.number_input("max_ionic_time_step", 150.0)
        misc_lines += 'max_ionic_time_step = "%f"  \n'%max_ionic_time_step
        system_charge= col0.number_input("system_charge", 0.0)
        misc_lines += 'system_charge = "%f"  \n'%system_charge
        unoccupied_tol_factor= col1.number_input("unoccupied_tol_factor", 1000.0)
        misc_lines += 'unoccupied_tol_factor = "%f"  \n'%unoccupied_tol_factor
        projector_expansion_factor= col2.number_input("projector_expansion_factor", 1.0)
        misc_lines += 'projector_expansion_factor = "%f"  \n'%projector_expansion_factor
        folded_spectrum_width= col0.number_input("folded_spectrum_width", 0.3)
        misc_lines += 'folded_spectrum_width = "%f"  \n'%folded_spectrum_width
        ecutrho= col1.number_input("ecutrho", 0.0)
        misc_lines += 'ecutrho = "%f"  \n'%ecutrho
        ecutwfc= col2.number_input("ecutwfc", 0.0)
        misc_lines += 'ecutwfc = "%f"  \n'%ecutwfc
        test_energy= col0.number_input("test_energy", 0.0)
        misc_lines += 'test_energy = "%f"  \n'%test_energy
        test_energy_tolerance= col1.number_input("test_energy_tolerance", 1.0e-7)
        misc_lines += 'test_energy_tolerance = "%f"  \n'%test_energy_tolerance
        test_bond_length= col2.number_input("test_bond_length", 0.0)
        misc_lines += 'test_bond_length = "%f"  \n'%test_bond_length
        test_bond_length_tolerance= col0.number_input("test_bond_length_tolerance", 1.0e-3)
        misc_lines += 'test_bond_length_tolerance = "%f"  \n'%test_bond_length_tolerance
        relax_max_force= col1.number_input("relax_max_force", 2.5E-3)
        misc_lines += 'relax_max_force = "%f"  \n'%relax_max_force
        stress_convergence_criterion= col2.number_input("stress_convergence_criterion", 0.5)
        misc_lines += 'stress_convergence_criterion = "%f"  \n'%stress_convergence_criterion
        gw_residual_convergence_criterion= col0.number_input("gw_residual_convergence_criterion", 1.0e-6)
        misc_lines += 'gw_residual_convergence_criterion = "%f"  \n'%gw_residual_convergence_criterion
        gw_residual_fraction= col1.number_input("gw_residual_fraction", 0.90)
        misc_lines += 'gw_residual_fraction = "%f"  \n'%gw_residual_fraction
        hartree_rms_ratio= col2.number_input("hartree_rms_ratio", 100000.0)
        misc_lines += 'hartree_rms_ratio = "%f"  \n'%hartree_rms_ratio
        electric_field_magnitude= col0.number_input("electric_field_magnitude", 0.0)
        misc_lines += 'electric_field_magnitude = "%f"  \n'%electric_field_magnitude
        wannier90_scdm_mu= col1.number_input("wannier90_scdm_mu", 0.0)
        misc_lines += 'wannier90_scdm_mu = "%f"  \n'%wannier90_scdm_mu
        wannier90_scdm_sigma= col2.number_input("wannier90_scdm_sigma", 1.0)
        misc_lines += 'wannier90_scdm_sigma = "%f"  \n'%wannier90_scdm_sigma
        stress= col0.checkbox("stress",False)
        misc_lines += 'stress = "%s"  \n'%str(stress)
        cell_relax= col1.checkbox("cell_relax",False)
        misc_lines += 'cell_relax = "%s"  \n'%str(cell_relax)
        dipole_moment= col2.checkbox("dipole_moment",False)
        misc_lines += 'dipole_moment = "%s"  \n'%str(dipole_moment)
        use_gpu_fd= col0.checkbox("use_gpu_fd",False)
        misc_lines += 'use_gpu_fd = "%s"  \n'%str(use_gpu_fd)
        laplacian_offdiag= col1.checkbox("laplacian_offdiag",False)
        misc_lines += 'laplacian_offdiag = "%s"  \n'%str(laplacian_offdiag)
        laplacian_autocoeff= col2.checkbox("laplacian_autocoeff",False)
        misc_lines += 'laplacian_autocoeff = "%s"  \n'%str(laplacian_autocoeff)
        use_cpdgemr2d= col0.checkbox("use_cpdgemr2d",True)
        misc_lines += 'use_cpdgemr2d = "%s"  \n'%str(use_cpdgemr2d)
        use_symmetry= col1.checkbox("use_symmetry",True)
        misc_lines += 'use_symmetry = "%d"  \n'%(use_symmetry)
        frac_symmetry= col2.checkbox("frac_symmetry",True)
        misc_lines += 'frac_symmetry = "%s"  \n'%str(frac_symmetry)
        rmg2bgw= col0.checkbox("rmg2bgw",False)
        misc_lines += 'rmg2bgw = "%s"  \n'%str(rmg2bgw)
        pin_nonlocal_weights= col1.checkbox("pin_nonlocal_weights",False)
        misc_lines += 'pin_nonlocal_weights = "%s"  \n'%str(pin_nonlocal_weights)
        use_cublasxt= col2.checkbox("use_cublasxt",False)
        misc_lines += 'use_cublasxt = "%s"  \n'%str(use_cublasxt)
        use_bessel_projectors= col0.checkbox("use_bessel_projectors",False)
        misc_lines += 'use_bessel_projectors = "%s"  \n'%str(use_bessel_projectors)
        write_orbital_overlaps= col1.checkbox("write_orbital_overlaps",False)
        misc_lines += 'write_orbital_overlaps = "%s"  \n'%str(write_orbital_overlaps)
        kohn_sham_ke_fft= col2.checkbox("kohn_sham_ke_fft",False)
        misc_lines += 'kohn_sham_ke_fft = "%s"  \n'%str(kohn_sham_ke_fft)
        fast_density= col0.checkbox("fast_density",True)
        misc_lines += 'fast_density = "%s"  \n'%str(fast_density)
        lcao_use_empty_orbitals= col1.checkbox("lcao_use_empty_orbitals",False)
        misc_lines += 'lcao_use_empty_orbitals = "%s"  \n'%str(lcao_use_empty_orbitals)
        write_qmcpack_restart_localized= col2.checkbox("write_qmcpack_restart_localized",False)
        misc_lines += 'write_qmcpack_restart_localized = "%s"  \n'%str(write_qmcpack_restart_localized)
        alt_laplacian= col0.checkbox("alt_laplacian",True)
        misc_lines += 'alt_laplacian = "%s"  \n'%str(alt_laplacian)
        use_alt_zgemm= col1.checkbox("use_alt_zgemm",False)
        misc_lines += 'use_alt_zgemm = "%s"  \n'%str(use_alt_zgemm)
        filter_dpot= col2.checkbox("filter_dpot",False)
        misc_lines += 'filter_dpot = "%s"  \n'%str(filter_dpot)
        sqrt_interpolation= col0.checkbox("sqrt_interpolation",False)
        misc_lines += 'sqrt_interpolation = "%s"  \n'%str(sqrt_interpolation)
        renormalize_forces= col1.checkbox("renormalize_forces",True)
        misc_lines += 'renormalize_forces = "%s"  \n'%str(renormalize_forces)
        coalesce_states= col2.checkbox("coalesce_states",False)
        misc_lines += 'coalesce_states = "%s"  \n'%str(coalesce_states)
        equal_initial_density= col0.checkbox("equal_initial_density",False)
        misc_lines += 'equal_initial_density = "%s"  \n'%str(equal_initial_density)
        write_pdos= col1.checkbox("write_pdos",False)
        misc_lines += 'write_pdos = "%s"  \n'%str(write_pdos)
        folded_spectrum= col2.checkbox("folded_spectrum",False)
        misc_lines += 'folded_spectrum = "%s"  \n'%str(folded_spectrum)
        use_numa= col0.checkbox("use_numa",True)
        misc_lines += 'use_numa = "%s"  \n'%str(use_numa)
        use_hwloc= col1.checkbox("use_hwloc",False)
        misc_lines += 'use_hwloc = "%s"  \n'%str(use_hwloc)
        use_async_allreduce= col2.checkbox("use_async_allreduce",True)
        misc_lines += 'use_async_allreduce = "%s"  \n'%str(use_async_allreduce)
        mpi_queue_mode= col0.checkbox("mpi_queue_mode",True)
        misc_lines += 'mpi_queue_mode = "%s"  \n'%str(mpi_queue_mode)
        spin_manager_thread= col1.checkbox("spin_manager_thread",True)
        misc_lines += 'spin_manager_thread = "%s"  \n'%str(spin_manager_thread)
        spin_worker_threads= col2.checkbox("spin_worker_threads",True)
        misc_lines += 'spin_worker_threads = "%s"  \n'%str(spin_worker_threads)
        require_huge_pages= col0.checkbox("require_huge_pages",False)
        misc_lines += 'require_huge_pages = "%s"  \n'%str(require_huge_pages)
        relax_dynamic_timestep= col1.checkbox("relax_dynamic_timestep",False)
        misc_lines += 'relax_dynamic_timestep = "%s"  \n'%str(relax_dynamic_timestep)
        freeze_occupied= col2.checkbox("freeze_occupied",False)
        misc_lines += 'freeze_occupied = "%s"  \n'%str(freeze_occupied)
        md_randomize_velocity= col0.checkbox("md_randomize_velocity",True)
        misc_lines += 'md_randomize_velocity = "%s"  \n'%str(md_randomize_velocity)
        time_reversal= col1.checkbox("time_reversal",True)
        misc_lines += 'time_reversal = "%s"  \n'%str(time_reversal)
        wannier90= col2.checkbox("wannier90",False)
        misc_lines += 'wannier90 = "%s"  \n'%str(wannier90)
        processor_grid= st.text_input("processor_grid", "1 1 1")
        misc_lines += 'processor_grid = "%s"  \n'%processor_grid
        dipole_correction= st.text_input("dipole_correction", "0  0  0")
        misc_lines += 'dipole_correction = "%s"  \n'%dipole_correction
        cell_movable= st.text_input("cell_movable", "0 0 0 0 0 0 0 0 0")
        misc_lines += 'cell_movable = "%s"  \n'%cell_movable
        atomic_orbital_type= st.radio("atomic_orbital_type", ["delocalized","localized"])
        misc_lines += 'atomic_orbital_type = "%s"  \n'%atomic_orbital_type
        electric_field_vector= st.text_input("electric_field_vector", "0  0  1")
        misc_lines += 'electric_field_vector = "%s"  \n'%electric_field_vector
        states_count_and_occupation_spin_up= st.text_input("states_count_and_occupation_spin_up", "")
        if states_count_and_occupation_spin_up != "":
            misc_lines += 'states_count_and_occupation_spin_up = "%s"  \n'%states_count_and_occupation_spin_up
        states_count_and_occupation_spin_down= st.text_input("states_count_and_occupation_spin_down", "")
        if states_count_and_occupation_spin_down != "":
            misc_lines += 'states_count_and_occupation_spin_down = "%s"  \n'%states_count_and_occupation_spin_down
        states_count_and_occupation= st.text_input("states_count_and_occupation", "")
        if states_count_and_occupation != "":
            misc_lines += 'states_count_and_occupation = "%s"  \n'%states_count_and_occupation
        energy_output_units= st.radio("energy_output_units", ["Hartrees", "Rydbergs"])
        misc_lines += 'energy_output_units = "%s"  \n'%energy_output_units
        interpolation_type= st.radio("interpolation_type", ["FFT", "Cubic Polynomial", "prolong"])
        misc_lines += 'interpolation_type = "%s"  \n'%interpolation_type
    return misc_lines     
def add_orbital_info(species):  
    expand_ = st.expander("Localized Orbital information")
    orbital_dict={}
    with expand_:
        st.markdown("number of orbitals per atom and their radius")
        cstart, col1, col2 = st.columns([0.2,1,1])
        for sp in species:
            num_orb = 4
            if sp in num_orbitals_dict:
                num_orb = num_orbitals_dict[sp]
            num_orb = col1.number_input("number of orbital for %s:"%sp, num_orb)
            radius = col2.number_input("radius (bohr) for %s:"%sp, 6.5)
            orbital_dict[sp] = [num_orb, radius]
    return orbital_dict        

def add_lead_info():  
    with expand_:
        eq_left_right = st.checkbox("left lead = right lead?", True)
        a_lead1 = st.number_input("length of left lead (lead1)", 10.0)
        a_lead2 = a_lead1
        if eq_left_right:
            a_lead2 = a_lead1
        else:
            a_lead2 = st.number_input("length of right lead (lead2)", 10.0)
    return a_lead1, a_lead2, eq_left_right

def add_grid_negf(cell):
    expand_ = st.expander("REAL SPACE GRID for NEGF")
    with expand_:
        cs, col1, col2 = st.columns([0.1,2,2])
        cutoff = col1.number_input("equivalent cutoff energy(Ry) in plane wave", value=70.0, step=5.0,
                    help ="approximate equivalent cutoff energy in plane wave code, it may be different for different lattice types")
        grid_spacing_0 = 3.1415926/math.sqrt(cutoff)
        grid_spacing = col2.number_input("grid spacing(bohr)", value=grid_spacing_0,
                    help ="use grid spacing to determine the real space grid")
        if cell.unit == "angstrom" :
            grid_spacing = grid_spacing * 0.529177

            
        nx = int(round(cell.a/grid_spacing))
        ny = int(round(cell.b/grid_spacing))
        nz = int(round(cell.c/grid_spacing))
        i2 = 1
        for i in range(4):
            i2 *= 2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2
            h_max = max(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            h_min = min(cell.a/nx1, cell.b/ny1, cell.c/nz1)
            anisotropy = h_max/h_min
            if(anisotropy > 1.1): break
        if i2 == 2:
            st.markdown("reduce grid spacing, anisotropy too large %f"%anisotropy) 
        else:
            i2 = i2//2
            nx1 = (nx+i2-1)//i2 * i2
            ny1 = (ny+i2-1)//i2 * i2
            nz1 = (nz+i2-1)//i2 * i2

        hx = cell.a/int(nx1)
        hy = cell.b/int(ny1)
        hz = cell.c/int(nz1)
            
        eq_left_right = st.checkbox("left lead = right lead?", True)
        a_lead1 = st.number_input("length of left lead (lead1)", 8.156)
        a_lead2 = a_lead1
        if eq_left_right:
            a_lead2 = a_lead1
        else:
            a_lead2 = st.number_input("length of right lead (lead2)", 8.156)
        a_center = cell.a - a_lead1 - a_lead2

        nx_lead1 = int(round(a_lead1/hx))
        nx_lead2 = int(round(a_lead2/hx))
        nx_center = int(round( a_center/hx))
        nx_lead1 = (nx_lead1+3)//4 * 4
        nx_lead2 = (nx_lead2+3)//4 * 4
        nx_center = (nx_center+3)//4 * 4

        hx_lead1 = a_lead1/int(nx_lead1)
        hx_lead2 = a_lead2/int(nx_lead2)
        hx_center = a_center/int(nx_center)
        hy = cell.b/int(ny1)
        hz = cell.c/int(nz1)
        st.markdown("final grid spacing: hx =%f %f %f hy=%f hz=%f "%(hx_lead1, hx_lead2, hx_center,hy,hz) + cell.unit)
        anisotropy = max(hx_lead1, hx_lead2, hx_center,hy,hz)/min(hx_lead1, hx_lead2, hx_center,hy,hz)
        st.markdown("grid anisotropy =%f"%anisotropy)
        if(anisotropy >=1.1):
            st.markdown('<p style="color:red;">WARNGING: too big grid anisotropy, need to be <1.1 rmg wont run</p>', unsafe_allow_html=True)

        st.markdown("real space grid for lead1: %d %d %d"%( nx_lead1, ny1, nz1))
        st.markdown("real space grid for lead2: %d %d %d"%( nx_lead2, ny1, nz1))
        st.markdown("real space grid for center: %d %d %d"%( nx_center, ny1, nz1))

        pot_grid= col2.number_input("rho pot grid refinement", value=1)
        grid_lines ='#******** REAL SPACE GRID ********   \n'
        grid_lines += 'potential_grid_refinement="%d"  \n'%pot_grid
        grid_lines += '  \n'

    return grid_lines, nx_lead1, nx_lead2, nx_center, ny1, nz1, a_lead1, a_lead2, a_center, eq_left_right
