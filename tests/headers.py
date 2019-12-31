# A very simple example header
header1 = """
SDDS1
!# little-endian
&parameter name=par1, type=double, &end
&parameter name=par2, type=long, &end
&column name=col1, type=double,  &end
&column name=col2, type=double, units="m", symbol="&n", description="A test description",  &end
&data mode=ascii, &end
"""

# Header for typical .sig file from elegant
header2 = """
SDDS1
!# little-endian
&description text="sigma matrix--input: elegant.ele  lattice: elegant.lte", contents="sigma matrix", &end
&parameter name=Step, description="Simulation step", type=long, &end
&column name=s, units=m, description=Distance, type=double,  &end
&column name=ElementName, description="Element name", format_string=%10s, type=string,  &end
&column name=ElementOccurence, description="Occurence of element", format_string=%6ld, type=long,  &end
&column name=ElementType, description="Element-type name", format_string=%10s, type=string,  &end
&column name=s1, symbol="$gs$r$b1$n", units=m, description="sqrt(<x*x>)", type=double,  &end
&column name=s12, symbol="$gs$r$b12$n", units=m, description="<x*xp'>", type=double,  &end
&column name=s13, symbol="$gs$r$b13$n", units="m$a2$n", description="<x*y>", type=double,  &end
&column name=s14, symbol="$gs$r$b14$n", units=m, description="<x*y'>", type=double,  &end
&column name=s15, symbol="$gs$r$b15$n", units="m$a2$n", description="<x*s>", type=double,  &end
&column name=s16, symbol="$gs$r$b16$n", units=m, description="<x*delta>", type=double,  &end
&column name=s17, symbol="$gs$r$b17$n", units="m*s", description="<x*t>", type=double,  &end
&column name=s2, symbol="$gs$r$b2$n", description="sqrt(<x'*x'>)", type=double,  &end
&column name=s23, symbol="$gs$r$b23$n", units=m, description="<x'*y>", type=double,  &end
&column name=s24, symbol="$gs$r$b24$n", description="<x'*y'>", type=double,  &end
&column name=s25, symbol="$gs$r$b25$n", units=m, description="<x'*s>", type=double,  &end
&column name=s26, symbol="$gs$r$b26$n", description="<x'*delta>", type=double,  &end
&column name=s27, symbol="$gs$r$b27$n", units=s, description="<x'*t>", type=double,  &end
&column name=s3, symbol="$gs$r$b3$n", units=m, description="sqrt(<y*y>)", type=double,  &end
&column name=s34, symbol="$gs$r$b34$n", units=m, description="<y*y'>", type=double,  &end
&column name=s35, symbol="$gs$r$b35$n", units="m$a2$n", description="<y*s>", type=double,  &end
&column name=s36, symbol="$gs$r$b36$n", units=m, description="<y*delta>", type=double,  &end
&column name=s37, symbol="$gs$r$b37$n", units="m*s", description="<y*t>", type=double,  &end
&column name=s4, symbol="$gs$r$b4$n", description="sqrt(<y'*y')>", type=double,  &end
&column name=s45, symbol="$gs$r$b45$n", units=m, description="<y'*s>", type=double,  &end
&column name=s46, symbol="$gs$r$b46$n", description="<y'*delta>", type=double,  &end
&column name=s47, symbol="$gs$r$b47$n", units=s, description="<y'*t>", type=double,  &end
&column name=s5, symbol="$gs$r$b5$n", units=m, description="sqrt(<s*s>)", type=double,  &end
&column name=s56, symbol="$gs$r$b56$n", units=m, description="<s*delta>", type=double,  &end
&column name=s57, symbol="$gs$r$b57$n", units="m*s", description="<s*t>", type=double,  &end
&column name=s6, symbol="$gs$r$b6$n", description="sqrt(<delta*delta>)", type=double,  &end
&column name=s67, symbol="$gs$r$b67$n", units=s, description="<delta*t>", type=double,  &end
&column name=s7, symbol="$gs$r$b7$n", description="sqrt(<t*t>)", type=double,  &end
&column name=ma1, symbol="max$sb$ex$sb$e", units=m, description="maximum absolute value of x", type=double,  &end
&column name=ma2, symbol="max$sb$ex'$sb$e", description="maximum absolute value of x'", type=double,  &end
&column name=ma3, symbol="max$sb$ey$sb$e", units=m, description="maximum absolute value of y", type=double,  &end
&column name=ma4, symbol="max$sb$ey'$sb$e", description="maximum absolute value of y'", type=double,  &end
&column name=ma5, symbol="max$sb$e$gD$rs$sb$e", units=m, description="maximum absolute deviation of s", type=double,  &end
&column name=ma6, symbol="max$sb$e$gd$r$sb$e", description="maximum absolute value of delta", type=double,  &end
&column name=ma7, symbol="max$sb$e$gD$rt$sb$e", units=s, description="maximum absolute deviation of t", type=double,  &end
&column name=minimum1, symbol="x$bmin$n", units=m, type=double,  &end
&column name=minimum2, symbol="x'$bmin$n", units=m, type=double,  &end
&column name=minimum3, symbol="y$bmin$n", units=m, type=double,  &end
&column name=minimum4, symbol="y'$bmin$n", units=m, type=double,  &end
&column name=minimum5, symbol="$gD$rs$bmin$n", units=m, type=double,  &end
&column name=minimum6, symbol="$gd$r$bmin$n", units=m, type=double,  &end
&column name=minimum7, symbol="t$bmin$n", units=s, type=double,  &end
&column name=maximum1, symbol="x$bmax$n", units=m, type=double,  &end
&column name=maximum2, symbol="x'$bmax$n", units=m, type=double,  &end
&column name=maximum3, symbol="y$bmax$n", units=m, type=double,  &end
&column name=maximum4, symbol="y'$bmax$n", units=m, type=double,  &end
&column name=maximum5, symbol="$gD$rs$bmax$n", units=m, type=double,  &end
&column name=maximum6, symbol="$gd$r$bmax$n", units=m, type=double,  &end
&column name=maximum7, symbol="t$bmax$n", units=s, type=double,  &end
&column name=Sx, symbol="$gs$r$bx$n", units=m, description=sqrt(<(x-<x>)^2>), type=double,  &end
&column name=Sxp, symbol="$gs$r$bx'$n", description=sqrt(<(x'-<x'>)^2>), type=double,  &end
&column name=Sy, symbol="$gs$r$by$n", units=m, description=sqrt(<(y-<y>)^2>), type=double,  &end
&column name=Syp, symbol="$gs$r$by'$n", description=sqrt(<(y'-<y'>)^2>), type=double,  &end
&column name=Ss, symbol="$gs$r$bs$n", units=m, description=sqrt(<(s-<s>)^2>), type=double,  &end
&column name=Sdelta, symbol="$gs$bd$n$r", description=sqrt(<(delta-<delta>)^2>), type=double,  &end
&column name=St, symbol="$gs$r$bt$n", units=s, description=sqrt(<(t-<t>)^2>), type=double,  &end
&column name=ex, symbol="$ge$r$bx$n", units=m, description="geometric horizontal emittance", type=double,  &end
&column name=enx, symbol="$ge$r$bx,n$n", units=m, description="normalized horizontal emittance", type=double,  &end
&column name=ecx, symbol="$ge$r$bx,c$n", units=m, description="geometric horizontal emittance less dispersive contributions", type=double,  &end
&column name=ecnx, symbol="$ge$r$bx,cn$n", units=m, description="normalized horizontal emittance less dispersive contributions", type=double,  &end
&column name=ey, symbol="$ge$r$by$n", units=m, description="geometric vertical emittance", type=double,  &end
&column name=eny, symbol="$ge$r$by,n$n", units=m, description="normalized vertical emittance", type=double,  &end
&column name=ecy, symbol="$ge$r$by,c$n", units=m, description="geometric vertical emittance less dispersive contributions", type=double,  &end
&column name=ecny, symbol="$ge$r$by,cn$n", units=m, description="normalized vertical emittance less dispersive contributions", type=double,  &end
&column name=betaxBeam, symbol="$gb$r$bx,beam$n", units=m, description="betax for the beam, excluding dispersive contributions", type=double,  &end
&column name=alphaxBeam, symbol="$ga$r$bx,beam$n", description="alphax for the beam, excluding dispersive contributions", type=double,  &end
&column name=betayBeam, symbol="$gb$r$by,beam$n", units=m, description="betay for the beam, excluding dispersive contributions", type=double,  &end
&column name=alphayBeam, symbol="$ga$r$by,beam$n", description="alphay for the beam, excluding dispersive contributions", type=double,  &end
&data mode=binary, &end
"""

# Header for a .stat file from OPAL 2.0.1
header3 = """
SDDS1
&description
        text="Statistics data 'spectrometer.in' 31/12/2019 18:03:12",
        contents="stat parameters"
&end
&parameter
        name=processors,
        type=long,
        description="Number of Cores used"
&end
&parameter
        name=revision,
        type=string,
        description="git revision of opal"
&end
&parameter
        name=flavor,
        type=string,
        description="OPAL flavor that wrote file"
&end
&column
        name=t,
        type=double,
        units=ns,
        description="1 Time"
&end
&column
        name=s,
        type=double,
        units=m,
        description="2 Average Longitudinal Position"
&end
&column
        name=numParticles,
        type=long,
        units=1,
        description="3 Number of Macro Particles"
&end
&column
        name=charge,
        type=double,
        units=1,
        description="4 Bunch Charge"
&end
&column
        name=energy,
        type=double,
        units=MeV,
        description="5 Mean Energy"
&end
&column
        name=rms_x,
        type=double,
        units=m,
        description="6 RMS Beamsize in x"
&end
&column
        name=rms_y,
        type=double,
        units=m,
        description="7 RMS Beamsize in y"
&end
&column
        name=rms_s,
        type=double,
        units=m,
        description="8 RMS Beamsize in s"
&end
&column
        name=rms_px,
        type=double,
        units=1,
        description="9 RMS Momenta in x"
&end
&column
        name=rms_py,
        type=double,
        units=1,
        description="10 RMS Momenta in y"
&end
&column
        name=rms_ps,
        type=double,
        units=1,
        description="11 RMS Momenta in s"
&end
&column
        name=emit_x,
        type=double,
        units=m,
        description="12 Normalized Emittance x"
&end
&column
        name=emit_y,
        type=double,
        units=m,
        description="13 Normalized Emittance y"
&end
&column
        name=emit_s,
        type=double,
        units=m,
        description="14 Normalized Emittance s"
&end
&column
        name=mean_x,
        type=double,
        units=m,
        description="15 Mean Beam Position in x"
&end
&column
        name=mean_y,
        type=double,
        units=m,
        description="16 Mean Beam Position in y"
&end
&column
        name=mean_s,
        type=double,
        units=m,
        description="17 Mean Beam Position in s"
&end
&column
        name=ref_x,
        type=double,
        units=m,
        description="18 x coordinate of reference particle in lab cs"
&end
&column
        name=ref_y,
        type=double,
        units=m,
        description="19 y coordinate of reference particle in lab cs"
&end
&column
        name=ref_z,
        type=double,
        units=m,
        description="20 z coordinate of reference particle in lab cs"
&end
&column
        name=ref_px,
        type=double,
        units=1,
        description="21 x momentum of reference particle in lab cs"
&end
&column
        name=ref_py,
        type=double,
        units=1,
        description="22 y momentum of reference particle in lab cs"
&end
&column
        name=ref_pz,
        type=double,
        units=1,
        description="23 z momentum of reference particle in lab cs"
&end
&column
        name=max_x,
        type=double,
        units=m,
        description="24 Max Beamsize in x"
&end
&column
        name=max_y,
        type=double,
        units=m,
        description="25 Max Beamsize in y"
&end
&column
        name=max_s,
        type=double,
        units=m,
        description="26 Max Beamsize in s"
&end
&column
        name=xpx,
        type=double,
        units=1,
        description="27 Correlation xpx"
&end
&column
        name=ypy,
        type=double,
        units=1,
        description="28 Correlation ypyy"
&end
&column
        name=zpz,
        type=double,
        units=1,
        description="29 Correlation zpz"
&end
&column
        name=Dx,
        type=double,
        units=m,
        description="30 Dispersion in x"
&end
&column
        name=DDx,
        type=double,
        units=1,
        description="31 Derivative of dispersion in x"
&end
&column
        name=Dy,
        type=double,
        units=m,
        description="32 Dispersion in y"
&end
&column
        name=DDy,
        type=double,
        units=1,
        description="33 Derivative of dispersion in y"
&end
&column
        name=Bx_ref,
        type=double,
        units=T,
        description="34 Bx-Field component of ref particle"
&end
&column
        name=By_ref,
        type=double,
        units=T,
        description="35 By-Field component of ref particle"
&end
&column
        name=Bz_ref,
        type=double,
        units=T,
        description="36 Bz-Field component of ref particle"
&end
&column
        name=Ex_ref,
        type=double,
        units=MV/m,
        description="37 Ex-Field component of ref particle"
&end
&column
        name=Ey_ref,
        type=double,
        units=MV/m,
        description="38 Ey-Field component of ref particle"
&end
&column
        name=Ez_ref,
        type=double,
        units=MV/m,
        description="39 Ez-Field component of ref particle"
&end
&column
        name=dE,
        type=double,
        units=MeV,
        description="40 energy spread of the beam"
&end
&column
        name=dt,
        type=double,
        units=ns,
        description="41 time step size"
&end
&column
        name=partsOutside,
        type=double,
        units=1,
        description="42 outside n*sigma of the beam"
&end
&data
        mode=ascii,
        no_row_counts=1
&end
"""