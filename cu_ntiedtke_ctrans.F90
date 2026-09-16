!==================================================================================================================
 module cu_ntiedtke_ctrans
 use ccpp_kind_types,only: kind_phys
 use cu_ntiedtke_common,only: cmfcmin
 use mpas_log, only : mpas_log_write

 implicit none
 private
 public:: cu_ntiedtke_ctrans_run,     &
          cu_ntiedtke_ctrans_init,    &
          cu_ntiedtke_ctrans_finalize


 contains


!==================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_init
!!\html\include cu_ntiedtke_ctrans_init.html
!!
 subroutine cu_ntiedtke_ctrans_init(errmsg,errflg)
!==================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!------------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_init

!==================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_finalize
!!\html\include cu_ntiedtke_ctrans_finalize.html
!!
 subroutine cu_ntiedtke_ctrans_finalize(errmsg,errflg)
!==================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!------------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_finalize

!==================================================================================================================
 subroutine cu_ntiedtke_ctrans_run(klon, klev, nchem, ldcum, lddraf, kctype, kcbot, kctop, kdtop, &
                                   grav, ztmst, do_scav, fscav, is_aerosol, scav_data_gas, &
                                   chem, thv, ptenc, ghti, paph, pmfu, pmfd, pmfude_rate, &
                                   pmfdde_rate, errmsg, errflg, itimestep, rate_incloud)
!==================================================================================================================

!--- input arguments:
 logical,intent(in):: do_scav

 integer,intent(in):: klon,klev
 integer,intent(in):: nchem

 integer,intent(in):: itimestep

 real(kind=kind_phys),intent(in),dimension(klon,klev):: thv  ! ADDED: 2D Temperature slice
 logical,intent(in),dimension(nchem):: is_aerosol
 real(kind=kind_phys),intent(in),dimension(6,nchem):: scav_data_gas
 real(kind=kind_phys), intent(inout), dimension(klon, nchem) :: rate_incloud

 integer,intent(in),dimension(klon):: kctype,kcbot,kctop,kdtop

 logical,intent(in),dimension(klon):: ldcum,lddraf

 real(kind=kind_phys),intent(in):: grav,ztmst
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfu,pmfude_rate
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfd,pmfdde_rate
 real(kind=kind_phys),intent(in),dimension(klon,klev+1):: paph,ghti

 real(kind=kind_phys),intent(in),dimension(nchem):: fscav
 real(kind=kind_phys),intent(in),dimension(klon,klev,nchem):: chem

!--- inout arguments:
 real(kind=kind_phys),dimension(klon,klev,nchem):: ptenc

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).
 integer,intent(out):: &
    errflg      ! output error flag (-).

!--- local variables and arrays:
 logical,dimension(klon):: lldcum,llddraf3

 integer:: ic,jl,jk

 real(kind=kind_phys):: zcons,zmfmax
 real(kind=kind_phys),dimension(klon):: zmfs
 real(kind=kind_phys),dimension(klon,klev):: zmfuus,zmfudr
 real(kind=kind_phys),dimension(klon,klev):: zmfdus,zmfddr

!------------------------------------------------------------------------------------------------------------------

 zcons=1./(grav*ztmst)

 do ic = 1,nchem
    ptenc(:,:,ic) = 0._kind_phys
 enddo


!--- no convective transport for mid-level convection:
 do jl = 1,klon
    if(ldcum(jl) .and. kctype(jl) /= 3 .and. kcbot(jl)-kctop(jl) >= 1 ) then
       lldcum(jl)   = .true.
       llddraf3(jl) = lddraf(jl)
    else
       lldcum(jl) = .false.
       llddraf3(jl) = .false.
    endif
 enddo


!--- check and correct mass fluxes for CFL criterium:
 zmfs(:) = 1.
 do jk = 2,klev
    do jl = 1,klon
       if(lldcum(jl) .and. jk >= kctop(jl) ) then
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*0.8*zcons

          if(pmfu(jl,jk) > zmfmax) then
             zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
          endif
       endif
    enddo
 enddo

 do jk = 1,klev
    do jl = 1,klon
       if(lldcum(jl) .and. jk >= kctop(jl)-1) then
          zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
          zmfudr(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
       else
          zmfuus(jl,jk) = 0.
          zmfudr(jl,jk) = 0.
       endif

       if(llddraf3(jl) .and. jk >= kdtop(jl)-1) then
          zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          zmfddr(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
       else
          zmfdus(jl,jk) = 0.
          zmfddr(jl,jk) = 0.
       endif
    enddo
 enddo


!--- call to subroutine that computes the convective transport of chemical species:
 call cuctracer(klon, klev, nchem, kctop, kdtop, ldcum, lddraf, grav, ztmst, do_scav, &
                ghti, paph, thv, zmfuus, zmfdus, zmfudr, zmfddr, chem, fscav, &
                is_aerosol, scav_data_gas, ptenc, itimestep, rate_incloud)

!--- output error flag and message:
 errflg = 0
 errmsg = " "


 end subroutine cu_ntiedtke_ctrans_run

!==================================================================================================================
 subroutine cuctracer(klon, klev, ktrac, kctop, kdtop, ldcum, lddraf, grav, ztmst, do_scav, &
                      ght, paph, thv, pmfu, pmfd, pudrate, pddrate, pcen, fscav, &
                      is_aerosol, scav_data_gas, ptenc, itimestep, rate_incloud)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: klon,klev,ktrac
 integer,intent(in),dimension(klon):: kctop,kdtop

 logical,intent(in):: do_scav
 logical,intent(in),dimension(klon):: ldcum,lddraf

 real(kind=kind_phys),intent(in):: grav,ztmst
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfu,pudrate
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfd,pddrate
 real(kind=kind_phys),intent(in),dimension(klon,klev+1):: ght,paph
 real(kind=kind_phys),intent(in),dimension(klon,klev,ktrac):: pcen

 real(kind=kind_phys),intent(in),dimension(klon,klev):: thv ! Local Temperature
 real(kind=kind_phys),intent(in),dimension(ktrac):: fscav
 logical,intent(in),dimension(ktrac):: is_aerosol
 real(kind=kind_phys),intent(in),dimension(6,ktrac):: scav_data_gas
 real(kind=kind_phys),intent(inout),dimension(klon,ktrac):: rate_incloud ! NEW: [kg/m2/s]
 real(kind=kind_phys) :: scav_rate
 real(kind=kind_phys) :: dz
 integer, intent(in)  :: itimestep
 real(kind=kind_phys) :: tfac, h_phys, h_star, k1, k2, h_ion, scav_eff

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(klon,klev,ktrac):: ptenc

!--- variables and arrays:
 integer:: ik,jk,jl,jn

 logical,dimension(klon,klev):: llcumask,llcumbas

 real(kind=kind_phys):: zzp,zmfa,zerate,zposi, chem_before
 real(kind=kind_phys),dimension(klon,klev):: zdp
 real(kind=kind_phys),dimension(klon,klev,ktrac):: zcen,zcu,zcd,zmfc,ztenc,sink


!------------------------------------------------------------------------------------------------------------------

!--- initialization:
 do jk = 2,klev
    do jl = 1,klon
        llcumask(jl,jk) = .false.
        if(ldcum(jl)) then
           zdp(jl,jk) = grav/(paph(jl,jk+1)-paph(jl,jk))
           if( jk >= kctop(jl)-1) llcumask(jl,jk) = .true.
        endif
    enddo
 enddo

! Rain pH = 5.0
 h_ion = 1.0e-5_kind_phys
 rate_incloud(:,:) = 0.0_kind_phys
 sink = 0.0_kind_phys
 zcen = pcen
 ztenc = 0.0_kind_phys

!--- loop over all chemical species:
 do jn = 1,ktrac
    !define chemical species at half levels:
    do jk = 2,klev
       ik = jk-1
       do jl = 1,klon
          zcen(jl,jk,jn) = pcen(jl,jk,jn)
          zcd(jl,jk,jn)  = pcen(jl,ik,jn)
          zcu(jl,jk,jn)  = pcen(jl,ik,jn)
          zmfc(jl,jk,jn) = 0.
          ztenc(jl,jk,jn)= 0.
       enddo
    enddo
    do jl = 1,klon
       zcu(jl,klev,jn) = pcen(jl,klev,jn)
    enddo

    ! Ascend from the lowest interface, passing scavenged air to the next level.
    do jk = klev,2,-1
       ik = min(jk+1,klev)
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             if(jk < klev .and. jk >= 3 .and. jk >= kctop(jl)) then
                zerate = pmfu(jl,jk)-pmfu(jl,ik)+pudrate(jl,jk)
                zmfa = 1./max(cmfcmin,pmfu(jl,jk))
                zcu(jl,jk,jn) = (pmfu(jl,ik)*zcu(jl,ik,jn) + &
                     zerate*pcen(jl,jk,jn)-pudrate(jl,jk)*zcu(jl,ik,jn))*zmfa
             endif
             if(do_scav .and. pmfu(jl,jk)>0.0_kind_phys) then
                dz = abs(ght(jl, jk+1) - ght(jl, jk))  ! layer thickness

                ! Define scavenging efficiency based on your species-specific fscav for aerosols and HLC for gases
                if (is_aerosol(jn)) then
                   scav_eff = fscav(jn)
                else
                ! THERMODYNAMIC BRANCH: Calculate H_star at local temperature
                   tfac = (1.0_kind_phys / thv(jl,jk)) - (1.0_kind_phys / 298.15_kind_phys)
                   h_phys = scav_data_gas(1, jn) * exp(scav_data_gas(2, jn) * tfac)
                   
                   if (scav_data_gas(3, jn) > 0.0_kind_phys) then
                      k1 = scav_data_gas(3, jn) * exp(scav_data_gas(4, jn) * tfac)
                      if (scav_data_gas(5, jn) > 0.0_kind_phys) then
                         k2 = scav_data_gas(5, jn) * exp(scav_data_gas(6, jn) * tfac)
                         
                         ! NH3 (Base) vs SO2 (Acid)
                         if (scav_data_gas(1, jn) > 10.0_kind_phys) then
                            h_star = h_phys * (1.0_kind_phys + (k1 * h_ion) / k2) ! NH3
                         else
                            h_star = h_phys * (1.0_kind_phys + k1/h_ion + (k1*k2)/(h_ion**2)) ! SO2
                         end if
                      else
                         h_star = h_phys * (1.0_kind_phys + k1/h_ion) ! MSA
                      end if
                   else
                      h_star = h_phys ! DMS
                   end if

               ! --- FULL SPECTRUM LOG-LINEAR RAMP ---
                   ! Interp log10(H*) from -1.0 to 6.0 across 7 orders of magnitude
                   if (h_star <= 0.1_kind_phys) then
                      scav_eff = 0.0_kind_phys
                   else if (h_star >= 1000000.0_kind_phys) then
                      scav_eff = 1.0_kind_phys
                   else
                      scav_eff = (log10(max(1.0e-10_kind_phys, h_star)) + 1.0_kind_phys) / 7.0_kind_phys
                   end if
                end if

                ! Since fscav has units of 1/km, we multiply it with dz/1000. dz is in meters.
                scav_rate = scav_eff * (dz/1000.0_kind_phys)
                chem_before = zcu(jl,jk,jn) ! Store for diagnostic use
                
                ! Apply exponential decay to the updraft concentration (zcu)
                zcu(jl,jk,jn) = zcu(jl,jk,jn) * exp(-min(scav_rate, 10.0_kind_phys))

                ! Flux [kg/m2/s] = Updraft Air Mass Flux [kg_air/m2/s] * Change in mixing ratio [kg_chem/kg_air]
                sink(jl,jk,jn) = max(0.0_kind_phys,pmfu(jl,jk)) * (chem_before-zcu(jl,jk,jn))
                rate_incloud(jl,jn) = rate_incloud(jl,jn) + sink(jl,jk,jn)

             endif
          endif
       enddo
    enddo

    !compute downdraft values:
    do jk = 3,klev
       ik = jk - 1
       do jl = 1,klon
          if(lddraf(jl) .and. jk == kdtop(jl)) then
             !note: in order to avoid final negative tracer values at LFS
             !the allowed value of ZCD depends on the jump in mass flux
             !at the LFS:
             zcd(jl,jk,jn) = 0.1*zcu(jl,jk,jn) + 0.9*pcen(jl,ik,jn)
          elseif(lddraf(jl).and.jk>kdtop(jl)) then
             zerate = -pmfd(jl,jk) + pmfd(jl,ik) + pddrate(jl,jk)
             zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
             zcd(jl,jk,jn) = (pmfd(jl,ik)*zcd(jl,ik,jn) - &
             zerate*pcen(jl,ik,jn)+pddrate(jl,jk)*zcd(jl,ik,jn))*zmfa
          endif
       enddo
    enddo

    !in order to avoid negative tracer at KLEV, then adjust ZCD:
    jk = klev
    ik = jk - 1
    do jl = 1,klon
       if(lddraf(jl)) then
          zposi = -zdp(jl,jk) *(pmfu(jl,jk)*zcu(jl,jk,jn) + &
                  pmfd(jl,jk)*zcd(jl,jk,jn)-(pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn))
          if(pcen(jl,jk,jn)+zposi*ztmst < 0.) then
             zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
             zcd(jl,jk,jn) = ((pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn) - &
                  pmfu(jl,jk)*zcu(jl,jk,jn)+pcen(jl,jk,jn) / &
                  (ztmst*zdp(jl,jk)))*zmfa
          endif
       endif
    enddo
 enddo


 do jn = 1,ktrac

    !compute fluxes:
    do jk = 2,klev
       ik = jk - 1
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             zmfa = pmfu(jl,jk) + pmfd(jl,jk)
             zmfc(jl,jk,jn) = pmfu(jl,jk)*zcu(jl,jk,jn) + &
                              pmfd(jl,jk)*zcd(jl,jk,jn) - zmfa*zcen(jl,ik,jn)
          endif
       enddo
    enddo

    !compute tendencies:
    do jk = 2 , klev - 1
       ik = jk + 1
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             ztenc(jl,jk,jn) = zdp(jl,jk)*(zmfc(jl,ik,jn)-zmfc(jl,jk,jn))
          endif
       enddo
    enddo
    jk = klev
    do jl = 1,klon
       if(ldcum(jl)) ztenc(jl,jk,jn) = -zdp(jl,jk)*zmfc(jl,jk,jn)
    enddo
 enddo


 do jn = 1,ktrac
    !update tendencies:
    do jk = 2,klev
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             ! Flux divergence alone conserves mass; precipitated mass is a separate sink.
             ptenc(jl,jk,jn) = ptenc(jl,jk,jn)+ztenc(jl,jk,jn)-zdp(jl,jk)*sink(jl,jk,jn)
          endif
       enddo
    enddo
 enddo

 end subroutine cuctracer

!==================================================================================================================
 end module cu_ntiedtke_ctrans
!==================================================================================================================
