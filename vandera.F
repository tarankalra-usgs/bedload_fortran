      PROGRAM VANDERA 

      USE mod_kinds
      USE mod_scalars 
!
      implicit none 

      integer  :: iunit
      real(r8) :: Hs, Td, depth, d50, d90
      real(r8) :: urms, umag_curr, phi_curwave 
      real(r8) :: y, kh, uhat, ahat 
      real(r8) :: k_wn, c_w
      real(r8) :: rhos 
      real(r8) :: smgd, osmgd
      real(r8) :: bed_frac 
!
      real(r8) :: r, phi, Su, Au
      real(r8) :: Sk, Ak
      real(r8) :: T_cu, T_tu, T_t, T_c
      real(r8) :: umax, umin, RR
!
      real(r8) :: uc_r, ut_r
      real(r8) :: uhat_c, uhat_t
      real(r8) :: udelta, mag_uc, mag_ut
!
      real(r8) :: fd, ksw, eta, alpha, tau_wRe                
      real(r8) :: dsf_c, dsf_t               
      real(r8) :: theta_c, theta_t
      real(r8) :: theta_cx, theta_cy, theta_tx, theta_ty
      real(r8) :: mag_theta_c, mag_theta_t
!
      real(r8) :: om_cc, om_tt, om_ct, om_tc
!        
      real(r8) :: cff, cff1, cff2, cff3   
      real(r8) :: smgd_3
      real(r8) :: bedld_cx, bedld_cy
      real(r8) :: bedld_tx, bedld_ty
      real(r8) :: bedld_x, bedld_y
     
      real(r8) :: wavecycle 
      real(r8) :: Zref, delta
!            
!     Input wave parameters 
!     Hs=significant wave height
!     Td=dominant wave period 
!     depth=water depth(m)
!   
      Hs = 1.6_r8
      Td = 6.5_r8
      depth =3.5_r8
      d50 = 0.245e-3
      d90 = 1.3_r8*d50 
      rhos= 2650.0_r8
! Zref is a reference height for computing current friction factor 
      Zref=0.04_r8
! delta is the reference height at which current velocity is computed (Wave boundary layer thickness) 
      delta=0.2_r8
!
! umag_curr is the current velocity magnitude (MAKE Sure its the magnitude)
! direction of currents is counter-clockwise from wave direction
!
      umag_curr=0.0_r8! ABS(0.5_r8)
      phi_curwave=0.0_r8*deg2rad
!
! uhat and ahat for the entire wave cycle  uhat is the wave orbital velocity defiend for the entire wave cycle
!
      urms=0.65_r8 
      uhat=SQRT(2.0_r8)*urms !0.5472_r8           
!
      ahat=uhat*Td/(2.0_r8*pi)
      k_wn=kh(Td,depth)/depth     ! Wave number 
      c_w=2*pi/(k_wn*Td)          ! Wave speed
! 
! Open the file to write code output
!
      open(iunit,file='vandera_output.dat',form='formatted')
      write(iunit,*)"****************************************"
      write(iunit,*)"Write input data"
      write(iunit,*)"****************************************"
      write(iunit,*) "Signfcnt wave height(m)=", Hs 
      write(iunit,*) "Surface wave period(s)=", Td 
      write(iunit,*) "Water depth(m)=", Depth
      write(iunit,*) "Current velocity(m/s)=", umag_curr
      write(iunit,*) "Current dir(degrees)=",phi_curwave*rad2deg
      write(iunit,*) "Sediment d50(m)=", d50
      write(iunit,*)"****************************************"
      write(iunit,*)"Write output data"
      write(iunit,*)"****************************************"
!
! VA-2013 equation 1 is solved in 3 sub-steps
!
!----------------------------------------------------------------------
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "bb" and "r"
! COMMON TO both crest and trough
!-----------------------------------------------------------------------
!
      CALL skewness_params(iunit, Hs, Td, depth, r, phi)
!        
!-----------------------------------------------------------------------
! Abreu et al. use skewness params to get representative critical orbital
! velocity for crest and trough cycles , get r and phi from above
!-----------------------------------------------------------------------
! 
      CALL abreu_points(iunit, r, phi, uhat, Td,  T_c, T_t,             &
     &                                    T_cu, T_tu, umax, umin, RR)
      write(iunit,*)"****************************************"
      write(iunit,*)"umax=", umax,"umin=",umin
      write(iunit,*)"****************************************"
!
!-----------------------------------------------------------------------
!           Crest half cycle
!-----------------------------------------------------------------------
! Get the "representative crest half cycle water particle velocity
!    as well as full cycle orbital velocity and excursion 
!-----------------------------------------------------------------------
!
! from Abreu points 
!
      uhat_c=umax
      uhat_t=umin 
!
!-----------------------------------------------------------------------
! VA2013 Equation 10, 11
!-----------------------------------------------------------------------
!
      uc_r=0.5_r8*sqrt(2.0_r8)*uhat_c
      ut_r=0.5_r8*sqrt(2.0_r8)*uhat_t
!
      smgd=(rhos/rho0-1.0_r8)*g*d50
      osmgd=1.0_r8/smgd
! 
! Full wave cycle 
! 
      CALL full_wave_cycle_stress_factors( d50, d90, osmgd,             &
     &                                                 Td, depth,       &
     &                                    umag_curr, phi_curwave,       &
     &                                              Zref,  delta,       &
     &                                            RR, uhat, ahat,       &
     &                                                umax, umin,       &
     &                                      T_c, T_t, T_cu, T_tu,       &
     &                           eta, udelta, alpha, fd, ksw, tau_wRe )
!
      write(iunit,*)"****************************************"
      write(iunit,*)" Full-wave cycle output"
      write(iunit,*)"****************************************"
      write(iunit,*) "Current velocity at the wave BL.",  udelta
      write(iunit,*) "Wave Reynolds Stress tau_wRe", tau_wRe
      write(iunit,*) "Ripple height (m)",  eta
      write(iunit,*)"****************************************"
      write(iunit,*)" Full-wave cycle output modified time periods"
      write(iunit,*) "Time period in crest cycle",  T_c, T_cu
      write(iunit,*) "Time period in trough cycle", T_t, T_tu
      write(iunit,*)"****************************************"
!
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Crest half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
      CALL half_wave_cycle_stress_factors( T_cu, T_c, ahat, ksw,        &
     &                                           fd, alpha,             &
     &                                           d50, osmgd,            &
     &                                uc_r, uhat_c, udelta, phi_curwave,&
     &                                           tau_wRe,               &
     &                          dsf_c, theta_cx, theta_cy, mag_theta_c )
!
      write(iunit,*)"****************************************"
      write(iunit,*)"2. CREST cycle - Bed shear stress output"
      write(iunit,*)"****************************************"
      write(iunit,*) "Shields par. x,y dir", theta_cx, theta_cy
      write(iunit,*) "Magnitude of Shields par.", mag_theta_c
      write(iunit,*) "Sheet flow thickness (m)", dsf_c
      write(iunit,*)"****************************************"
!    
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each crest half cycle
!-----------------------------------------------------------------------
!
      wavecycle=1.0_r8
      CALL sandload_vandera(iunit, wavecycle,                           &
     &                              Hs, Td,  depth, RR,                 &
     &                              d50, rhos, c_w,                     &
     &                              eta, dsf_c,                         &
     &                              T_c, T_cu, uhat_c, mag_theta_c,     &
     &                              om_cc, om_tc)
!
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Trough half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
      CALL half_wave_cycle_stress_factors( T_tu, T_t, ahat, ksw,        &
     &                                           fd, alpha,             &
     &                                           d50, osmgd,            &
     &                                ut_r, uhat_t, udelta, phi_curwave,&
     &                                           tau_wRe,               &
     &                          dsf_t, theta_tx, theta_ty, mag_theta_t )
!
      write(iunit,*)"****************************************"
      write(iunit,*)"2.TROUGH cycle - Bed shear stress output"
      write(iunit,*)"****************************************"
      write(iunit,*) "Shields par. x,y dir", theta_tx, theta_ty
      write(iunit,*) "Magnitude of Shields par.", mag_theta_t
      write(iunit,*) "Sheet flow thickness (m)", dsf_t
      write(iunit,*)"****************************************"
!
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each trough half cycle
!-----------------------------------------------------------------------
!
      wavecycle=-1.0_r8
      CALL sandload_vandera(iunit, wavecycle,                           &
     &                              Hs, Td,  depth, RR,                 &
     &                              d50, rhos, c_w,                     &
     &                              eta, dsf_t,                         &
     &                              T_t, T_tu, uhat_t, mag_theta_t,     &
     &                              om_tt, om_ct)
!      
      cff1=0.5_r8*T_c/(T_cu)
      cff2=sqrt(mag_theta_c)*(om_cc+cff1*om_tc)
!
      cff3=theta_cx/mag_theta_c
      bedld_cx=cff2*cff3

      cff3=theta_cy/mag_theta_c
      bedld_cy=cff2*cff3

      cff1=0.5_r8*T_t/(T_tu)
      cff2=sqrt(mag_theta_t)*(om_tt+cff1*om_ct)
      
      cff3=theta_tx/mag_theta_t      
      bedld_tx=cff2*cff3
!
      cff3=theta_ty/mag_theta_t
      bedld_ty=cff2*cff3
!
! The units of these are m2 sec-1
! bed_frac, rhos multiplied 
!
      bed_frac=1.0_r8
!
!-----------------------------------------------------------------------
! VA2013  Use the velocity-load equation 1. 
! Non-dimensional net transport rate
!-----------------------------------------------------------------------
!
      smgd_3=sqrt((rhos/rho0-1.0_r8)*g*d50**3.0_r8)

      bedld_x=bed_frac*smgd_3*(bedld_cx*T_c+bedld_tx*T_t)/Td 
      bedld_y=bed_frac*smgd_3*(bedld_cy*T_c+bedld_ty*T_t)/Td 
!
      write(iunit,*)"****************************************"
      write(iunit,*)"Net transport rate in x dir.", bedld_x
      write(iunit,*)"Net transport rate in y dir.", bedld_y
      write(iunit,*)"****************************************"
      close(iunit)
      END PROGRAM VANDERA
 
      SUBROUTINE sandload_vandera(iunit, wavecycle,                     &
     &                              Hs, Td,  depth, RR,                 &
     &                              d50, rhos, c_w,                     &
     &                              eta, dsf,                           &
     &                              T_i, T_iu, uhat_i, mag_theta_i,     &
     &                              om_ii, om_iy)
!
      USE mod_kinds
      USE mod_scalars

      implicit none 

      integer, intent(in)  :: iunit 
      real(r8), intent(in) :: wavecycle 
      real(r8), intent(in) :: Hs, Td, depth, RR 
      real(r8), intent(in) :: d50, rhos, c_w
      real(r8), intent(in) :: eta, dsf
      real(r8), intent(in) :: T_i, T_iu
      real(r8), intent(in) :: uhat_i, mag_theta_i
      real(r8), intent(out):: om_ii, om_iy

      ! local variables
! 
! VA2013 Text under equation 37
! 
      real(r8), parameter :: m=11.0_r8, n=1.2_r8, alpha=8.2_r8
      real(r8), parameter :: xi=1.7_r8 ! Based on Santoss_core.m
      real(r8) :: eps_eff
      real(r8) :: om_i
      real(r8) :: theta_diff, theta_ieff, theta_cr, theta_cr_calc 
      real(r8) :: w_s, w_s_calc
      real(r8) :: ws_eta, ws_dsf, w_sc_calc
      real(r8) :: w_sc_eta, w_sc_dsf
      real(r8) :: cff, cff1_eta, cff1_dsf
      real(r8) :: P 
! 
! Find settling velocity based on Soulsby (1997). 
! VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
!
      w_s=w_s_calc(0.8_r8*d50, rhos)    
!
! VA2013 Equation 29, for crest cycle
!
      ws_eta=w_sc_calc(Hs, Td, depth, RR, w_s, eta)
      ws_dsf=w_sc_calc(Hs, Td, depth, RR, w_s, dsf)
      IF(wavecycle.eq.1.0_r8) THEN 
        w_sc_eta=MAX(w_s+ws_eta,0.0_r8)
        w_sc_dsf=MAX(w_s+ws_dsf,0.0_r8)
      ENDIF 
!
! VA2013 Equation 30, for trough cycle
!
      IF(wavecycle.eq.-1.0_r8) THEN 
        w_sc_eta=MAX(w_s-ws_eta,0.36*w_s)
        w_sc_dsf=MAX(w_s-ws_dsf,0.36*w_s)
      ENDIF 
!
! VA2013 Equation 33, Phase lag parameter
!
      cff=1.0_r8-(wavecycle*xi*uhat_i/c_w)
!
      cff1_eta=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_eta))
      cff1_dsf=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_dsf))
!
! For ripple regime 
!
      IF(eta.gt.0.0_r8) THEN 
        P=alpha*eta*cff*cff1_eta
        IF(wavecycle.eq.1.0_r8) THEN
         print*, "Ripple regime during crest"
        ELSE
         print*, "Ripple regime during trough"
        ENDIF 
      ELSEIF(eta.eq.0.0_r8)THEN
!
! For sheet flow regime 
!
        P=alpha*dsf*cff*cff1_dsf
        IF(wavecycle.eq.-1.0_r8) THEN
         print*, "Sheet flow regime during crest"
        ELSE
         print*, "Sheet flow regime during trough"
        ENDIF 
      ENDIF 
!
!      eps_eff=(dsf/d50)**0.25_r8 
      eps_eff=1.0_r8
      theta_ieff=eps_eff*mag_theta_i
! 
! Find critical Shields parameters based on Soulsby (1997).
!
      theta_cr=theta_cr_calc(d50, rhos)
!
! Sand load entrained in the flow during each half-cycle
!
      theta_diff=MAX((theta_ieff-theta_cr),0.0_r8)
      om_i=m*(theta_diff)**n
!
! VA2013 Equation 23-26, Sandload entrained during half cycle 
! 
              
      IF(P.le.1.0_r8) THEN 
        om_ii=om_i
        om_iy=0.0_r8
      ELSE
        om_ii=om_i/P 
        cff=1.0_r8/P
        om_iy=om_i*(1.0_r8-cff) 
      ENDIF 
!
      IF(wavecycle.eq.1.0_r8) THEN 
        write(iunit,*)"****************************************"
        write(iunit,*)"3. Sand load calculations"
        write(iunit,*)"****************************************"
        write(iunit,*)"Crest cycle"
        write(iunit,*)"****************************************"
        write(iunit,*)"Phase lag parameter", P
        write(iunit,*)"O_cc (sandload entrained and carried)", om_ii
        write(iunit,*)"O_tc (sandload carried in trough period)",om_iy
      ENDIF
      IF(wavecycle.eq.-1.0_r8) THEN 
        write(iunit,*)"****************************************"
        write(iunit,*)"Trough cycle"
        write(iunit,*)"****************************************"
        write(iunit,*)"Phase lag parameter", P
        write(iunit,*)"O_tt (sandload entrained and carried)", om_ii
        write(iunit,*)"O_ct (sandload carried in crest period)",om_iy
      ENDIF
      RETURN 
      END SUBROUTINE sandload_vandera
!
      REAL(r8) FUNCTION w_sc_calc(Hs, Td, depth, RR, w_s, zws)
! 
! Second order Stokes theory to get vertical velocity of water particle
! at a given elevation based on santoss_core.m
! 
      USE mod_kinds
      USE mod_scalars

      implicit none 
!     
      real(r8) :: Hs, Td, depth, RR, zws, w_s
      real(r8) :: cff, worb1, worb2, worb
!      
      worb1=pi*Hs*zws/(Td*depth)
      worb2=worb1*2.0_r8*(RR+RR-1.0_r8)
!
!  Using the SANTOSS model formulation 
! 
      cff=1.0_r8/8.0_r8
      worb=cff*worb1*SQRT(64.0_r8-(-worb1+                              &
     &          SQRT(worb1**2.0_r8+32.0_r8*                             &
     &               worb2**2.0_r8))**2.0_r8/(worb2**2.0_r8))+          &
     &               worb2*SIN(2.0_r8*ACOS(cff*(-worb1+                 &
     &          SQRT(worb1**2.0_r8+32.0_r8*worb2**2.0_r8))/worb2))
!
      w_sc_calc=worb
!
      RETURN       
      END FUNCTION w_sc_calc
!
      SUBROUTINE full_wave_cycle_stress_factors( d50, d90, osmgd,       &
     &                                                 Td, depth,       &
     &                              umag_curr, phi_curwave,             &
     &                                            Zref,  delta,         &
     &                                            RR, uhat, ahat,       &
     &                                                umax, umin,       &
     &                                      T_c, T_t, T_cu, T_tu,       &
     &                         eta, udelta, alpha, fd, ksw, tau_wRe) 
!
!**********************************************************************
!  This subroutine returns the following:
!  eta                 : ripple height
!  udelta              : current velocity at the wave boundary layer
!  fd                  : current friction factor  
!  tau_wRe             : Wave averaged Reynolds stress
!  T_c, T_t, T_cu, T_tu: Updated time periods in half cycles 
!                        based on current velocity
!**********************************************************************
! 
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!  
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles 
!       
      real(r8), intent(in) :: d50, d90, osmgd
      real(r8), intent(in) :: Td, depth
      real(r8), intent(in) :: umag_curr, phi_curwave, Zref, delta 
      real(r8), intent(in) :: RR, uhat, ahat
      real(r8), intent(in) :: umax, umin
      real(r8), intent(inout) :: T_c, T_t, T_cu, T_tu
      real(r8), intent(inout) :: eta, udelta, fd, ksw, tau_wRe
!
!  Local variables
! 
      integer  :: iter
      integer, parameter :: total_iters=10
      real, parameter :: tol=0.001_r8
      real(r8) :: mu, mu_calc 
      real(r8) :: theta_timeavg_old, theta_timeavg, theta_hat_i

      real(r8) :: psi
      real(r8) :: rlen  ! ripple dimensions
      real(r8) :: ksd, ksd_calc
      real(r8) :: ksw_calc
      real(r8) :: fw, fw_calc
      real(r8) :: fd_calc
      real(r8) :: fw_i, fwi_calc, alpha, fwd_i 
      real(r8) :: dsf_calc
      real(r8) :: alpha_w, fwd, kh, k_wn, c_w
      real(r8) :: ustarc
!
! Iterative solution to obtain current and wave related bed roughness
! VA2013 Apendix A, Shields parameter (Stress) depends on bed roughness 
! Bed roughness computed from converged Shields parameter
!
! Maximum mobility number at crest and trough 
! For irregular waves, use Rayleigh distributed maximum value
! VA, text under equation Appendix B.4 
!
      psi=(1.27_r8*uhat)**2*osmgd
!
! Use Appendix B eqn B.1 and B.2 to get ripple height and length 
!   
      CALL ripple_dim(psi, d50, eta, rlen)
!
      eta=eta*ahat
      rlen=rlen*ahat

      alpha=umag_curr/(umag_curr+uhat)
!
! Initiliaze with theta_timeavg=0 and theta_hat_i=theta_timeavg
!
      theta_timeavg=0.0_r8
      theta_timeavg_old=0.0_r8
      DO iter=1,total_iters
! 
! Calculate wave related bed roughness from VA2013 A.5     
! 
        ksw=ksw_calc(d50, mu_calc(d50), theta_timeavg, eta, rlen)
!
! Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4 
!
        fw=fw_calc(ahat, ksw)
!
! Calculate current-related bed roughness from VA2013 Appendix A.1
!
        ksd=ksd_calc(d50, d90, mu_calc(d50), theta_timeavg,             &
     &                                                  eta, rlen)    
!
! Calculate full-cycle current friction factor from VA2013 Eqn. 20
! 
!      Within COAWST use bustr, bvstr to get delta 
!    % In the standalone matlab code, delta is an input 
!    %   
!    % Use fd (full cycle current friction to get udelta)
!    % udelta=current velocity at the top of the wave boundary layer
!
        fd=fd_calc(umag_curr, Zref, ksd)
        ustarc=SQRT(0.5_r8*fd)*umag_curr;  !    %friction velocity [m/s]
        udelta=ustarc/(0.4_r8*log(30.0_r8*delta/ksd));
!
! Calculate Time-averaged absolute Shields stress VA2013 Appendix Eq. A.3
! 
        theta_timeavg=osmgd*(0.5_r8*fd*udelta**2.0_r8+                  &
     &                       0.25_r8*fw*uhat**2.0_r8)
!
        IF(ABS(theta_timeavg-theta_timeavg_old).lt.tol) THEN 
          EXIT 
        ENDIF
        theta_timeavg_old=theta_timeavg 
      END DO 
!
! Calculate wave Reynolds stress from full cycle wave and friction factor
! that were formed from the iterative cycle, VA2013, Eqn.22
! 
      alpha=udelta/(udelta+uhat);
      fwd=alpha*fd+(1.0-alpha)*fw;
!
      k_wn=kh(Td,depth)/depth     ! Wave number 
      c_w=2*pi/(k_wn*Td)          ! Wave speed
      alpha_w=0.424_r8
!
      tau_wRe=rho0*fwd*alpha_w*uhat**3.0_r8/(2.0_r8*c_w)
!        
! Compute the change in time period based on converged udelta 
! (current velocity at wave boundary layer)
!
      CALL current_timeperiod(udelta, phi_curwave, umax, umin, RR,      &
     &                        T_c, T_t, Td)
!
! Calculate the effect of surface waves 
!
!     CALL surface_wave_mod(Td, depth, uhat, T_c, T_cu, T_t, T_tu)
! 
      END SUBROUTINE full_wave_cycle_stress_factors
!
      SUBROUTINE half_wave_cycle_stress_factors( T_iu, T_i, ahat, ksw,  &
     &                                           fd, alpha,             &
     &                                           d50, osmgd,            &
     &                                ui_r, uhat_i, udelta, phi_curwave,&
     &                                           tau_wRe,               &
     &                            dsf, theta_ix, theta_iy, mag_theta_i )
!
!**********************************************************************
!  This subroutine returns the following:
!  dsf                 : sheetflow thickness
!  theta_ix, theta_iy  : Shields parameter in x and y dir. 
!  mag_theta_i         : Magnitude of Shields parameter for half cycle
!**********************************************************************
! 
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!  
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles 
!       
      real(r8), intent(in) :: T_iu, T_i, ahat, ksw
      real(r8), intent(in) :: fd, alpha
      real(r8), intent(in) :: d50, osmgd
      real(r8), intent(in) :: ui_r, uhat_i, udelta, phi_curwave
      real(r8), intent(in) :: tau_wRe
      real(r8), intent(inout) :: dsf, theta_ix, theta_iy, mag_theta_i
!
!  Local variables
! 
!      real(r8) :: psi
!      real(r8) :: lambda  ! ripple dimensions
!      real(r8) :: ksd, ksd_calc
!      real(r8) :: ksw, ksw_calc
!      real(r8) :: fw, fw_calc
!      real(r8) :: fd_calc
      real(r8) :: fw_i, fwi_calc, fwd_i 
      real(r8) :: dsf_calc
      real(r8) :: alpha_w, fwd, kh, k, c_w
      real(r8) :: theta_hat_i
      real(r8) :: ui_rx, ui_ry, mag_ui
!        
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
! 
      fw_i=fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
! Different for crest and trough 
!
      fwd_i=alpha*fd+(1.0_r8-alpha)*fw_i
!
! VA2013-Magnitude of Shields parameter Eqn. 17
! 
      theta_hat_i=0.5_r8*fwd_i*uhat_i**2*osmgd
!
! Sheet flow thickness VA2013 Appendix C C.1 
! Update from initial value 
!
      dsf=dsf_calc(d50, theta_hat_i) !this dsf is in m 
!
! Calculated the velocity magnitude based on representative velocities
!  equation 12 from Van der A, 2013 
!
!-----------------------------------------------------------------------
! Get the representative trough half cycle water particle velocity
!    as well as full cycle orbital velocity and excursion
!-----------------------------------------------------------------------
!
      ui_rx=ui_r+udelta*COS(phi_curwave)
      ui_ry=udelta*SIN(phi_curwave)
      mag_ui=SQRT(ui_rx*ui_rx+ui_ry*ui_ry)
!
! VA2013-Magnitude of Shields parameter Eqn. 17
! 
      mag_theta_i=0.5_r8*fwd_i*mag_ui**2*osmgd
!
!-----------------------------------------------------------------------
! Shields parameter in crest cycle
! rho0 is required for non-dimensionalizing 
!-----------------------------------------------------------------------
!
      theta_ix=ABS(mag_theta_i)*ui_rx/(mag_ui)+tau_wRe*osmgd/rho0
      theta_iy=ABS(mag_theta_i)*ui_ry/(mag_ui)
      mag_theta_i=sqrt(theta_ix*theta_ix+theta_iy*theta_iy)
! 
      END SUBROUTINE half_wave_cycle_stress_factors
!
      SUBROUTINE ripple_dim(psi, d50, eta, rlen)
! 
!**********************************************************************
!  This subroutine returns the following:
!  eta, rlen  : Ripple dimensions: (height and length) 
!**********************************************************************
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
!        
      USE mod_kinds
      USE mod_scalars 

      implicit none 
!
      real(r8), intent(in)  :: psi, d50
      real(r8), intent(out) :: eta, rlen
!
      real(r8) :: d50_mm 
      real(r8) :: m_eta, m_lambda, n_eta, n_lambda 
!     
      d50_mm=0.001_r8*d50
      IF(d50_mm.lt.0.22_r8) THEN
        m_eta=0.55_r8
        m_lambda=0.73_r8
      ELSEIF(d50_mm.ge.0.22_r8.and.d50_mm.lt.0.30_r8) THEN
        m_eta=0.55_r8+(0.45_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
        m_lambda=0.73_r8+(0.27_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
      ELSE
        m_eta=1.0_r8
        m_lambda=1.0_r8
      ENDIF
! 
! Smooth transition between ripple regime and bed sheet flow regime 
!
      IF(psi.le.190.0_r8) THEN
        n_eta=1.0_r8
      ELSEIF(psi.gt.190.0_r8.and.psi.lt.240.0_r8) THEN
        n_eta=0.5_r8*(1.0_r8+cos(pi*(psi-190.0_r8)/(50.0_r8)))
      ELSEIF(psi.ge.240.0_r8) THEN
        n_eta=0.0_r8
      ENDIF
      n_lambda=n_eta
!
      eta=MAX(0.0_r8,m_eta*n_eta*(0.275_r8-0.022*psi**0.42_r8))
      rlen=MAX(0.0_r8,m_lambda*n_lambda*                                &
     &                             (1.97_r8-0.44_r8*psi**0.21_r8))
!
      RETURN
      END SUBROUTINE ripple_dim
!
      SUBROUTINE current_timeperiod(unet, phi_curwave, umax, umin,      &
     &                              RR, T_c, T_t, Td)
!
!**********************************************************************
!  This subroutine returns the following:
!  T_c, T_t  : Time period in crest and trough cycle
!**********************************************************************
!
! Modify the crest and trough time periods based on current velocity
! This function was developed by Chris Sherwood and Tarandeep Kalra
!
! The basis for these formulations are formed from Appendix A.3 
! in SANTOSS report. 
! Report number: SANTOSS_UT_IR3 Date: January 2010
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!  
      real(r8), intent(in) :: unet, phi_curwave
      real(r8), intent(in) :: umax, umin
      real(r8), intent(in) :: RR, Td
      real(r8), intent(inout) :: T_c, T_t
!
!  Local variables
!
      real(r8) :: unet_xdir, udelta, delt
!
      unet_xdir=unet*cos(phi_curwave/180*pi);

      IF(RR.eq.0.5_r8) THEN
        T_c=0.5_r8*Td
        T_t=0.5_r8*Td
        IF(unet_xdir.ge.umax) THEN
          T_c=Td
          T_t=0.0_r8
        ELSEIF(unet_xdir.le.umin) THEN 
          T_c=0.0_r8
          T_t=Td 
        ELSEIF(unet_xdir.lt.0.0_r8.and.unet_xdir.gt.umin) THEN
          delt=ASIN(-unet/umin)/pi
          T_t=T_t*(1.0_r8-2.0_r8*delt)
          T_c=Td-T_t;
        ELSEIF(unet_xdir.gt.0.0_r8.and.unet_xdir.lt.umax) THEN 
          delt=ASIN(unet_xdir/(-umax))/pi
          T_c=T_c*(1.0_r8-2.0_r8*delt)
          T_t=Td-T_c
        ELSEIF(unet_xdir.eq.0.0_r8) THEN 
          T_c=T_c
          T_t=T_t 
        ENDIF 
      ELSEIF(RR.gt.0.5_r8) THEN 
        T_c=T_c
        T_t=T_t 
        IF(unet_xdir.ge.umax) THEN 
          T_c=Td
          T_t=0.0_r8
        ELSEIF(unet_xdir<=umin) THEN 
          T_c=0.0_r8
          T_t=Td
        ELSEIF(unet_xdir.lt.0.0_r8.and.unet_xdir.gt.umin) THEN 
          delt=ASIN(-unet_xdir/umin)/pi
          T_t=T_t*(1.0_r8-2.0_r8*delt)
          T_c=Td-T_t
        ELSEIF(unet_xdir.gt.0.0_r8.and.unet_xdir.lt.umax) THEN 
          delt=ASIN(unet_xdir/(-umax))/pi 
          T_c=T_c*(1.0_r8-2.0_r8*delt) 
          T_t=Td-T_c
        ELSEIF(unet_xdir.eq.0.0_r8) THEN 
          T_c=T_c 
          T_t=T_t
        ENDIF
      ENDIF   
!
      END SUBROUTINE current_timeperiod
!
      SUBROUTINE surface_wave_mod(Td, depth, uhat,                      &
     &                            T_c, T_cu, T_t, T_tu)
! 
!**********************************************************************
!  This subroutine returns the following:
!  T_c, T_cu, T_t, T_tu  : Change in time period in crest and 
!                          trough cycle due to particle displacement
!                          under surface waves. 
!**********************************************************************
!
! Crest period extension for horizontal particle displacement.
! Tuning factor eps_eff = 0.55 from measurements GWK Schretlen 2010         
! Equations in Appendix B of SANTOSS Report 
! Report number: SANTOSS_UT_IR3 Date: January 2010 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
! 
      real(r8), intent(in) :: Td, depth, uhat 
      real(r8), intent(inout) :: T_c, T_cu, T_t, T_tu
!
!  Local variables
!
      real(r8) :: k_wn, kh, eps_eff, c
      real(r8) :: delta_Tc, delta_Tt
      real(r8) :: T_c_new, T_cu_new
      real(r8) :: T_t_new, T_tu_new
!
      k_wn=kh(Td,depth)/depth
      c=2.0_r8*pi/(k_wn*Td) 
!
      eps_eff=0.55
!
      delta_Tc=eps_eff*uhat/(c*pi-eps_eff*2.0*uhat)
      T_c_new=T_c+delta_Tc
      T_cu_new=T_cu*T_c_new/T_c
!
      delta_Tt=eps_eff*uhat/(c*pi+eps_eff*2.0*uhat)
      T_t_new=T_t-delta_Tt
      T_tu_new=T_tu*T_t_new/T_t
!
      T_c=T_c_new
      T_cu=T_cu_new
      T_t=T_t_new
      T_tu=T_tu_new
!
      END SUBROUTINE surface_wave_mod
!
      REAL(r8) FUNCTION theta_cr_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
!
      USE mod_kinds
      USE mod_scalars

      implicit none 
! 
      real(r8) :: d50, rhos
      real(r8) :: s, dstar 
      real(r8) :: cff1, cff2 
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff1=0.30_r8/(1.0_r8+1.2_r8*dstar)
      cff2=0.055_r8*(1.0_r8-EXP(-0.020_r8*dstar))
      theta_cr_calc=cff1+cff2
!
      RETURN
      END FUNCTION theta_cr_calc
!
      REAL(r8) FUNCTION w_s_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
! Dynamics of Marine Sands 
!
      USE mod_kinds
      USE mod_scalars

      implicit none 
! 
      real(r8) :: d50, rhos
      real(r8) :: s, dstar 
      real(r8) :: cff, cff1
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff=nu/d50
      cff1=10.36_r8
      w_s_calc=cff*(sqrt(cff1*cff1+1.049_r8*dstar**3.0_r8)-cff1)
!
      RETURN
      END FUNCTION 
!
      REAL(r8) FUNCTION mu_calc(d50)
!
! Calculate bed roughness factor based on grain size
! VA2013 Appendix A., required for current related bed roughness
! and wave related bed roughness. 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!
      real(r8) :: d50, d50_mm
!
      d50_mm=d50*1000.0_r8 ! converting from m to mm 
!
      IF(d50_mm.le.0.15_r8) THEN
        mu_calc=6.0_r8
      ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8) THEN
        mu_calc=6.0_r8-5.0_r8*((d50_mm-0.15_r8)/(0.2_r8-0.15_r8))
      ELSEIF(d50_mm.gt.0.20_r8) THEN
        mu_calc=1.0_r8
      ENDIF
!
      RETURN
      END FUNCTION mu_calc
!
      REAL(r8) FUNCTION ksd_calc(d50, d90, mu, theta_timeavg,           &
     &                           eta, rlen)
!
! Calculate current-related bed roughness from VA2013 Appendix A.1.
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
      real(r8) :: d50, d90, mu, theta_timeavg, eta, rlen
      real(r8) :: ripple_fac
!
!      eta=MAX(eta,d50)
      rlen=MAX(rlen,d50)
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksd_calc=MAX( 3.0_r8*d90,                                         &
     &             d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )+            &
     &             ripple_fac
!
      RETURN 
      END FUNCTION ksd_calc
!
      REAL(r8) FUNCTION ksw_calc(d50, mu, theta_timeavg, eta, rlen)
!
! Calculate wave related bed roughness from VA2013 Eqn. A.5.
! 
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
      real(r8) :: d50, mu, theta_timeavg, eta, rlen
      real(r8) :: ripple_fac, ksw
!
!      eta=MAX(eta,d50)
      rlen=MAX(rlen,d50)
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksw_calc=MAX( d50,                                                &
     &                d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )          &
     &               +ripple_fac
!
      RETURN
      END FUNCTION ksw_calc
!
      REAL(r8) FUNCTION fw_calc(ahat, ksw)
!
! Calculate full-cycle wave friction factor from VA2013 Eqn. A.4. 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
      real(r8) :: ahat, ksw, ratio, fw
!
      ratio=ahat/ksw
      IF(ratio.gt.1.587_r8) THEN
        fw_calc=0.00251_r8*EXP(5.21_r8*(ratio)**(-0.19_r8))
      ELSE
        fw_calc=0.3_r8
      ENDIF
!
      RETURN
      END FUNCTION fw_calc 
!
      REAL(r8) FUNCTION fd_calc(umag_curr, dsf, ksd)
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!
! Calculate current related friction factor VA2013 Eqn. 20
! Assuming logarithmic velocity profile.
! 
      real(r8) :: umag_curr, dsf, ksd
      real(r8), parameter :: von_k=0.41_r8
! 
      IF(umag_curr.eq.0.0_r8) THEN
        fd_calc=0.0_r8
      ELSE
        fd_calc=2.0_r8*(von_k/LOG(30.0_r8*dsf/ksd))**2.0_r8
      ENDIF 
!     
      RETURN
      END FUNCTION fd_calc 
!
      REAL(r8) FUNCTION fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21.
!       
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
      real(r8) :: T_iu, T_i, ahat, ksw
      real(r8) :: c1, ratio, fwi
      real(r8) :: cff 
!
      c1=2.6_r8
      ratio=ahat/ksw
      IF(ratio.gt.1.587_r8) THEN
        cff=(2.0_r8*T_iu/T_i)**c1
        fwi_calc=0.00251_r8*EXP(5.21_r8*(cff*ratio)**(-0.19_r8))
      ELSE
        fwi_calc=0.3_r8
      ENDIF
!     
      RETURN
      END FUNCTION fwi_calc
!
      REAL(r8) FUNCTION dsf_calc(d50, theta_i)
!
! Sheet flow thickness VA2013 Appendix C.1.
!       
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
      real(r8) :: d50, theta_i
      real(r8) :: d50_mm
      real(r8) :: cff
!
      d50_mm=d50*1000.0_r8 !converting from m to mm 
      IF(d50_mm.le.0.15_r8)THEN
        cff=25.0_r8*theta_i
      ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8)THEN 
        cff=25.0_r8-(12.0_r8*(d50_mm-0.15_r8)/0.05_r8)
      ELSEIF(d50_mm.ge.0.20_r8)THEN 
        cff=13.0_r8*theta_i
      ENDIF 
      dsf_calc=MAX(d50*cff,d50)
!
      RETURN
      END FUNCTION dsf_calc
!
!  End of functions for step 2 for Shear stress formulation
!
      SUBROUTINE skewness_params(iunit, H_s, T, depth, r, phi)
!        
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "bb" and "r"
! Given input of H_s, T and depth 
! r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
! phi   - skewness/asymmetry parameter                         [value]
! Su     - umax/(umax-umin)                                    [value]
! Au   - amax/(amax-amin)                                      [value]
! alpha - tmax/pi                                              [value]
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!
      integer, intent(in) :: iunit
      real(r8), intent(in)  :: H_s, T, depth
      real(r8), intent(out) :: r, phi
!
! Local variables 
! 
      real(r8), parameter :: p1=0.0_r8
      real(r8), parameter :: p2=0.857_r8
      real(r8), parameter :: p3=-0.471_r8
      real(r8), parameter :: p4=0.297_r8
      real(r8), parameter :: p5=0.815_r8
      real(r8), parameter :: p6=0.672_r8
      real(r8) :: a_w, Ur
      real(r8) :: B, psi, bb 
      real(r8) :: k_wn, kh, cff
      real(r8) :: Su, Au
!
! Ruessink et al., 2012, Coastal Engineering 65:56-63.
!
! k is the local wave number computed with linear wave theory.
!
      k_wn=kh(T,depth)/depth       
!
!      H_s=sqrt(2.0_r8)*H_rms
      a_w=0.5_r8*H_s 
      Ur=0.75_r8*a_w*k_wn/((k_wn*depth)**3.0_r8)
!
! Ruessink et al., 2012 Equation 9.
!
      cff=EXP( (p3-log10(Ur)) /p4)
      B=p1+((p2-p1)/(1.0_r8+cff))
      psi=-90.0_r8*deg2rad*(1.0_r8-TANH(p5/Ur**p6))
! 
! Markaley and Davies, equation provides bb which is "b" in paper
! Check from where CRS found these equations
! 
      bb=sqrt(2.0_r8)*B/(sqrt(2.0_r8*B**2.0_r8+9.0_r8))
      r=2.0_r8*bb/(bb**2.0_r8+1.0_r8)
!
! Ruessink et al., 2012 under Equation 12.
!
      phi=-psi-0.5_r8*pi
!
! Where are these asymmetry Su, Au utilized 
! recreate the asymetry 
!          
      Su=B*cos(psi)
      Au=B*sin(psi)

      write(iunit,*)"****************************************"
      write(iunit,*)" Reuissink parameters"
      write(iunit,*) "b=", bb
      write(iunit,*) "r=",r
      write(iunit,*) "phi=",phi
      write(iunit,*)"****************************************"
      RETURN     
      END SUBROUTINE skewness_params

        SUBROUTINE abreu_points(iunit, r, phi, Uw, T, DTc, DTt,         &
     &                                 DTcu, DTtu, umax, umin, RR)
! 
!  Calculate umax, umin, and phases of asymmetrical wave orbital velocity 
!
!  Use the asymmetry parameters from Ruessink et al, 2012
!  to get the umax, umin and phases of asymettrical wave 
!  orbital velocity to be used by Van Der A. 
!  T_c is duration of crest
!  T_cu Duration of accerating flow within crest half cycle
!
        USE mod_kinds
        USE mod_scalars
!
        implicit none 
!
        integer, intent(in) :: iunit
        real(r8), intent(in)  :: r, phi, Uw, T
        real(r8), intent(out) :: Dtc, DTt, DTcu, DTtu 
        real(r8), intent(out) :: umax, umin, RR
!
! Local variables 
! 
        real(r8) :: b, c, ratio, tmt, tmc, tzd, tzu
        real(r8) :: omega, w, phi_new 
        real(r8) :: P, F0, betar_0, beta
        real(r8) :: T_tu, T_cu, T_c, T_t 
        real(r8) :: cff1, cff2, cff 
        real(r8) :: Sk, Ak
!
        omega=2.0_r8*pi/T
!
        phi_new=-phi

! Malarkey and Davies (Under equation 16b) 
        P=SQRT(1.0_r8-r*r)
!
! Malarkey and Davies (Under equation 16b) 
!
        b=r/(1.0_r8+P)
!
! Appendix E of Malarkey and Davies 
!
        c=b*SIN(phi_new)
!
        cff1=4.0_r8*c*(b*b-c*c)+(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
        cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
        ratio=cff1/cff2
!
! These if conditionals prevent ASIN to be between [-1,1] and prevent NaNs
! Not a problem in the MATLAB code
!
        IF(ratio.gt.1.0_r8)THEN 
          ratio=1.0_r8
        ENDIF
        IF(ratio.lt.-1.0_r8)THEN
          ratio=-1.0_r8
        ENDIF 
        tmc=ASIN(ratio)
!
        cff1=4.0_r8*c*(b*b-c*c)-(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
        cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
        ratio=cff1/cff2
        IF(ratio.gt.1.0_r8)THEN 
          ratio=1.0_r8
        ENDIF
        IF(ratio.lt.-1.0_r8)THEN
          ratio=-1.0_r8
        ENDIF 
        tmt=ASIN(ratio)
!         
        IF(tmc.lt.0.0_r8) THEN 
          tmc=tmc+2.0_r8*pi
        ENDIF
        IF(tmt.lt.0.0_r8) THEN 
          tmt=tmt+2.0_r8*pi
        ENDIF 
! 
! Non dimensional umax and umin, under E5 in Malarkey and Davies 
! 
        umax=1.0_r8+c
        umin=umax-2.0_r8
!
!       Dimensionalize
!
        umax=umax*Uw
        umin=umin*Uw
!
! phase of zero upcrossing and downcrossing (radians)
!
        tzu=ASIN(b*SIN(phi_new))
        tzd=2.0_r8*ACOS(c)+tzu 
! 
! MD, equation 17
!
        RR=0.5_r8*(1.0_r8+b*SIN(phi_new)) 
! 
! MD, under equation 18
! 
        IF(r.le.0.5_r8) THEN
          F0=1.0_r8-0.27_r8*(2.0_r8*r)**(2.1_r8)
        ELSE
          F0=0.59_r8+0.14_r8*(2.0_r8*r)**(-6.2_r8)
        ENDIF
!
! MD, Equation 15a,b 
!
        IF(r.ge.0.0_r8.and.r.lt.0.5)THEN
          betar_0=0.5_r8*(1.0_r8+r)
        ELSEIF(r.gt.0.5_r8.and.r.lt.1.0_r8)THEN
          cff1=4.0_r8*r*(1.0_r8+r)
          cff2=cff1+1.0_r8
          betar_0=cff1/cff2
        ENDIF
!
! MD, Equation 18
!
        cff=SIN((0.5_r8*pi-ABS(phi_new))*F0)/SIN(0.5_r8*pi*F0)
        beta=0.5_r8+(betar_0-0.5_r8)*cff
!
! MD, Table 1, get asymmetry parameterization
! using GSSO (10a,b)
!
        cff=SQRT(2.0_r8*(1.0_r8+b*b)**3.0_r8)
        Sk=3.0_r8*SIN(phi_new)/cff
        Ak=-3.0_r8*COS(phi_new)/cff
!
! These are the dimensional fractions of wave periods needed by Van der A eqn.
! TSK - Check source of these equations 
!
        w=1.0_r8/omega
        DTc=(tzd-tzu)*w
        DTt=T-DTc;
        DTcu=(tmc-tzu)*w
        DTtu=(tmt-tzd)*w
!
        T_tu=DTtu
        T_cu=DTcu
        T_c=DTc
        T_t=DTt
!
        write(iunit,*)"****************************************"
        write(iunit,*)" Abreu outputs"
        write(iunit,*) "Beta=", beta, "R=",RR
        write(iunit,*)"Crest cycle 'T_c='", T_c, "T_cu=", T_cu
        write(iunit,*)"Trough cycle 'T_t='", T_t, "T_tu=", T_tu
        write(iunit,*)"****************************************"
        
        RETURN
        END SUBROUTINE abreu_points

        REAL(r8) FUNCTION kh(Td,depth)
!
!  Calculate wave number from Wave period and depth 

! RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
! HR Wallingford Report TR 155, February 2006
!
        USE mod_kinds
        USE mod_scalars
!
        implicit none 
! 
        real(r8) :: Td, depth
        real(r8) :: cff
        real(r8) :: x, y, t, omega 
!
        omega=2.0_r8*pi/Td
        x=omega**2.0_r8*depth/g
! 
        IF(x.lt.1.0_r8) THEN
          y=SQRT(x)
        ELSE
          y=x
        ENDIF
!
! Iteratively solving 3 times for eqn.7 of Soulsby 1997 by using 
! eqns. (12a-14)
!      
        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff
!
        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff

        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff
        kh=y
!       
        RETURN 
        END FUNCTION 
!
