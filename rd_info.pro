;+
; NAME:
;	RD_INFO
;
; purpose:
;	
;
; CATEGORY:
;	
;
; CALLING SEQUENCE:
;       RD_INFO,info_file,simu_info,verbose=verbose
;
; INPUTS:
;	info_file: ramses info_xxxxx.txt information output
;       
; OUTPUTS:
;       simu_info: structure containing all scale and snapshot info
;            ;; code units ;;
;            .unit_l: box size in cgs
;            .unit_d: density convertion in cgs
;            .unit_v: velocity conversion in cgs 
;            .unit_t: time conversion in cgs 
;            .unit_T2: Temperature convertion in cgs
;            .boxlen: box size in code unit
;            .levmin: min level of refinment
;            .levmax: max level of refinment
;            ;; useful convertions ;;
;            .boxtompc: box size in Mpc
;            .tGyr: time in Gyr
;
; OPTIONAL INPUTS:
;      VERBOSE: if set verbose mode
;   
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;
;	       RD_INFO,"output_00381/info_00381.txt",simu_info,/verbose
;
; MODIFICATION HISTORY:
; 	Written by:	Dylan Tweed, 2010
;-


pro rd_info,info,nout=nout,file=file_info,verbose=verbose,silent=silent,ierr=ierr,path=path
 

  curr_dir='./'
  if keyword_set(path) then begin
     cd, path, curr=curr_dir
  endif 
  if keyword_set(nout) then begin
     snout=string(nout, format='(I05)')
     file_info='output_'+snout+'/info_'+snout+'.txt'
     if ~file_test(file_info) then file_info='../'+file_info
     if ~file_test(file_info) then file_info='Info.txt/info_'+snout+'.txt'
     if ~file_test(file_info) then file_info='../Info.txt/info_'+snout+'.txt'
  endif
  if keyword_set(silent) then verbose=0

  ierr=0
  if ~file_test(file_info) then begin
    ierr=1
    goto,program_end
  endif
 
  line          = "data je sais pas"
  dataname      = string("",format='(a13)')
  my_narr       = lonarr(6)
  my_narr(0:5)  = 0L
  my_rarr       = dblarr(11)
  my_rarr(0:10) = 0.


  if(file_test(file_info) ne 1) then begin
     print,"file ",file_info," does not exist"
     stop
  endif
  openr,2,file_info
  if keyword_set(verbose) then print, "reading : ",file_info
  value = 0L
  for il = 0,5 do begin
     readf,2,dataname,value,format='(a13,I11)'
     my_narr(il) = value
  end
  data="alors je sais vraiment pas qupi mettre pour que ca marche"
  readf,2,data,format='(a50)'
  value =0d0
  for il = 0,10 do begin
     readf,2,dataname,value,format='(a13,E23.15)'
     my_rarr(il) = value
  end
  close,2
 
  time=my_rarr(1)
 
  if keyword_set(verbose) then begin
     print,"ncpu          :",my_narr(0)
     print,"ndim          :",my_narr(1)
     print,"levmin,levmax :",my_narr(2),my_narr(3)
     print,"ngridmax      :",my_narr(4)
     print,"nstep_coarse  :",my_narr(5)
     print,"boxlen        :",my_rarr(0)
     print,"time          :",my_rarr(1)
     print,"aexp          :",my_rarr(2)
     print,"H0            :",my_rarr(3)
     print,"omega_m       :",my_rarr(4)
     print,"omega_l       :",my_rarr(5)
     print,"omega_k       :",my_rarr(6)
     print,"omega_b       :",my_rarr(7)
     ;; print,"unit_l        :",my_rarr(8)
     ;; print,"unit_d        :",my_rarr(9)
     ;; print,"unit_t        :",my_rarr(10)
  endif
  
;; unit in cgs
  kpc     = 3.086D21
  twopi   = 6.2831853d0
  hplanck = 6.6262000d-27
  eV      = 1.6022000d-12
  kB      = 1.38D-16
  clight  = 2.9979250d+10
  Gyr     = 3.1536000d+16
  X       = 0.76
  Y       = 0.24 
  rhoc    = 1.8800000d-29
  mH      = 1.6600000d-24
  mu_mol  = 1.2195d0
  G       = 6.67D-8
  m_sun   = 1.98892e33
  
  if keyword_set(verbose) then begin
     print,"time Gyr      :",my_rarr(1)*my_rarr(10)/Gyr
     print,"boxsize (kpc) :",my_rarr(0)*my_rarr(8)/kpc
     print,"res     (kpc) :",my_rarr(0)*my_rarr(8)*2^(-1.*my_narr(3))/kpc
     print,"res (boxunit) ;",my_rarr(0)*2^(-1.*my_narr(3))
     print,"Lf      (Mpc) ;",my_rarr(0)*my_rarr(8)/(1e3*kpc)*1/my_rarr(2)
     print,"redshift      :",1./my_rarr(2) - 1
  endif
  
  scale_l    = my_rarr(8)
  scale_d    = my_rarr(9)
  scale_t    = my_rarr(10)
;; scale_v convert velocity in user units into cm/s
  scale_v    = scale_l / scale_t
;; scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2   = mH/kB * scale_v^2.
  scale_nH   = X/mH * scale_d
;; scale covert mettalicty into solar metallicity Z_sol = 0.02
  scale_Z    = 1./0.02 
  scale_flux = scale_v*scale_d*kpc*kpc*Gyr/m_sun
;;  print,"scale_flux",scale_flux
  if(keyword_set(verbose)) then begin
     print,'unit_l :',scale_l
     print,'unit_d :',scale_d
     print,'unit_t :',scale_t
     print,'unit_v :',scale_v
     print,'unit_T2:',scale_T2
     print,'unit_nH:',scale_nH
  end
  
;;jl = sqrt(3.*!pi*kB /(32.*G)) / mH did not agree first hand
;; used  jl = sqrt(!pi*kB /(G)) / mH
  
  info={boxtokpc:my_rarr(0)*scale_l/kpc,$
        tGyr:my_rarr(1)*scale_t/Gyr,$
        boxlen:my_rarr(0),$
        ndim:my_narr(1),$
        lmin:my_narr(2),$
        lmax:my_narr(3),$
        unit_l:scale_l,$
        unit_d:scale_d,$
        unit_t:scale_t,$
        unit_v:scale_v,$
        unit_nH:scale_nH,$
        unit_T2:scale_T2,$
        unit_Z:scale_Z,$
        unit_B:sqrt(4*!dpi*scale_d*scale_v^2d0),$
        kms:scale_v/1d5,$
        unit_flux:scale_d*scale_v*(1e-9*Gyr)*(kpc)*(kpc)/m_sun,$
        punit_m:scale_d*scale_l^3d0/2d33, $
        pboxsize:my_rarr(0)*scale_l/(kpc*1000),$
        boxkpc:my_rarr(0)*scale_l/kpc,$
        kpc:my_rarr(0)*scale_l/kpc,$
        pc:my_rarr(0)*scale_l/kpc*1d3,$
        mpc:my_rarr(0)*scale_l/kpc/1d3,$
        time:time,$
        aexp:my_rarr(2),$
        zred:1/my_rarr(2)-1,$
        H0:my_rarr(3),$
        om:my_rarr(4),$
        ol:my_rarr(5),$
        ok:my_rarr(6),$
        ob:my_rarr(7),$
        msun:scale_d*scale_l^3d0/1.989d33,$
        yr:scale_t/365./24./3600}
  
  if keyword_set(verbose) then begin
     print,"zred    :",info.zred
     print,"time Gyr:",info.tGyr
     print,"box kpc :",info.boxtokpc
     print,"scale t :",info.unit_t
     ;; print,"scale l:",info.unit_l
     ;; print,"scale d:",info.unit_d
     ;; print,"scale nH:",info.unit_nH
  endif   

  program_end:
  cd, curr_dir
end
