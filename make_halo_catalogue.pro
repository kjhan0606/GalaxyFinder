pro make_halo_catalogue,nout=nout
	
	rd_info,info,nout=nout
	snout0=string(nout,format='(i03)')
	snout1=string(nout,format='(i05)')
	file0='FoF_Data/FoF.'+snout1+'/GALCATALOG.LIST.'+snout1
	file1='FoF_Data/FoF.'+snout1+'/background_ptl.'+snout1
	output='FoF_Data/halo_'+snout1+'.sav'
	sw_background=1
	thalo={nsub:0L, ndm:0L, nstar:0L, nsink: 0L, ngas: 0L, npall: 0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tsub={ndm:0L, ngas:0L, nsink: 0L, nstar: 0L, npall: 0L, dum:0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tpart={pos:dblarr(3),vel:dblarr(3),mass:0.d,mass0:0.d,tp:0.d,zp:0.d,ID:0L, family:0b, tag:0b, dum:0,level:0L,dum1:0L}
	tsink={pos:dblarr(3),vel:dblarr(3),mass:0.d,tbirth:0.d,angm:dblarr(3),ang:dblarr(3),dmsmbh:dblarr(3),esave:0.d,smag:0.d,eps:0.d,id:0L,dum0:0L}
	tgas={pos:dblarr(3),dx:0.d,vel:fltarr(3),dum0:0.,density:0.d,temp:0.,metal:0.,mass:0.,dum1:0.,potential:0.d,f:dblarr(3)}

	bpart=n_tags(tpart,/length)
	bsink=n_tags(tsink,/length)
	bgas=n_tags(tgas,/length)

	nhalo=0L
	nsub=0L	
	openr,10,file0
	while ~eof(10) do begin
		nhalo++
		readu,10,thalo
		for i=0L,thalo.nsub-1 do begin
			readu,10,tsub
			nsub++
		endfor
	endwhile
	close,10
	halo=replicate(thalo,nhalo)
	subhalo=replicate(tsub,nsub)
	sub_host=lonarr(nsub)
	unbound=replicate(tsub,nhalo)
	openr,10,file0
	k=0L
	for i=0L,nhalo-1 do begin
		readu,10,thalo
		halo[i]=thalo
		for j=0L,thalo.nsub-1 do begin
			readu,10,tsub
			subhalo[k]=tsub
			sub_host[k]=i
			k++
		endfor
	endfor
	print,'Making halo catalogue'
	print,'n_halo=',nhalo
	close,10
	openr,10,file1
	for i=0L,nhalo-1 do begin
;		nhalo=nhalo+1
		readu,10,tsub
		unbound[i]=tsub
		sbdm=bpart*tsub.ndm
		sbstar=bpart*tsub.nstar
		sbgas=bgas*tsub.ngas
		sbsink=bsink*tsub.nsink
		skip_lun,10,sbdm
		skip_lun,10,sbgas
		skip_lun,10,sbsink
		skip_lun,10,sbstar
	
	endfor
	close,10
	save,file=output,halo,subhalo,unbound,info,sub_host



end
