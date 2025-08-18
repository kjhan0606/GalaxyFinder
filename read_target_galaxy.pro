pro read_target_galaxy,nout=nout,list=list,dump_file=dump_file,data=data,info=info


	

	out_path='catalogue/'
	snout1=string(nout,format='(i05)')

		

	info=1.
	
	snout0=string(nout,format='(i03)')

	path='catalogue/sn'+snout0+'/halo_data'
	ft=file_test(path)
	if(ft[0] eq 0) then spawn,'mkdir -p '+path

	print,'Extracting data from snapshot '+snout0
	file0='FoF_Data/FoF.'+snout1+'/GALFIND.DATA.'+snout1
	file1='FoF_Data/FoF.'+snout1+'/GALCATALOG.LIST.'+snout1
	file2='FoF_Data/halo_'+snout1+'.sav'
	tp=file_test(file2)
	if(tp[0] eq 0) then begin
		 make_halo_catalouge,nout=nout	
	endif

	restore,file2

	thalo={nsub:0L, ndm:0L, nstar:0L, nsink: 0L, ngas: 0L, npall: 0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tsub={ndm:0L, ngas:0L, nsink: 0L, nstar: 0L, npall: 0L, dum:0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tpart={pos:dblarr(3),vel:dblarr(3),mass:0.d,mass0:0.d,tp:0.d,zp:0.d,ID:0L, family:0b, tag:0b, dum:0,level:0L,dum1:0L}
	tsink={pos:dblarr(3),vel:dblarr(3),mass:0.d,tbirth:0.d,angm:dblarr(3),ang:dblarr(3),dmsmbh:dblarr(3),esave:0.d,smag:0.d,eps:0.d,id:0L,dum0:0L}
	tgas={pos:dblarr(3),dx:0.d,vel:fltarr(3),dum0:0.,density:0.d,temp:0.,metal:0.,mass:0.,dum1:0.,potential:0.d,f:dblarr(3)}
	bhalo=n_tags(thalo,/length)
	bsub=n_tags(tsub,/length)
	bpart=n_tags(tpart,/length)
	bsink=n_tags(tsink,/length)
	bgas=n_tags(tgas,/length)
	

	nhalo=0L
	nsub=0L
	tnstar=total(halo.nstar,/integer)
	nhalo=n_elements(halo)	
	nsub=n_elements(subhalo)	

	openr,10,file0



	print,'Nhalo total: ',nhalo
	print,'Nsub total: ',nsub



	idx=where(list ge 0)
	ctnstar=total(subhalo[0:nsub-2].nstar,/integer,/cumulative)
	ctndm=total(subhalo[0:nsub-2].ndm,/integer,/cumulative)
	ctngas=total(subhalo[0:nsub-2].ngas,/integer,/cumulative)
	ctnsink=total(subhalo[0:nsub-2].nsink,/integer,/cumulative)
;	ctnsub=total(halo[0:nhalo-2].nsub,/integer,/cumulative)

	ctnstar=[0,ctnstar]
	ctndm=[0,ctndm]
	ctngas=[0,ctngas]
	ctnsink=[0,ctnsink]


	data=ptrarr(n_elements(idx))

	if(idx[0] ne -1) then begin
		for i=0L,n_elements(idx)-1 do begin
			gid=list[idx[i]]
			hid=sub_host[gid]
			output=path+'/gid_'+string(gid,format='(i07)')+'.sav'
			tnstar=ctnstar[gid]
			tndm=ctndm[gid]
			tngas=ctngas[gid]
			tnsink=ctnsink[gid]

		
			skip=long64(hid+1)*bhalo+long64(gid+1)*bsub+(tnstar+tndm)*bpart+tngas*bgas+tnsink*bsink
			print,'Galaxy ID : '+string(gid,format='(i07)')+' Skip bytes : '+string(skip,format='(i15)')
			dm=-1
			star=-1
			gas=-1
			mbh=-1
			if(subhalo[gid].ndm gt 0) then dm=replicate(tpart,subhalo[gid].ndm)
			if(subhalo[gid].nstar gt 0) then star=replicate(tpart,subhalo[gid].nstar)
			if(subhalo[gid].ngas gt 0) then gas=replicate(tgas,subhalo[gid].ngas)
			if(subhalo[gid].nsink gt 0) then mbh=replicate(tsink,subhalo[gid].nsink)


	
			point_lun,10,skip
			if(subhalo[gid].ndm gt 0) then begin
				readu,10,dm
			endif

			if(subhalo[gid].ngas gt 0) then begin
				readu,10,gas
			endif	

			if(subhalo[gid].nsink gt 0) then begin
				readu,10,mbh
			endif
			if(subhalo[gid].nstar gt 0) then begin
				readu,10,star
			endif

			ahalo={id:gid,dm:dm,gas:gas,star:star,mbh:mbh}					
			data[i]=ptr_new(ahalo)
			if keyword_set(dump_file) then save,file=output,dm,star,gas,sink,info

		endfor
	endif
	
	
	close,10
end
