pro read_target_halo,nout=nout,list=list,dump_file=dump_file,data=data,info=info


	

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
	file3='FoF_Data/FoF.'+snout1+'/background_ptl.'+snout1
	file4='FoF_Data/halo_'+snout1+'.sav'
	tp=file_test(file4)
	if(tp[0] eq 0) then begin
		 make_halo_catalouge,nout=nout	
	endif

	restore,file4

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
	openr,20,file3



	print,'Nhalo total: ',nhalo
	print,'Nsub total: ',nsub



	idx=where(list ge 0)
	ctnstar=total(subhalo[0:nsub-2].nstar,/integer,/cumulative)
	ctndm=total(subhalo[0:nsub-2].ndm,/integer,/cumulative)
	ctngas=total(subhalo[0:nsub-2].ngas,/integer,/cumulative)
	ctnsink=total(subhalo[0:nsub-2].nsink,/integer,/cumulative)
	ctnsub=total(halo[0:nhalo-2].nsub,/integer,/cumulative)

	ctnstar=[0,ctnstar]
	ctndm=[0,ctndm]
	ctngas=[0,ctngas]
	ctnsink=[0,ctnsink]
	ctnsub=[0,ctnsub]


	bndm=unbound.ndm
	bnstar=unbound.nstar
	bngas=unbound.ngas
	bnsink=unbound.nsink


	bctndm=total(bndm[0:nhalo-2],/integer,/cumulative)
	bctnstar=total(bnstar[0:nhalo-2],/integer,/cumulative)
	bctngas=total(bngas[0:nhalo-2],/integer,/cumulative)
	bctnsink=total(bnsink[0:nhalo-2],/integer,/cumulative)
	bctndm=[0,bctndm]
	bctnstar=[0,bctnstar]
	bctngas=[0,bctngas]
	bctnsink=[0,bctnsink]
	data=ptrarr(n_elements(idx))

	if(idx[0] ne -1) then begin
		for i=0L,n_elements(idx)-1 do begin
			hid=list[idx[i]]
			nsub=halo[hid].nsub
			gid=ctnsub[hid]
			output=path+'/hid_'+string(hid,format='(i07)')+'.sav'
			tnstar=ctnstar[gid]
			tndm=ctndm[gid]
			tngas=ctngas[gid]
			tnsink=ctnsink[gid]

			btndm=bctndm[hid]
			btnstar=bctnstar[hid]
			btnsink=bctnsink[hid]
			btngas=bctngas[hid]

			skip=long64(hid+1)*bhalo+long64(gid)*bsub+(tnstar+tndm)*bpart+tngas*bgas+tnsink*bsink
			skip1=long64(hid+1)*bsub+(btnstar+btndm)*bpart+btngas*bgas+btnsink*bsink
			print,'Halo ID : '+string(hid,format='(i07)')+' Skip bytes : '+string(skip,format='(i15)')
			adm=-1
			astar=-1
			agas=-1
			ambh=-1
			if(halo[hid].ndm gt 0) then adm=replicate(tpart,halo[hid].ndm)
			if(halo[hid].nstar gt 0) then astar=replicate(tpart,halo[hid].nstar)
			if(halo[hid].ngas gt 0) then agas=replicate(tgas,halo[hid].ngas)
			if(halo[hid].nsink gt 0) then ambh=replicate(tsink,halo[hid].nsink)
			id_dm=lonarr(2,nsub+1)
			id_star=lonarr(2,nsub+1)
			id_gas=lonarr(2,nsub+1)
			id_mbh=lonarr(2,nsub+1)
			id_galaxy=lonarr(nsub+1)
			id_galaxy[0]=-1
			id_dm[*]=-1
			id_star[*]=-1
			id_mbh[*]=-1
			id_gas[*]=-1

			if(unbound[hid].ndm gt 0) then begin
				id_dm[0,0]=0
				id_dm[1,0]=unbound[hid].ndm-1
			endif
			if(unbound[hid].nstar gt 0) then begin
				id_star[0,0]=0
				id_star[1,0]=unbound[hid].nstar-1
			endif
			if(unbound[hid].ngas gt 0) then begin
				id_gas[0,0]=0
				id_gas[1,0]=unbound[hid].ngas-1
			endif
			if(unbound[hid].nsink gt 0) then begin
				id_mbh[0,0]=0
				id_mbh[1,0]=unbound[hid].nsink-1
			endif
			lidm=0L
			listar=0L
			ligas=0L
			lisink=0L

			point_lun,20,skip1
			if(unbound[hid].ndm gt 0) then begin
				dm=replicate(tpart,unbound[hid].ndm)
				readu,20,dm
				adm[0:id_dm[1,0]]=dm
				lidm=id_dm[1,0]+1
			endif
			if(unbound[hid].ngas gt 0) then begin
				gas=replicate(tgas,unbound[hid].ngas)
				readu,20,gas
				agas[0:id_gas[1,0]]=gas
				ligas=id_gas[1,0]+1
			endif
			if(unbound[hid].nsink gt 0) then begin
				sink=replicate(tsink,unbound[hid].nsink)
				readu,20,sink
				ambh[0:id_mbh[1,0]]=sink
				lisink=id_mbh[1,0]+1
			endif
			if(unbound[hid].nstar gt 0) then begin
				star=replicate(tpart,unbound[hid].nstar)
				readu,20,star
				astar[0:id_star[1,0]]=star
				listar=id_star[1,0]+1
			endif

	
			point_lun,10,skip
			for j=0,nsub-1 do begin
				id_galaxy[j+1]=gid
				readu,10,tsub
		;		skip=skip+bsub
				if(subhalo[gid].ndm gt 0) then begin
					dm=replicate(tpart,subhalo[gid].ndm)
					readu,10,dm
					adm[lidm:lidm+subhalo[gid].ndm-1]=dm
					id_dm[0,j+1]=lidm
					lidm=lidm+subhalo[gid].ndm
					id_dm[1,j+1]=lidm-1
				endif

				if(subhalo[gid].ngas gt 0) then begin
					gas=replicate(tgas,subhalo[gid].ngas)
					readu,10,gas
					gas.temp=gas.temp/gas.density
					agas[ligas:ligas+subhalo[gid].ngas-1]=gas
					id_gas[0,j+1]=ligas
					ligas=ligas+subhalo[gid].ngas
					id_gas[1,j+1]=ligas-1
				endif	

				if(subhalo[gid].nsink gt 0) then begin
					sink=replicate(tsink,subhalo[gid].nsink)
					readu,10,sink
					ambh[lisink:lisink+subhalo[gid].nsink-1]=sink
					id_mbh[0,j+1]=lisink
					lisink=lisink+subhalo[gid].nsink
					id_mbh[1,j+1]=lisink-1
				endif
				if(subhalo[gid].nstar gt 0) then begin
					star=replicate(tpart,subhalo[gid].nstar)
					readu,10,star
					astar[listar:listar+subhalo[gid].nstar-1]=star
					id_star[0,j+1]=listar
					listar=listar+subhalo[gid].nstar
					id_star[1,j+1]=listar-1
				endif
;				skip_lun,10,bsub
				gid=gid+1
			endfor

			ahalo={id:hid,id_galaxy:id_galaxy,id_dm:id_dm,id_gas:id_gas,id_star:id_star,id_mbh:id_mbh,dm:adm,gas:agas,star:astar,mbh:ambh}					
			data[i]=ptr_new(ahalo)
			if keyword_set(dump_file) then save,file=output,adm,astar,agas,ambh,info,hid,id_dm,id_gas,id_mbh,id_star,id_galaxy

		endfor
	endif
	
	
	close,10
	close,20
end
