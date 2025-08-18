import os
import glob
import sys
import numpy as np
import time



class read_PGalF():

	def __init__(self,path):
		self.fllist=glob.glob(path+'GALCATALOG.LIST*')
		self.flcen=glob.glob(path+'GALFIND.CENTER*')
		self.fldat=glob.glob(path+'GALFIND.DATA*')
		self.flbytename='./byte2skip4subhalo.txt'
		self.flbyte=glob.glob(self.flbytename)

		if not bool(self.fllist)*bool(self.flcen)*bool(self.fldat):
			print('Not all files ready')
			return

		self.chunksize=4*6+8*11
		self.chunksize1=8*10+4*3+2*1+1*2
		self.chunksize2=8*9+4*8
		self.chunksize3=8*20+4*2


	def read_list(self):
		filesize=os.path.getsize(self.fllist[0])
		f=open(self.fllist[0],'rb')

		nhalo=0
		nsubhalo=0
		byteleft=filesize*1
		while byteleft>0:
			[nsub,ndm,nst,nsink,ngas,nall]=np.fromfile(f,dtype=np.int32,count=6)
			[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
			pos=np.fromfile(f,dtype=np.float64,count=3)
			vel=np.fromfile(f,dtype=np.float64,count=3)

			byteleft-=self.chunksize

			nhalo+=1
			nsubhalo+=nsub
			for i in range(nsub):
				[ndm,ngas,nsink,nst,npall,dummy4]=np.fromfile(f,dtype=np.int32,count=6)
				[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
				pos=np.fromfile(f,dtype=np.float64,count=3)
				vel=np.fromfile(f,dtype=np.float64,count=3)

			byteleft-=self.chunksize*nsub

		self.nhalo=nhalo
		self.nsubhalo=nsubhalo
		print(f'nhalo={nhalo}, nsubhalo={nsubhalo}')
		f.close()


	def create_byte2skip_file(self):
		if not bool(self.flbyte):
			filesize=os.path.getsize(self.fllist[0])
			f=open(self.fllist[0],'rb')

			f1=open(self.flbytename,'w')
			f1.write('subhalo-id\tbyte-to-skip\n')

			nsubhalo=0
			byteleft=filesize*1
			byte2skip=0
			while byteleft>0:
				[nsub,ndm,nst,nsink,ngas,nall]=np.fromfile(f,dtype=np.int32,count=6)
				[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
				pos=np.fromfile(f,dtype=np.float64,count=3)
				vel=np.fromfile(f,dtype=np.float64,count=3)

				byte2skip+=self.chunksize

				for i in range(nsub):
					[ndm,ngas,nsink,nst,npall,dummy4]=np.fromfile(f,dtype=np.int32,count=6)
					[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
					pos=np.fromfile(f,dtype=np.float64,count=3)
					vel=np.fromfile(f,dtype=np.float64,count=3)

					f1.write(f'{nsubhalo+i}\t{byte2skip}\n')

					byte2skip+=self.chunksize+self.chunksize1*(ndm+nst)+self.chunksize2*ngas+self.chunksize3*nsink

				nsubhalo+=nsub
				byteleft-=self.chunksize*(nsub+1)

			f.close()
			f1.close()


	def read_data(self):
		filesize=os.path.getsize(self.fldat[0])
		f=open(self.fldat[0],'rb')

		nhalo=0
		nsubhalo=0
		byteleft=filesize*1
		tstt=time.time()
		while byteleft>0:
			[nsub,ndm,nst,nsink,ngas,nall]=np.fromfile(f,dtype=np.int32,count=6)
			[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
			pos=np.fromfile(f,dtype=np.float64,count=3)
			vel=np.fromfile(f,dtype=np.float64,count=3)

			nhalo+=1

			ndm1=0
			ngas1=0
			nsink1=0
			nst1=0
			for i in range(nsub):
				[ndm,ngas,nsink,nst,npall,dum]=np.fromfile(f,dtype=np.int32,count=6)
				[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
				pos=np.fromfile(f,dtype=np.float64,count=3)
				vel=np.fromfile(f,dtype=np.float64,count=3)

				ndm1+=ndm
				ngas1=+ngas
				nsink1+=nsink
				nst1+=nst

				for j in range(ndm):
					pos=np.fromfile(f,dtype=np.float64,count=3)
					vel=np.fromfile(f,dtype=np.float64,count=3)
					[mass,mass0,tp,zp]=np.fromfile(f,dtype=np.float64,count=4)
					idx=np.fromfile(f,dtype=np.int32,count=1)
					[family,tag]=np.fromfile(f,dtype=np.int8,count=2)
					dum=np.fromfile(f,dtype=np.int16,count=1)
					[level,dum1]=np.fromfile(f,dtype=np.int32,count=2)
					#print(j,pos,vel,mass,mass0,tp,zp,idx,family,tag,dum,level,dum1)
					#print(nsubhalo+i,j,pos,vel,mass,mass0,tp,zp,idx,family,tag,dum,level,dum1)

				for j in range(ngas):
					pos=np.fromfile(f,dtype=np.float64,count=3)
					dx=np.fromfile(f,dtype=np.float64,count=1)
					vel=np.fromfile(f,dtype=np.float32,count=3)
					dum0=np.fromfile(f,dtype=np.float32,count=1)
					density=np.fromfile(f,dtype=np.float64,count=1)
					[temp,metal,mass,dum1]=np.fromfile(f,dtype=np.float32,count=4)
					potential=np.fromfile(f,dtype=np.float64,count=1)
					force=np.fromfile(f,dtype=np.float64,count=3)
					#print(j,pos,dx,vel,dum0,density,temp,metal,mass,dum1,potential,force)

				for j in range(nsink):
					pos=np.fromfile(f,dtype=np.float64,count=3)
					vel=np.fromfile(f,dtype=np.float64,count=3)
					[mass,tbirth]=np.fromfile(f,dtype=np.float64,count=2)
					angm=np.fromfile(f,dtype=np.float64,count=3)
					ang=np.fromfile(f,dtype=np.float64,count=3)
					dmsmbh=np.fromfile(f,dtype=np.float64,count=3)
					[esave,smag,eps]=np.fromfile(f,dtype=np.float64,count=3)
					[idx,dum]=np.fromfile(f,dtype=np.int32,count=2)
					#print(j,pos,vel,mass,tbirth,angm,ang,dmsmbh,esave,smag,eps,idx,dum)

				for j in range(nst):
					pos=np.fromfile(f,dtype=np.float64,count=3)
					vel=np.fromfile(f,dtype=np.float64,count=3)
					[mass,mass0,tp,zp]=np.fromfile(f,dtype=np.float64,count=4)
					idx=np.fromfile(f,dtype=np.int32,count=1)
					[family,tag]=np.fromfile(f,dtype=np.int8,count=2)
					dum=np.fromfile(f,dtype=np.int16,count=1)
					[level,dum1]=np.fromfile(f,dtype=np.int32,count=2)
					#print(j,pos,vel,mass,mass0,tp,zp,idx,family,tag,dum,level,dum1)
			
			nsubhalo+=nsub
			byteread=self.chunksize*(1+nsub)+self.chunksize1*(ndm1+nst1)+self.chunksize2*ngas1+self.chunksize3*nsink1
			byteleft-=byteread
		tend=time.time()
		print(f'The time elapsed to read the data file is {(tend-tstt)/60} min.')
		f.close()


	def read_a_subhalo(self,subhaloid):
		if not bool(self.flbyte):
			print(f'{self.flbytename} does not exist. First run INSTANT-NAME.create_byte2skip_file().')
			return

		subhaloidarr,byte2skiparr=np.loadtxt(self.flbyte[0],dtype=int,unpack=True,skiprows=1)
		byteskip=byte2skiparr[subhaloid]

		f=open(self.fldat[0],'rb')
		f.read(byteskip)
		[ndm,ngas,nsink,nst,npall,dum]=np.fromfile(f,dtype=np.int32,count=6)
		[mtot,mdm,mgas,msink,mst]=np.fromfile(f,dtype=np.float64,count=5)
		pos=np.fromfile(f,dtype=np.float64,count=3)
		vel=np.fromfile(f,dtype=np.float64,count=3)
		#print(ndm,ngas,nsink,nst,npall,mtot,mdm,mgas,msink,mst,pos,vel)

		dtype=np.dtype([('ndm','i4'),('ngas','i4'),('nsink','i4'),('nst','i4'),('npall','i4'),\
						('mtot','f8'),('mdm','f8'),('mgas','f8'),('msink','f8'),('mst','f8'),\
						('pos','f8',(3)),('vel','f8',(3))])	
		SUMMARY=np.array([(ndm,ngas,nsink,nst,npall,mtot,mdm,mgas,msink,mst,pos,vel)],dtype=dtype)


		dtype=np.dtype([('pos','f8',(3)),('vel','f8',(3)),('mass','f8'),('mass0','f8'),('tp','f8'),('zp','f8'),\
						('id','i4'),('family','i1'),('tag','i1'),('level','i4')])
		DM=np.zeros(ndm,dtype=dtype)
		for j in range(ndm):
			pos=np.fromfile(f,dtype=np.float64,count=3)
			vel=np.fromfile(f,dtype=np.float64,count=3)
			[mass,mass0,tp,zp]=np.fromfile(f,dtype=np.float64,count=4)
			idx=np.fromfile(f,dtype=np.int32,count=1)
			[family,tag]=np.fromfile(f,dtype=np.int8,count=2)
			dum=np.fromfile(f,dtype=np.int16,count=1)
			[level,dum1]=np.fromfile(f,dtype=np.int32,count=2)
			#print(j,pos,vel,mass,mass0,tp,zp,idx,family,tag,dum,level,dum1)
			DM[j]=np.array([(pos,vel,mass,mass0,tp,zp,idx,family,tag,level)],dtype=dtype)


		dtype=np.dtype([('pos','f8',(3)),('dx','f8'),('vel','f4',(3)),('density','f8'),\
						('temp','f4'),('metal','f4'),('mass','f4'),\
						('potential','f8'),('force','f8',(3))])
		GAS=np.zeros(ngas,dtype=dtype)
		for j in range(ngas):
			pos=np.fromfile(f,dtype=np.float64,count=3)
			dx=np.fromfile(f,dtype=np.float64,count=1)
			vel=np.fromfile(f,dtype=np.float32,count=3)
			dum0=np.fromfile(f,dtype=np.float32,count=1)
			density=np.fromfile(f,dtype=np.float64,count=1)
			[temp,metal,mass,dum1]=np.fromfile(f,dtype=np.float32,count=4)
			potential=np.fromfile(f,dtype=np.float64,count=1)
			force=np.fromfile(f,dtype=np.float64,count=3)
			#print(j,pos,dx,vel,dum0,density,temp,metal,mass,dum1,potential,force)
			GAS[j]=np.array([(pos,dx,vel,density,temp,metal,mass,potential,force)],dtype=dtype)

		dtype=np.dtype([('pos','f8',(3)),('vel','f8',(3)),('mass','f8'),('tbirth','f8'),('angm','f8'),('ang','f8'),\
						('dmsmbh','f8',(3)),('esave','f8'),('smag','f8'),('eps','f8'),('id','i4')])
		SINK=np.zeros(ndm,dtype=dtype)
		for j in range(nsink):
			pos=np.fromfile(f,dtype=np.float64,count=3)
			vel=np.fromfile(f,dtype=np.float64,count=3)
			[mass,tbirth]=np.fromfile(f,dtype=np.float64,count=2)
			angm=np.fromfile(f,dtype=np.float64,count=3)
			ang=np.fromfile(f,dtype=np.float64,count=3)
			dmsmbh=np.fromfile(f,dtype=np.float64,count=3)
			[esave,smag,eps]=np.fromfile(f,dtype=np.float64,count=3)
			[idx,dum]=np.fromfile(f,dtype=np.int32,count=2)
			#print(j,pos,vel,mass,tbirth,angm,ang,dmsmbh,esave,smag,eps,idx,dum)
			SINK[j]=np.array([(pos,vel,mass,tbirth,angm,ang,dmsmbh,esave,smag,eps,idx)],dtype=dtype)

		dtype=np.dtype([('pos','f8',(3)),('vel','f8',(3)),('mass','f8'),('mass0','f8'),('tp','f8'),('zp','f8'),\
						('id','i4'),('family','i1'),('tag','i1'),('level','i4')])
		STAR=np.zeros(ndm,dtype=dtype)
		for j in range(nst):
			pos=np.fromfile(f,dtype=np.float64,count=3)
			vel=np.fromfile(f,dtype=np.float64,count=3)
			[mass,mass0,tp,zp]=np.fromfile(f,dtype=np.float64,count=4)
			idx=np.fromfile(f,dtype=np.int32,count=1)
			[family,tag]=np.fromfile(f,dtype=np.int8,count=2)
			dum=np.fromfile(f,dtype=np.int16,count=1)
			[level,dum1]=np.fromfile(f,dtype=np.int32,count=2)
			#print(j,pos,vel,mass,mass0,tp,zp,idx,family,tag,dum,level,dum1)
			STAR[j]=np.array([(pos,vel,mass,mass0,tp,zp,idx,family,tag,level)],dtype=dtype)

		f.close()

		return SUMMARY,DM,GAS,SINK,STAR



path='/scratch/jaehyun/Darwin/02_hydro_test/FoF_Data/FoF.00052/'
PGalF=read_PGalF(path)
PGalF.read_list()	# print nhalo and nsubhalo
PGalF.create_byte2skip_file()	# create an ascii file that includes subhaloid and byte to skip to jump to where the target subhalo information starts
PGalF.read_data()	# simple read for now; print the elapsed time to read the whole data
subhaloid=4
SUMMARY,DM,GAS,SINK,STAR=PGalF.read_a_subhalo(subhaloid)	# return information of a given subhaloid (subhaloid starts from 0)
