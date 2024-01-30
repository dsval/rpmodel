#' rpmodel.grid
#'
#' rpmodel looping through pixels in parallel and writing on the go 
#' @param  same as rpmodel, tc, vpd, fapar, ppfd, elev, soilm, meanalpha as grid type objects
#' @import raster  
#' @keywords rpmodel
#' @export
#' @examples
#' rpmodel.grid()	

rpmodel.grid<-function(tc, vpd, co2, fapar, ppfd, soilm,meanalpha,elev = NA,AI=NA ,
	beta = 146.0,	do_ftemp_kphio = TRUE, do_soilmstress = TRUE,inmem=FALSE,outdir=getwd()){
	
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster(ncores) if you want to run in parallel!!')
		
	}
	###############################################################################################
	# 00. create array for results
	###############################################################################################
	#rasterOptions(maxmemory=1e7, tmptime = 24, chunksize = 1e6,todisk = FALSE, overwrite=TRUE, tolerance = 0.5)
	y<-as.numeric(unique(format(getZ(tc),'%Y')))
	ny <- nlayers(tc)
	################ arrays for well watered
	if (!do_soilmstress){
		soilm=calc(tc,fun= function(x){x/x})
		meanalpha = soilm[[1]]
		#soilm=reclassify(tc, cbind(-Inf,+Inf , 1), right=NA);gc()
	}
	
	
	# gpp gC/m2
	gpp<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(gpp)<-extent(elev)
	gpp<-setZ(gpp,getZ(tc))
	# lue mmol C / mol photons
	lue<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(lue)<-extent(elev)
	lue<-setZ(lue,getZ(tc))
	# transpiration mm
	Tr<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(Tr)<-extent(elev)
	Tr<-setZ(Tr,getZ(tc))
	# wue mmolC / mol H2O
	wue<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(wue)<-extent(elev)
	wue<-setZ(wue,getZ(tc))
	gc()
	setwd(outdir)
	if(inmem){
		tc<-readAll(tc)
		vap<-readAll(vap)
		elev<-readAll(elev)
		fapar<-readAll(fa)
		ppfd<-readAll(ppfd)
		soilm<-readAll(soilm)
		meanalpha<-readAll(meanalpha)
		
	}
	
	###############################################################################################
	# 01. set the clusters for parallel computing
	###############################################################################################	
	cl <- getCluster()
	on.exit( returnCluster() )
	nodes <- length(cl)
	bs<- blockSize(tc, minblocks=nodes)
	parallel::clusterEvalQ(cl, library("rpmodel"))
	#pmodel<-Vectorize(rpmodel,c("tc","vpd","co2","fapar","ppfd","meanalpha","soilm"))
	parallel::clusterExport(cl, varlist=c("tc","vpd","co2",'elev','fapar',"ppfd",'beta','soilm','AI','meanalpha','do_soilmstress','ny','y','bs','calc_phi0'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 3)
	###############################################################################################
	# 02. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun <- function(i) {
		tcrow<-split(getValues(tc,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		vpdrow<-split(getValues(vpd,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		elevrow<-split(getValues(elev,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		AIrow<-split(getValues(AI,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		faparrow<-split(getValues(fapar,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		ppfdrow<-split(getValues(ppfd,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		soilmrow<-split(getValues(soilm,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		alpharow<-split(getValues(meanalpha,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		co2row<-lapply(1:(ncol(elev)*bs$nrows[i]),function(j) co2)
		# define function to get transpiration
		calc_transp<-function(mat,vpd,tc,ppfd ,elev){
			press<-calc_patm(elev)
			#assimilation in mmol
			A_mmol<-(mat$gpp/12.0107)*1000
			gs_mol<-as.numeric(mat$gs)
			gs_mol[gs_mol<0]<-NA
			T_mol<-1.6*gs_mol*(vpd)
			T_mm<-T_mol*(18/density_h2o(tc,press))
			T_mm[T_mm>1000]<-NA
			wue<-A_mmol/T_mol
			lue<-A_mmol/ppfd
			result<-t(as.matrix(data.frame(gpp=mat$gpp,T_mm=T_mm,wue=wue,lue=lue)))
			#list(gpp=mat$gpp,T_mm=T_mm,wue=wue,lue=lue)
			#return(rbind(mat,T_mm,wue))
			return(result)
		}
		
		# apply vectorized rpmodel		
		if(do_ftemp_kphio){
			kphiorow  =mapply(calc_phi0,tc=tcrow,AI=AIrow,SIMPLIFY = FALSE) 
		}
						
										
														
																						
		result<-mapply(rpmodel, 
			tc             = tcrow,           # temperature, deg C
			vpd            = vpdrow,         # Pa,
			co2            = co2row,          # ppm,
			elv            = elevrow,        # m.a.s.l.,
			kphio          = kphiorow,         # quantum yield efficiency,
			beta           = 146,         # unit cost ratio a/b,
			fapar          = faparrow,      # fraction  ,
			ppfd           = ppfdrow,      # mol/m2/d month?,
			soilm = soilmrow,
			meanalpha=alpharow, MoreArgs = list(do_ftemp_kphio = FALSE,do_soilmstress = do_soilmstress),SIMPLIFY = FALSE)
		# calc transpiration and wue
		result<-lapply(result,as.data.frame)
		result<-mapply(calc_transp,result,vpdrow,tcrow,ppfdrow,elevrow, SIMPLIFY = FALSE)
		return(result)
	}
	###############################################################################################
	# 03. send tasks to the nodes
	###############################################################################################
	for (i in 1:nodes) {
		parallel:::sendCall(cl[[i]], clFun, i, tag=i)
	}
	###############################################################################################
	# 04. write to the disk on the go, or save to the ram
	###############################################################################################
	if(!inmem){
		gpp<-writeStart(gpp,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'gpp',".","nc"),
			format="CDF",overwrite=TRUE,varname='gpp', varunit='gC/m2',longname='gross primary production',
			xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		Tr<-writeStart(Tr,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'Tr',".","nc"),
			format="CDF",overwrite=TRUE,varname='Tr', varunit='mm',longname='Transpiration',
			xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		wue<-writeStart(wue,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'wue',".","nc"),
			format="CDF",overwrite=TRUE,varname='wue', varunit='mmol/mol',longname='water use efficiency',
			xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		lue<-writeStart(lue,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",'lue',".","nc"),
			format="CDF",overwrite=TRUE,varname='wue', varunit='mmol/mol',longname='light use efficiency',
			xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		
		
	}else {
		matgpp <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
		matTr<- matrix(ncol=nlayers(tc), nrow=ncell(tc))
		matwue <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
		matlue <- matrix(ncol=nlayers(tc), nrow=ncell(tc))
		endind<-cumsum(bs$nrows*tc@ncols)
		startind<-c(1,endind+1)    
	}
	
	# mbl1<-do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat[10,])}))
	###############################################################################################
	# 05. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- parallel:::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('cluster error:',"\n",d$value$value)
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!inmem) {
			gpp <- writeValues(gpp,do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['gpp',])})), bs$row[b])
			Tr <- writeValues(Tr, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['T_mm',])})), bs$row[b])
			wue <- writeValues(wue, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['wue',])})), bs$row[b])
			lue <- writeValues(lue, do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['lue',])})), bs$row[b])
			
			
		} else {
			
			matgpp[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['gpp',])}))
			matTr[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['T_mm',])}))
			matwue[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['wue',])}))
			matlue[startind[b]:endind[b],] <- do.call(rbind,lapply(d$value$value,function(mat){as.numeric(mat['lue',])}))
			
		}
		
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
		}
		gc()
		setTxtProgressBar(pb,i)
	}
	###############################################################################################
	# 06. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!inmem) {
		gpp <- writeStop(gpp)
		Tr <- writeStop(Tr)
		wue <- writeStop(wue)
		lue <- writeStop(lue)
		
	} else {
		# gpp
		gpp<-setValues(gpp,matgpp)
		# transpiration
		Tr<-setValues(Tr,matTr)
		# water use efficiency
		wue<-setValues(wue,matwue)
		# light use efficiency
		lue<-setValues(lue,matlue)
		
	}
	close(pb)
	gc()
	return(list(gpp=gpp,lue=lue,Tr=Tr,wue=wue))
}
