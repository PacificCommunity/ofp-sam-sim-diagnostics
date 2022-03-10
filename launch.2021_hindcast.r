

# Nicholas Ducharme-Barth
# 09/07/2021
# 2021_hindcast
# apply hindcasting according to Kell et al. 2021
# https://academic.oup.com/icesjms/advance-article/doi/10.1093/icesjms/fsab104/6296435 

#________________________________________________________________________________________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ssh)
	library(FLR4MFCL)
	library(ggplot2)
	library(ggthemes)
	library(frqit)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# set working directory & define directories
	proj.dir = "C:/Users/nicholasd/HOME/SPC/SPC_SAM/2021_SC/"
	setwd(proj.dir)
	dir.model_runs = paste0(proj.dir,"SWO/Assessment/Model_Runs/")
	dir.condor = paste0(proj.dir,"SWO/Assessment/condor_files/")
	dir.stem = paste0(proj.dir,"SWO/Assessment/Model_Runs/2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp_selexExplore/v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1/")
	dir.frq = paste0(proj.dir,"SWO/Assessment/Data_Prep/frq_file/")
	dir.launch = "/home/nicholasd/"
	dir.mfcl.launch = "/home/nicholasd/mfcl/"
	submit.realm = "NOU"
	supress_mfcl_IO = TRUE

#________________________________________________________________________________________________________________________________________________________________________________________________________
# create array of testing runs
	testing.options = expand.grid(model = "v2080",endYear = seq(from=2019,to=2009,by=-1),type=c("cpue"), stringsAsFactors = FALSE)

#____________________________________________________________________________________________________________
# launch condor runs
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}
		jobs.group = "2021_hindcast/"

		# create directory on launch machine
			ssh_exec_wait(session, command = paste0("mkdir -p ",dir.launch,jobs.group))

		# iterate across testing runs
			for(i in 1:nrow(testing.options))
			{
				runname = paste0(testing.options[i,],collapse = "_")
				model.run.dir = paste0(dir.model_runs,jobs.group,runname,"/")
				# create new directory for model run
					if (! dir.exists(model.run.dir))dir.create(model.run.dir,recursive=TRUE)
					ssh_exec_wait(session, command = paste0("mkdir ",dir.launch,jobs.group,runname))
				# transfer files
					FileList=c("condor_doitall.swo","mfcl.cfg","selblocks.dat","swo.frq","swo.ini")
					file.copy(paste0(dir.stem,FileList),model.run.dir,overwrite=TRUE)
				# rename files
					# file.rename(paste0(model.run.dir,c("doitall.swo")),paste0(model.run.dir,c("condor_doitall.swo")))

				# modify frq, downweight data after end year
					if(testing.options$endYear[i] < 2019)
					{
						tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
					  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[fishery>13 & year>testing.options$endYear[i],penalty:=0.001]
						cateffpen(tmp.frqit) = as.data.frame(catch.dt)
					  	writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
					  	tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
					  	tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
					  	writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))	
					}

				# copy condor files	
					file.copy(paste0(dir.condor,"DAG/",c("condor.dag","mfcl.X86_64.LINUX.bat","DAGcondorLIN.sub","changeTars.sh","change0.sh")),model.run.dir,overwrite=TRUE)
					if(submit.realm != "SUV")
					{
						if(submit.realm == "NOU")
						{
							tmp.sub = readLines(paste0(model.run.dir,"DAGcondorLIN.sub"),warn=FALSE)
							pointer = grep("Requirements",tmp.sub)
							tmp.sub[pointer] = 'requirements = ((Realm == "NOU") && (OpSys == "LINUX") && (Machine =!= "NOUOFPCALC02.corp.spc.int") && (Machine =!= "NOUOFPCAND01") && (Machine =!= "NOUOFPCAND02") && (Machine =!= "NOUOFPCAND03") && (Machine =!= "NOUOFPCAND04") && (Machine =!= "NOUOFPCAND05") && (Machine =!= "NOUOFPCAND06") && (Machine =!= "NOUOFPCAND07") && (Machine =!= "NOUOFPCAND08") && (Machine =!= "nouofpcand09.corp.spc.int") && (Machine =!= "nouofpcand10") && (Machine =!= "nouofpcand11.corp.spc.int") && (Machine =!= "nouofpcand12"))'
							writeLines(tmp.sub,con=paste0(model.run.dir,"DAGcondorLIN.sub"))
							rm(list=c("tmp.sub","pointer"))
						} else {
							tmp.sub = readLines(paste0(model.run.dir,"DAGcondorLIN.sub"),warn=FALSE)
							pointer = grep("Requirements",tmp.sub)
							tmp.sub[pointer] = 'requirements = ((OpSys == "LINUX") && (Machine =!= "NOUOFPCALC02.corp.spc.int") && (Machine =!= "NOUOFPCAND01") && (Machine =!= "NOUOFPCAND02") && (Machine =!= "NOUOFPCAND03") && (Machine =!= "NOUOFPCAND04") && (Machine =!= "NOUOFPCAND05") && (Machine =!= "NOUOFPCAND06") && (Machine =!= "NOUOFPCAND07") && (Machine =!= "NOUOFPCAND08") && (Machine =!= "nouofpcand09.corp.spc.int") && (Machine =!= "nouofpcand10") && (Machine =!= "nouofpcand11.corp.spc.int") && (Machine =!= "nouofpcand12"))'
							writeLines(tmp.sub,con=paste0(model.run.dir,"DAGcondorLIN.sub"))
							rm(list=c("tmp.sub","pointer"))
						}
					}

				# modify .bat file
					tmp.bat = readLines(paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"),warn=FALSE)
					pointer = grep("./doitall.condor",tmp.bat,fixed=TRUE)
					tmp.bat = c(tmp.bat[1:(pointer-1)],"start=`date +%s`",tmp.bat[pointer],"end=`date +%s`","runtime=$((end-start))","echo $runtime","",'if [ -f "runtime.txt" ]',"then","  touch runtime.txt","  echo $runtime >> runtime.txt","else","  echo $runtime > runtime.txt","fi",tmp.bat[(pointer+1):length(tmp.bat)])
					pointer = grep("# Repack select files that are needed for the next DAG phase",tmp.bat) + 1
					tmp.bat[pointer] = "        tar -czf End.tar.gz 'doitall.condor' 'mfclo64' 'mfcl.cfg' 'selblocks.dat' 'swo.frq' 'swo.ini' 'runtime.txt' *.par"
					pointer = grep("# Repack select files that are needed for the final DAG phase",tmp.bat) + 1
					tmp.bat[pointer] = "        tar -czf End.tar.gz 'doitall.condor' 'mfclo64' 'mfcl.cfg' 'selblocks.dat' 'swo.frq' 'swo.ini' 'plot-16.par.rep' 'runtime.txt' 'length.fit' 'weight.fit' *.par"
					pointer = grep("# Repack select files so that this is all that needs to be exported",tmp.bat) + 1
					tmp.bat[pointer] = "    tar -czf End.tar.gz '16.par' 'plot-16.par.rep' 'test_plot_output' 'sorted_gradient.rpt' 'xinit.rpt' 'new_cor_report' 'neigenvalues' 'runtime.txt' 'length.fit' 'weight.fit' swo.var"
					if(supress_mfcl_IO)
					{
						tmp.bat = gsub("./doitall.condor","./doitall.condor &>/dev/null",tmp.bat,fixed=TRUE)
					}
					writeLines(tmp.bat,con=paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"))
					rm(list=c("tmp.bat","pointer"))

				# make tar with the necessary files
					TarList=c("condor_doitall.swo","mfcl.cfg","selblocks.dat","swo.frq","swo.ini")
					shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -czf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
				# send all files to the launch machine
			       	scp_upload(session,files=paste0(model.run.dir,c("Start.tar.gz","condor.dag","mfcl.X86_64.LINUX.bat","DAGcondorLIN.sub","changeTars.sh","change0.sh")),to=paste0(dir.launch,jobs.group,runname))
				# copy mfcl to launch machine run directory
			       	if(testing.options$model[i] == "v2070")
			       	{
			       		ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'2070/mfclo64 ',paste0(dir.launch,jobs.group,runname)))
			       	} else if (testing.options$model[i] == "v2082"){
			       		ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'v2082/mfclo64 ',paste0(dir.launch,jobs.group,runname)))
			       	} else {
						ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'v2080_dev_20210503/mfclo64 ',paste0(dir.launch,jobs.group,runname)))			       		
			       	}
			    # remake tar file to include mfcl
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,runname),'; tar -xzf Start.tar.gz; rm Start.tar.gz'))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,runname),'; tar -czf Start.tar.gz ',paste(c(TarList,"mfclo64"),collapse=' ')))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,runname),'; rm ',paste(c(TarList,"mfclo64"),collapse=' ')))
				# launch condor job
			        ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,runname),'; dos2unix change*.sh; chmod 700 change*.sh; dos2unix mfcl.X86_64.LINUX.bat; chmod 700 mfcl.X86_64.LINUX.bat; condor_submit_dag condor.dag'))   
	
				# clean-up
		        	rm(list=c("runname","model.run.dir","FileList","TarList"))
			}

			ssh_disconnect(session)


#____________________________________________________________________________________________________________
# pull down finished models, untar, and clean directories
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	dir.list = dir( paste0(dir.model_runs,jobs.group),full.names = FALSE,recursive=FALSE,no..=TRUE)
	if("LLprof" %in% dir.list){dir.list=dir.list[-which(dir.list=="LLprof")]}

	for(i in 1:length(dir.list))
	{	

		model.run.dir = paste0(dir.model_runs,jobs.group,dir.list[i],"/")
		condor.dir = paste0(dir.launch,jobs.group,dir.list[i],"/")

		# download files
		scp_download(session, paste0(condor.dir,"End.tar.gz"), to = model.run.dir, verbose = FALSE)

		# untar
		shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -xzf End.tar.gz"),translate=TRUE)

		# delete and clean-up workspace
		shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -xzf End.tar.gz"),translate=TRUE)
		file.remove(paste0(model.run.dir,c("Start.tar.gz","End.tar.gz")))
		rm(list=c("model.run.dir","condor.dir"))			
	}

	ssh_disconnect(session)

#____________________________________________________________________________________________________________
# parse the output
	model.dir.list = paste0(dir.model_runs,jobs.group,dir.list,"/")

	# data.structures
	test.dt = data.table() %>%
						   .[,model_name:=dir.list] %>%
						   .[,model:=as.character(NA)] %>%
						   .[,endYear:=as.character(NA)] %>%
						   .[,type:=as.numeric(NA)] %>%
						   .[,model_complete:=0] %>%
						   .[,runtime_min:=as.numeric(NA)] %>%
						   .[,dep_latest:=as.numeric(NA)] %>%
						   .[,b_init:=as.numeric(NA)] %>%
						   .[,b_latest:=as.numeric(NA)] %>%
						   .[,sb_init:=as.numeric(NA)] %>%
						   .[,sb_latest:=as.numeric(NA)] %>%
						   .[,f_latest:=as.numeric(NA)] %>%
						   .[,Fmsy:=as.numeric(NA)] %>%
						   .[,F_Fmsy:=as.numeric(NA)] %>%
						   .[,MSY:=as.numeric(NA)] %>%
						   .[,SBmsy:=as.numeric(NA)] %>%
						   .[,SB_SBmsy:=as.numeric(NA)] %>%
						   .[,mgc:=as.numeric(NA)] %>%
						   .[,hessian:=NA] %>%
						   .[,tot_lik:=as.numeric(NA)] %>%
						   .[,bh_lik:=as.numeric(NA)] %>%
						   .[,edev_lik:=as.numeric(NA)] %>%
						   .[,qdev_lik:=as.numeric(NA)] %>%
						   .[,length_lik:=as.numeric(NA)] %>%
						   .[,weight_lik:=as.numeric(NA)] %>%
						   .[,catch_lik:=as.numeric(NA)] %>%
						   .[,npar:=as.numeric(NA)] %>%
						   .[,SBlatest_SB0:=as.numeric(NA)] %>%
						   .[,dep_inst_latest:=as.numeric(NA)] %>%
						   .[,dep_inst_initial:=as.numeric(NA)] %>%
						   .[,avg_annual_delta_rec:=as.numeric(NA)] %>%
						   .[,rel_avg_rec_rate:=as.numeric(NA)] %>%
						   .[,avg_reg_rec_prop:=as.numeric(NA)] %>%
						   .[,n_fish:=as.numeric(NA)] %>%
						   .[,bound_SDA:=as.numeric(NA)] %>%
						   .[,bound_SDB:=as.numeric(NA)] %>%
						   .[,bound_edev:=as.numeric(NA)] %>%
						   .[,SDA:=as.numeric(NA)] %>%
						   .[,SDB:=as.numeric(NA)]


	bio.ts = matrix(NA,nrow=length(model.dir.list),ncol=68)
	dep.ts = matrix(NA,nrow=length(model.dir.list),ncol=68)
	rownames(bio.ts) = rownames(dep.ts) = dir.list
	colnames(bio.ts) = colnames(dep.ts) = seq(from=1952,by=1,length.out=68)

	runtime.mat = matrix(NA,nrow=length(model.dir.list),ncol=18)

	# iterate over conventional models
	for(i in 1:length(model.dir.list))
	{	
		test.dt$model[i] = as.character(strsplit(dir.list[i],"_")[[1]][1])
		test.dt$endYear[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][2])
		test.dt$type[i] = as.character(strsplit(dir.list[i],"_")[[1]][3])
		tmp.runtime = scan(paste0(model.dir.list[i],"runtime.txt"))
		test.dt$runtime_min[i] = sum(tmp.runtime)/60
		
		# read
		if(length(list.files(model.dir.list[i]))==20)
		{
			tmp.rep = read.MFCLRep(paste0(model.dir.list[i],"plot-16.par.rep"))
			tmp.par = read.MFCLPar(paste0(model.dir.list[i],"16.par"))
			tmp.lik = read.MFCLLikelihood(paste0(model.dir.list[i],"test_plot_output"))
			tmp.hes = scan(paste0(model.dir.list[i],"neigenvalues"))
			
			# tmp_rmse = calc_index_RMSE(model.dir.list[i],par_stem="16",frq_stem = "swo",index_fishery = 14:19,nrmse_group=c(rep(1,13),rep(2,6)))
			# tmp_fit_length = parse_fit_file(model.dir.list[i],type="length",frq_stem = "swo",quantile_probs = c(0.25,0.5,0.75))
			# tmp_fit_weight = parse_fit_file(model.dir.list[i],type="weight",frq_stem = "swo",quantile_probs = c(0.25,0.5,0.75))

			# extract & store
			test.dt$model_complete[i] = 1
			test.dt$dep_latest[i] = tail(as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=10, sbf0_lag_nyears=1)),n=1)
			test.dt$b_init[i] = as.data.table(totalBiomass(tmp.rep))[,.(biomass=sum(value)),by=year]$biomass[1]
			test.dt$b_latest[i] = tail(as.data.table(totalBiomass(tmp.rep))[,.(biomass=sum(value)),by=year]$biomass,n=1)
			test.dt$sb_init[i] = head(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)),n=1)
			test.dt$sb_latest[i] = tail(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)),n=1)
			test.dt$f_latest[i] = mean(head(tail(as.vector(AggregateF(tmp.rep)),n=5),n=4))
			test.dt$Fmsy[i] = FMSY(tmp.rep)
			test.dt$F_Fmsy[i] = mean(head(tail(as.vector(AggregateF(tmp.rep)),n=5),n=4))/FMSY(tmp.rep)
			test.dt$MSY[i] =  MSY(tmp.rep)
			test.dt$SBmsy[i] = BMSY(tmp.rep)
			test.dt$SB_SBmsy[i] = ABBMSY(tmp.rep)
			test.dt$mgc[i] = max_grad(tmp.par)
			test.dt$hessian[i] = tmp.hes[1]
			test.dt$tot_lik[i] = obj_fun(tmp.par) # from par
			test.dt$bh_lik[i] = sum(bh_steep_contrib(tmp.lik))
			test.dt$edev_lik[i] = sum(effort_dev_penalty(tmp.lik))
			test.dt$qdev_lik[i] = sum(q_dev_pen_fish_grp(tmp.lik))
			test.dt$length_lik[i] = sum(total_length_fish(tmp.lik))
			test.dt$weight_lik[i] = sum(total_weight_fish(tmp.lik))
			test.dt$catch_lik[i] = sum(unlist(catch_fish(tmp.lik)))
			test.dt$npar[i] = n_pars(tmp.par)
			test.dt$SBlatest_SB0[i] = tail(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)),n=1)/head(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)),n=1)
			test.dt$dep_inst_latest[i] =  tail(as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=1, sbf0_lag_nyears=0)),n=1)
			test.dt$dep_inst_initial[i] =  head(as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=1, sbf0_lag_nyears=0)),n=1)
			test.dt$avg_annual_delta_rec[i] = lm(rec~yy,data=as.data.table(popN(tmp.rep))[age=="1"][,.(rec=sum(value),yy=as.numeric(year)),by=year])$coefficients["yy"]
			test.dt$rel_avg_rec_rate[i] = test.dt$avg_annual_delta_rec[i]/mean(as.data.table(popN(tmp.rep))[age=="1"][,.(rec=sum(value),yy=as.numeric(year)),by=year]$rec)
			tmp.dt=as.data.table(popN(tmp.rep))[age=="1"][,prop:=value/sum(value),by=year]
			test.dt$avg_reg_rec_prop[i] = mean(tmp.dt[area=="1"]$prop)
			test.dt$n_fish[i] = length(unique(as.data.table(q_fishery(tmp.rep))$unit))
			test.dt$SDA[i] = growth_var_pars(tmp.par)[1,1]
			test.dt$SDB[i] = growth_var_pars(tmp.par)[2,1]
			test.dt$bound_SDA[i] =  as.numeric(growth_var_pars(tmp.par)[1,1] < (growth_var_pars(tmp.par)[1,2]+diff(growth_var_pars(tmp.par)[1,2:3])*0.005) | growth_var_pars(tmp.par)[1,1] > (growth_var_pars(tmp.par)[1,3]-diff(growth_var_pars(tmp.par)[1,2:3])*0.005))
			test.dt$bound_SDB[i] =  as.numeric(growth_var_pars(tmp.par)[2,1] < (growth_var_pars(tmp.par)[2,2]+diff(growth_var_pars(tmp.par)[2,2:3])*0.005) | growth_var_pars(tmp.par)[2,1] > (growth_var_pars(tmp.par)[2,3]-diff(growth_var_pars(tmp.par)[2,2:3])*0.005))

			tmp.par = readLines(paste0(model.dir.list[i],"/16.par"))
			n.fish = grep("# correlation in selectivity deviations",tmp.par,fixed=TRUE)-grep("# effort deviation coefficients",tmp.par,fixed=TRUE)-2
			tmp.age_flags = scan(paste0(model.dir.list[i],"/16.par"),comment.char = "#",skip=grep("# age flags",tmp.par,fixed=TRUE),nlines=1,quiet=TRUE)
			tmp.edevs = scan(paste0(model.dir.list[i],"/16.par"),comment.char = "#",skip=grep("# effort deviation coefficients",tmp.par,fixed=TRUE),nlines=n.fish,quiet=TRUE)
			tmp.edevs.nozero = tmp.edevs[-which(tmp.edevs==0)]
			test.dt$bound_edev[i] = length(which(tmp.edevs.nozero >= tmp.age_flags[35]-0.01*tmp.age_flags[35] | tmp.edevs.nozero <= -tmp.age_flags[35]+0.01*tmp.age_flags[35]))


			bio.ts[i,rev(seq(from=68,by=-1,length.out=length(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)))))] = as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0))
			dep.ts[i,rev(seq(from=68,by=-1,length.out=length(as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)))))] = as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=1, sbf0_lag_nyears=0))
		}
		


		runtime.mat[i,] = tmp.runtime
		# clean-up
		rm(list=c("tmp.rep","tmp.par","tmp.lik","tmp.hes","tmp.runtime","tmp.dt"))
	}

	runtime_minutes.mat = t(apply(runtime.mat,1,cumsum))/60
	rownames(runtime_minutes.mat) = dir.list
	colnames(runtime_minutes.mat) = 1:18

#____________________________________________________________________________________________________________
# save and plot the output
	if (! dir.exists(paste0(proj.dir,"SWO/Assessment/Model_Summaries/",jobs.group)))dir.create(paste0(proj.dir,"SWO/Assessment/Model_Summaries/",jobs.group),recursive=TRUE)
	if (! dir.exists(paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group)))dir.create(paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),recursive=TRUE)

	# quick model selection
	test.dt = test.dt %>% .[,delta_totlik := abs(max(tot_lik)-tot_lik)] %>% .[,delta_lflik := abs(min(length_lik)+length_lik)] %>% .[,delta_wflik := abs(min(weight_lik)+weight_lik)] %>% .[,aic:=2*npar - 2*tot_lik]

	fwrite(test.dt,file=paste0(proj.dir,"SWO/Assessment/Model_Summaries/",jobs.group,"run_summaries.csv"))
	
	source(paste0(proj.dir,"Utilities/turbo.r"))
	turbo_pal=function(selected.model.names, all.model.names=selected.model.names)
	{
		n=length(all.model.names)-1
		out = turbo_vec(n,start = 0.1, end = 0.8)
		out = c("black",out)[1:length(all.model.names)]
	    names(out) = all.model.names
	    out = out[selected.model.names]
	    return(out)
	}


	tmp.dep = as.data.table(dep.ts) %>% .[,model:=rownames(dep.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","depletion")) %>% na.omit(.) %>%  
				merge(.,test.dt[,.(model_name,endYear)],by.x="model",by.y="model_name") %>% 
				.[,year:=as.numeric(as.character(year))] %>% .[year<=endYear]
	tmp.dep2 = as.data.table(dep.ts) %>% .[,model:=rownames(dep.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","depletion")) %>% na.omit(.) %>%  
				merge(.,test.dt[,.(model_name,endYear)],by.x="model",by.y="model_name") %>% 
				.[,year:=as.numeric(as.character(year))] %>% .[year>endYear] %>% setnames(.,"depletion","forecast")
	tmp.dep = merge(tmp.dep[,.(model,year,depletion)],tmp.dep2[,.(model,year,forecast)],by=c("model","year"),all=TRUE)
	p = tmp.dep %>%			
				ggplot() +
	     		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=depletion, color=model),size=2) + 
	  			geom_point(aes(x=year, y=forecast, fill=model),size=3,color="black",shape=21) +
	  			expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(rev(unique(tmp.dep$model)))) +
	  			scale_fill_manual("Model", values=turbo_pal(rev(unique(tmp.dep$model))))
	ggsave(filename="model.dep.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	tmp.bio = as.data.table(bio.ts) %>% .[,model:=rownames(bio.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","bio")) %>% .[,bio:=bio/1000] %>% na.omit(.) %>%  
				merge(.,test.dt[,.(model_name,endYear)],by.x="model",by.y="model_name") %>% 
				.[,year:=as.numeric(as.character(year))] %>% .[year<=endYear]
	tmp.bio2 = as.data.table(bio.ts) %>% .[,model:=rownames(bio.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","bio")) %>% .[,bio:=bio/1000] %>% na.omit(.) %>%  
				merge(.,test.dt[,.(model_name,endYear)],by.x="model",by.y="model_name") %>% 
				.[,year:=as.numeric(as.character(year))] %>% .[year>endYear] %>% setnames(.,"bio","forecast")
	tmp.bio = merge(tmp.bio[,.(model,year,bio)],tmp.bio2[,.(model,year,forecast)],by=c("model","year"),all=TRUE)
	p = tmp.bio %>%			
				ggplot() +
	     		xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=bio, color=model),size=2) + 
	  			geom_point(aes(x=year, y=forecast, fill=model),size=3,color="black",shape=21) +
	  			expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(rev(unique(tmp.bio$model)))) +
	  			scale_fill_manual("Model", values=turbo_pal(rev(unique(tmp.bio$model))))
	ggsave(filename="model.bio.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

#____________________________________________________________________________________________________________
# package and plot CPUE


	cpue_dt.list = as.list(rep(NA,nrow(test.dt)))
	idx_fsh = 14:16
	fishery_names = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1","05_EU_1","06_Other_1","07_DW_2N","08_DW_2C","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU")

	for(i in 1:length(cpue_dt.list))
	{
		tmp.dir = paste0(proj.dir,"SWO/Assessment/Model_Runs/",jobs.group,test.dt$model_name[i],"/")
		tmp_dt = cbind(cateffpen(readfrq(paste0(tmp.dir,"swo.frq"))),unlist(effort_dev_coffs(read.MFCLPar(paste0(tmp.dir,"16.par")))))
		colnames(tmp_dt)[8] = "edev"
		cpue_dt.list[[i]] = as.data.table(tmp_dt) %>% .[,model:=test.dt$endYear[i]] %>% .[fishery %in% idx_fsh]
		rm(list=c("tmp.dir","tmp_dt"))
	}

	cpue_dt = rbindlist(cpue_dt.list) %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,Fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names[idx_fsh])] %>% 
					 .[,Quarter:=c(1,1,1,2,2,2,3,3,3,4,4,4)[month]] 
	hindcast_dt_1 = rbindlist(cpue_dt.list) %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,Fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names[idx_fsh])] %>% 
					 .[,Quarter:=c(1,1,1,2,2,2,3,3,3,4,4,4)[month]] %>% .[year == as.numeric(as.character(model)) + 1]
	hindcast_dt_3 = rbindlist(cpue_dt.list) %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,Fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names[idx_fsh])] %>% 
					 .[,Quarter:=c(1,1,1,2,2,2,3,3,3,4,4,4)[month]] %>% .[year == as.numeric(as.character(model)) + 3]

	hindcast_dt_5 = rbindlist(cpue_dt.list) %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,Fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names[idx_fsh])] %>% 
					 .[,Quarter:=c(1,1,1,2,2,2,3,3,3,4,4,4)[month]] %>% .[year == as.numeric(as.character(model)) + 5]
		
	obs.dt = cpue_dt %>% .[model==2019,.(Fishery,Quarter,year,ts,cpue_obs,l_se,u_se,penalty)] %>% unique(.)
				
			# plots 
				p = cpue_dt %>% .[year<=as.numeric(as.character(model)) + 1] %>% .[order(model)] %>%
				ggplot() +
	     		xlab("Year") + ylab("Standardized CPUE") + facet_grid(Fishery~Quarter,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(data=obs.dt,aes(x=ts,xend=ts, y=l_se,yend=u_se),color="gray65",size=0.75) + 
	  			geom_point(data=obs.dt, aes(x=ts, y=cpue_obs),size=3,fill="gray90",color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_line( aes(x=ts, y=cpue_pred,color=model),size=1.25) +
	  			geom_point(data=hindcast_dt_1, aes(x=ts, y=cpue_pred,fill=model),size=3,color="black",shape=21) + 
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(rev(test.dt$endYear))) +
	  			scale_fill_manual("Model", values=turbo_pal(rev(test.dt$endYear)))
				ggsave(filename="hindcast_1.fit2cpue.fishery.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

			# plots 
				p = cpue_dt %>% .[year<=as.numeric(as.character(model)) + 3] %>% .[order(model)] %>%
				ggplot() +
	     		xlab("Year") + ylab("Standardized CPUE") + facet_grid(Fishery~Quarter,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(data=obs.dt,aes(x=ts,xend=ts, y=l_se,yend=u_se),color="gray65",size=0.75) + 
	  			geom_point(data=obs.dt, aes(x=ts, y=cpue_obs),size=3,fill="gray90",color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_line( aes(x=ts, y=cpue_pred,color=model),size=1.25) +
	  			geom_point(data=hindcast_dt_3, aes(x=ts, y=cpue_pred,fill=model),size=3,color="black",shape=21) + 
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(rev(test.dt$endYear))) +
	  			scale_fill_manual("Model", values=turbo_pal(rev(test.dt$endYear)))
				ggsave(filename="hindcast_3.fit2cpue.fishery.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

			# plots 
				p = cpue_dt %>% .[year<=as.numeric(as.character(model)) + 5] %>% .[order(model)] %>%
				ggplot() +
	     		xlab("Year") + ylab("Standardized CPUE") + facet_grid(Fishery~Quarter,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(data=obs.dt,aes(x=ts,xend=ts, y=l_se,yend=u_se),color="gray65",size=0.75) + 
	  			geom_point(data=obs.dt, aes(x=ts, y=cpue_obs),size=3,fill="gray90",color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_line( aes(x=ts, y=cpue_pred,color=model),size=1.25) +
	  			geom_point(data=hindcast_dt_5, aes(x=ts, y=cpue_pred,fill=model),size=3,color="black",shape=21) + 
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(rev(test.dt$endYear))) +
	  			scale_fill_manual("Model", values=turbo_pal(rev(test.dt$endYear)))
				ggsave(filename="hindcast_5.fit2cpue.fishery.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

#____________________________________________________________________________________________________________
# MASE calcs
		mase_dt = as.data.table(expand.grid(Fishery=fishery_names[idx_fsh],h_window=1:5)) %>% .[,mase:=as.numeric(NA)]

		for(f in 1:length(idx_fsh))
		{
				for(h in 1:5)
				{
					peel_dt = rbindlist(cpue_dt.list) %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,Fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names[idx_fsh])] %>% 
					 .[year == as.numeric(as.character(model)) + h] %>% .[Fishery==fishery_names[idx_fsh[f]]] %>%
					 .[,.(Fishery,year,ts,cpue_pred)] %>% merge(.,obs.dt[,.(Fishery,year,ts,cpue_obs)],by=c("Fishery","year","ts")) %>% .[order(ts)]
					 
					 tmp_mase = yardstick::mase_vec(peel_dt$cpue_obs,peel_dt$cpue_pred, m = 4*h)

					 mase_dt = mase_dt %>% .[Fishery==fishery_names[idx_fsh[f]]&h_window==h,mase:=tmp_mase]
					 rm(list=c("peel_dt","tmp_mase"))
				}
		}
	
#  mase_dt
#       Fishery h_window      mase
#  1: 14_idx_AU        1 1.0403823
#  2: 15_idx_NZ        1 0.7697709
#  3: 16_idx_EU        1 1.2925303
#  4: 14_idx_AU        2 0.9860213
#  5: 15_idx_NZ        2 0.7489998
#  6: 16_idx_EU        2 1.2713509
#  7: 14_idx_AU        3 0.7466104
#  8: 15_idx_NZ        3 0.5847295
#  9: 16_idx_EU        3 1.2273702
# 10: 14_idx_AU        4 0.7820844
# 11: 15_idx_NZ        4 0.3882792
# 12: 16_idx_EU        4 0.9955504
# 13: 14_idx_AU        5 1.0502938
# 14: 15_idx_NZ        5 0.2560910
# 15: 16_idx_EU        5 0.9153885

		out.tab = round(as.matrix(as.data.table(mase_dt[,!"Fishery"])),digits=3)
		colnames(out.tab) = c("$h_{window}$","MASE")

		out.tab = cbind(rep(c("14 AU idx","15 NZ idx","16 EU idx"),5),out.tab)
		colnames(out.tab) =  c("Index fishery","$h_{window}$","MASE")

		per.xtab <- xtable::xtable(out.tab, align=c("l","c",rep("c", dim(out.tab)[2]-1)), caption=paste0("MASE from \\enquote{Model-free} hindcasting with 1--5 year prediction periods ($h_{window}$) for each index fishery in the 2021 diagnostic case model"), label="tab:diag_hindcast",digits=matrix(c(0,0,0,3),15,4,byrow=FALSE))
		xtable::print.xtable(per.xtab, file=paste0(proj.dir,"Reports/SWO_SA/Tables/tbl.diag_hindcast.tex"), include.rownames=TRUE, sanitize.text.function = function(x) x, size="normalsize", caption.placement="top")


#  mase_dt[,.(mean(mase)),by=h_window]
#    h_window        V1
# 1:        1 1.0342279
# 2:        2 1.0021240
# 3:        3 0.8529034
# 4:        4 0.7219713
# 5:        5 0.7405911
# > mase_dt[,.(mean(mase)),by=Fishery]
#      Fishery        V1
# 1: 14_idx_AU 0.9210784
# 2: 15_idx_NZ 0.5495741
# 3: 16_idx_EU 1.1404381

