

# Nicholas Ducharme-Barth
# 12/07/2021
# 2021_diagnostic_retrospective_w_LLprof
# combine retrospective with likelihood profile on terminal biomass and depletion

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
	# dir.stem = paste0(proj.dir,"SWO/Assessment/Model_Runs/2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp_selexExplore/v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1/")
	dir.stem = paste0(dir.model_runs,"2021_diagnostic_case","/")
	dir.frq = paste0(proj.dir,"SWO/Assessment/Data_Prep/frq_file/")
	dir.launch = "/home/nicholasd/"
	dir.mfcl.launch = "/home/nicholasd/mfcl/"
	submit.realm = "NOU"
	supress_mfcl_IO = TRUE

#________________________________________________________________________________________________________________________________________________________________________________________________________
# create array of testing runs
	testing.options = expand.grid(model = "v2080",endYear = seq(from=2019,to=2009,by=-1), stringsAsFactors = FALSE)

#____________________________________________________________________________________________________________
# launch condor runs
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}
		jobs.group = "2021_diagnostic_retrospective_w_LLprof/"

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

				# modify doitall, change time period for SRR calcs
					if(testing.options$endYear[i] < 2019)
					{
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						tmp.doitall = gsub("2 199 64       # B-H SRR calculation begins in first model period; 1952",paste0("2 199 ",length(1952:testing.options$endYear[i]),"       # B-H SRR calculation begins in first model period; 1952"),tmp.doitall,fixed=TRUE)
						if(testing.options$endYear[i] < 2008) # remove time-block selectivity if before defined in selblocks.dat
						{
							tmp.doitall = gsub("-4 71 1","-4 71 0",tmp.doitall)
						}
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
					}

				# modify frq, downweight data after end year
					if(testing.options$endYear[i] < 2019)
					{
						tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
					  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[year<=testing.options$endYear[i]]
						cateffpen(tmp.frqit) = as.data.frame(catch.dt)
						lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
						# update lf
							ln.dt = as.data.table(lnfrq(tmp.frqit))  %>% .[year<=testing.options$endYear[i]] 
							lnfrq(tmp.frqit) = as.data.frame(ln.dt)
						# update wf
							wt.dt = as.data.table(wtfrq(tmp.frqit))  %>% .[year<=testing.options$endYear[i]] 
							wtfrq(tmp.frqit) = as.data.frame(wt.dt)
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
						   .[,version:=as.character(NA)] %>%
						   .[,endYear	:=as.numeric(NA)] %>%
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
	rec.ts = matrix(NA,nrow=length(model.dir.list),ncol=68)
	rownames(rec.ts) = rownames(bio.ts) = rownames(dep.ts) = dir.list
	colnames(rec.ts) = colnames(bio.ts) = colnames(dep.ts) = seq(from=1952,by=1,length.out=68)

	runtime.mat = matrix(NA,nrow=length(model.dir.list),ncol=18)

	# iterate over conventional models
	for(i in 1:length(model.dir.list))
	{	
		test.dt$model[i] = as.character(strsplit(dir.list[i],"_")[[1]][1])
		test.dt$endYear[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][2])
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


			bio.ts[i,(1:which(seq(from=1952,by=1,length.out=68)==test.dt$endYear[i]))] = as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0))
			dep.ts[i,(1:which(seq(from=1952,by=1,length.out=68)==test.dt$endYear[i]))] = as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=1, sbf0_lag_nyears=0))
			rec.ts[i,(1:which(seq(from=1952,by=1,length.out=68)==test.dt$endYear[i]-2))] = as.data.table(popN(tmp.rep))[age=="1"][,.(rec=sum(value),yy=as.numeric(year)),by=year][year%in%1952:(test.dt$endYear[i]-2)]$rec
			runtime.mat[i,] = tmp.runtime
		}
		# clean-up
		rm(list=c("tmp.rep","tmp.par","tmp.lik","tmp.hes","tmp.runtime"))
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

	tmp.dep = as.data.table(dep.ts) %>% .[,model:=rownames(dep.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","depletion")) %>% na.omit(.) %>% .[,endYear:=sapply(model,function(x)strsplit(x,"_")[[1]][3])]

	p = tmp.dep  %>% .[,.(model,year,depletion)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
				.[,model:=factor(model,levels=test.dt$model_name)] %>%
				ggplot() +
	     		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=depletion, color=model),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(test.dt$model_name))
	ggsave(filename="model.dep.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# p = tmp.dep  %>% .[,.(model,year,depletion)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
	# 			.[PDH=="1"] %>% 
	# 			ggplot() +
	#      		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	#   			geom_line( aes(x=year, y=depletion, color=model,linetype=PDH),size=2) + expand_limits(y=0) +
	#   			theme_few(base_size = 20)
	# ggsave(filename="model.dep.lines.subset.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	#   			scale = 1.25, width = 16, height = 9, units = c("in"),
	#   			dpi = 300, limitsize = TRUE)

	# test.dt[hessian == 0 & aic <= -49892.4]
	# p = tmp.dep  %>% .[,.(model,year,depletion)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian,aic)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
	# 			.[PDH=="1" & aic <= -49892.4] %>% 
	# 			ggplot() +
	#      		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	#   			geom_line( aes(x=year, y=depletion, color=model,linetype=PDH),size=2) + expand_limits(y=0) +
	#   			theme_few(base_size = 20)
	# ggsave(filename="model.dep.lines.subset.v2.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	#   			scale = 1.25, width = 16, height = 9, units = c("in"),
	#   			dpi = 300, limitsize = TRUE)

	tmp.bio = as.data.table(bio.ts) %>% .[,model:=rownames(bio.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","bio"))  %>% na.omit(.) %>% .[,endYear:=sapply(model,function(x)strsplit(x,"_")[[1]][3])]

	p = tmp.bio  %>% .[,.(model,year,bio)] %>% .[,bio:=bio/1000] %>% .[,year:=as.numeric(as.character(year))] %>% merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
				.[,model:=factor(model,levels=test.dt$model_name)] %>%
				ggplot() +
	     		xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=bio, color=model),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(test.dt$model_name))
	ggsave(filename="model.bio.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# p = tmp.bio  %>% .[,.(model,year,bio)] %>%  .[,bio:=bio/1000] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
	# 			.[PDH=="1"] %>% 
	# 			ggplot() +
	#      		xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	#   			geom_line( aes(x=year, y=bio, color=model,linetype=PDH),size=2) + expand_limits(y=0) +
	#   			theme_few(base_size = 20)
	# ggsave(filename="model.bio.lines.subset.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	#   			scale = 1.25, width = 16, height = 9, units = c("in"),
	#   			dpi = 300, limitsize = TRUE)

	# p = tmp.bio  %>% .[,.(model,year,bio)] %>% .[,bio:=bio/1000] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian,aic)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
	# 			.[PDH=="1" & aic <= -49892.4] %>% 
	# 			ggplot() +
	#      		xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	#   			geom_line( aes(x=year, y=bio, color=model,linetype=PDH),size=2) + expand_limits(y=0) +
	#   			theme_few(base_size = 20)
	# ggsave(filename="model.bio.lines.subset.v2.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	#   			scale = 1.25, width = 16, height = 9, units = c("in"),
	#   			dpi = 300, limitsize = TRUE)

	tmp.rec = as.data.table(rec.ts) %>% .[,model:=rownames(rec.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","rec")) %>% na.omit(.) %>% .[,endYear:=sapply(model,function(x)strsplit(x,"_")[[1]][3])]

	p = tmp.rec  %>% .[,.(model,year,rec)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>% 
				.[,model:=factor(model,levels=test.dt$model_name)] %>%
				ggplot() +
	     		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=rec, color=model),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(test.dt$model_name))
	ggsave(filename="model.rec.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	tmp.time = as.data.table(runtime_minutes.mat) %>% .[,model:=rownames(runtime_minutes.mat)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("phase","time")) 

	p = tmp.time  %>% .[,.(model,phase,time)] %>% .[,phase:=as.numeric(as.character(phase))] %>% .[,pf50:= as.character(-as.numeric(lapply(model,function(x)strsplit(x,"_")[[1]][2])))] %>% .[,pf387:=lapply(model,function(x)strsplit(x,"_")[[1]][3])] %>%
				.[,model:=factor(model,levels=test.dt$model_name)] %>%
				ggplot() +
	     		xlab("Estimation phase") + ylab("Elapsed time (minutes)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=phase, y=time, color=model),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(test.dt$model_name))
	ggsave(filename="model.etime.lines.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
#___________________________________________________________________________________________________________
# plot selectivity...
	rep.list= lapply(paste0(model.dir.list,"plot-",16,".par.rep"),read.MFCLRep)
	names(rep.list) = sapply(model.dir.list,function(x)tail(strsplit(x,"/")[[1]],n=1))

	agg.years = TRUE
	agg.regions=TRUE
	biomass.type = "SSB"
	biomass.units=1000
	recdist.year_range=NULL
	yaxis.free=FALSE
	LRP=0.20
	TRP=NULL

	fsh.lab = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1 - early","04_AU_1 - late","05_EU_1","06_Other_1","07_DW_2N","08_DW_2C","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU")
	
	# selectivity - age
		p = diags4MFCL::plot.selectivity(rep.list,names(rep.list),sel.basis="AGE",palette.func=turbo_pal,fsh.lab=fsh.lab,all.model.names=sapply(model.dir.list,function(x)tail(strsplit(x,"/")[[1]],n=1)))
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("sel.age.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# selectivity - length
		p = diags4MFCL::plot.selectivity(rep.list,names(rep.list),sel.basis="LENGTH",palette.func=turbo_pal,fsh.lab=fsh.lab,all.model.names=sapply(model.dir.list,function(x)tail(strsplit(x,"/")[[1]],n=1)))
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("sel.length.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)


#___________________________________________________________________________________________________________
# plot fit to length & weight
	source(paste0(proj.dir,"Utilities/parse.fit_file.r"))
	for(i in 1:nrow(test.dt))
	{
		p = plot_fit_file(parse_fit_file(paste0(dir.model_runs,jobs.group,test.dt$model_name[i],"/"),type="length","swo"),type="length",fishery_names=fishery_names,save_dir=paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),save_name_stem=test.dt$model_name[i])
		p = plot_fit_file(parse_fit_file(paste0(dir.model_runs,jobs.group,test.dt$model_name[i],"/"),type="weight","swo"),type="weight",fishery_names=fishery_names,save_dir=paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),save_name_stem=test.dt$model_name[i])

	}


#___________________________________________________________________________________________________________
# CPUE; goodness of fit, RMSE


	   rc_pal = function(n,alpha=1)
	   {
	  	  # taken from 'https://github.com/jabbamodel/ss3diags/blob/master/R/SSplotModelcomp.R'
	   	  # their comment >>>
	      # a subset of rich.colors by Arni Magnusson from the gregmisc package
	      # a.k.a. rich.colors.short, but put directly in this function
	      # to try to diagnose problem with transparency on one computer
	      x = seq(0, 1, length = n)
	      r = 1/(1 + exp(20 - 35 * x))
	      g = pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
	      b = dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
	      rgb.m = matrix(c(r, g, b), ncol = 3)
	      rich.vector = apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha=alpha))
	    }

	cpue_rmse.model_dir_list = paste0(dir.model_runs,jobs.group,test.dt$model_name,"/")
	cpue_met.dt_list = as.list(rep(NA,length(cpue_rmse.model_dir_list)))

	for(m in 1:length(cpue_rmse.model_dir_list))
	{
			# if(m==1)
			# {
			# 	fishery_names = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1","05_Other_1","06_DW_2N","07_DW_2C_pre2001","08_DW_2C_post2001","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU")

			# } else {
			# 	fishery_names = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1","05_Other_1","06_DW_2N","07_DW_2C_pre2001","08_DW_2C_post2001","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU","17_idx_JPearly1","18_idx_JPearly2","19_idx_JPmid1","20_idx_JPmid2","21_idx_JPlate1","22_idx_JPlate2")
			# }
			tmp_dt = cbind(cateffpen(readfrq(paste0(cpue_rmse.model_dir_list[m],"swo.frq"))),unlist(effort_dev_coffs(read.MFCLPar(paste0(cpue_rmse.model_dir_list[m],"16.par")))))
			colnames(tmp_dt)[8] = "edev"
			tmp.dt = as.data.table(tmp_dt) %>% .[,mean_penalty := mean(penalty),by=fishery] %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[mean_penalty != 0.001 & !is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=fishery] %>% .[,norm_effort:=effort/mean(effort)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,residual := cpue_obs - cpue_pred] %>% .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,n:=.N,by=fishery] %>% .[,mean_cpue_obs := mean(cpue_obs),by=fishery] %>% .[,devsq_cpue_obs := sum((cpue_obs-mean_cpue_obs)^2),by=fishery] %>% .[,rse:=sd(residual),by=fishery] %>% .[,leverage:=1/n+(cpue_obs-mean_cpue_obs)^2/devsq_cpue_obs] %>%
					 .[,std_residual := residual/(rse*sqrt(1-leverage))] %>% .[,fishery:=factor(as.character(fishery),levels=as.character(1:length(unique(tmp_dt[,4]))),labels=fishery_names)] %>% .[,runs:=TSA::runs(residual)$pvalue,by=fishery] %>%
					 .[,max_u_se:=max(u_se,cpue_pred),by=fishery] %>% .[,max_res := max(residual),by=fishery] %>% .[,max_std_res  := max(std_residual),by=fishery]
		
			# metrics
				tmp_met.dt = tmp.dt %>% .[,.(mse=mean((cpue_obs - cpue_pred)^2),rmse=sqrt(mean((cpue_obs - cpue_pred)^2)),nrmse=sqrt(mean((cpue_obs - cpue_pred)^2))/mean(cpue_obs),mean_penalty=median(penalty),mean_cv=median(input_cv),runs=mean(runs),max_u_se=mean(max_u_se), max_res=mean(max_res),max_std_res=mean(max_std_res)),by=fishery] %>%
				 			            .[,runs_status := ifelse(runs <= 0.05,"Fail","Pass")] %>%
				 			            .[,x:=min(tmp.dt$ts)] %>% .[,x2:=max(tmp.dt$ts)] %>% .[,model:=tail(strsplit(cpue_rmse.model_dir_list[m],"/")[[1]],n=1)]

				
			# plots 
				p = tmp.dt %>%  
				ggplot() +
	     		xlab("Year") + ylab("Standardized CPUE") + facet_wrap(~fishery,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(aes(x=ts,xend=ts, y=l_se,yend=u_se,color=penalty),size=0.75) + 
	  			geom_point( aes(x=ts, y=cpue_obs, fill=penalty),size=3,color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_line( aes(x=ts, y=cpue_pred),color="black",size=1) +
	  			theme_few(base_size = 20) +
	  			scale_colour_gradientn("Penalty",colours = turbo_vec(100,start = 0.1, end = 0.8)) +
	  			scale_fill_gradientn("Penalty",colours = turbo_vec(100,start = 0.1, end = 0.8))	+
	  			geom_rect(data=tmp_met.dt,aes(xmin = x, xmax = 0.4*(x2-x)+x, ymin = max_u_se*0.85, ymax = max_u_se*1.05),fill="white",alpha=0.85) +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max_u_se,label=paste0("Median penalty weight: ",round(mean_penalty,digits=2))),hjust="inward") +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max_u_se*0.95,label=paste0("Median input CV: ",round(mean_cv,digits=3))),hjust="inward") +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max_u_se*0.90,label=paste0("CV (RMSE): ",round(nrmse,digits=3))),hjust="inward")
				ggsave(filename=paste0(tail(strsplit(cpue_rmse.model_dir_list[m],"/")[[1]],n=1),".fit2cpue.fishery.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


				p = tmp.dt %>%  
				ggplot() +
	     		xlab("Year") + ylab("Standardized Residual") + facet_wrap(~fishery,scales="free_y") +
	     		geom_rect(aes(xmin = min(ts), xmax = max(ts), ymin = -3, ymax = 3),fill="gray95") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(aes(x=ts,xend=ts, y=0,yend=std_residual)) + 
	  			geom_point( aes(x=ts, y=std_residual, fill=penalty),size=3,color="black",shape=21) + expand_limits(y=c(-3.5,3.5)) +
	  			geom_smooth( aes(x=ts, y=std_residual),color="black",size=1.5,se=FALSE) +
	  			theme_few(base_size = 20) +
	  			scale_fill_gradientn("Penalty",colours = turbo_vec(100,start = 0.1, end = 0.8))	+	 
				geom_rect(data=tmp_met.dt,aes(xmin = x, xmax = 0.5*(x2-x)+x, ymin = max(c(max_std_res,3.5))*0.7, ymax = max(c(max_std_res,3.5))*1.05),fill="white",alpha=0.85) +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max(c(max_std_res,3.5)),label=paste0("Median input CV: ",round(mean_cv,digits=3))),hjust="inward") +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max(c(max_std_res,3.5))*0.9,label=paste0("CV (RMSE): ",round(nrmse,digits=3))),hjust="inward") +
	  			geom_text(data=tmp_met.dt,aes(x=x,y=max(c(max_std_res,3.5))*0.8,label=paste0("Runs test: ",runs_status,", p-value: ",round(runs,digits=3))),hjust="inward")
				ggsave(filename=paste0(tail(strsplit(cpue_rmse.model_dir_list[m],"/")[[1]],n=1),".stdres.fishery.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


				tmp2.dt = tmp.dt %>% .[,.(q25=quantile(std_residual,probs=0.25),q50=median(std_residual),q75=quantile(std_residual,probs=0.75),n=.N),by=ts] %>% .[n>1]
				tmp3.dt = tmp.dt %>% .[,.(min=min(std_residual),q50=median(std_residual),max=max(std_residual),n=.N),by=ts] %>% .[n==1,min:=ifelse(min<0,min,0)] %>% .[n==1,max:=ifelse(max>0,max,0)]
				p = tmp.dt %>%  
				ggplot() +
	     		xlab("Year") + ylab("Standardized Residual") + 
	     		geom_rect(aes(xmin = min(ts), xmax = max(ts), ymin = -3, ymax = 3),fill="gray95") +
			    geom_segment(data=tmp3.dt,aes(x=ts,xend=ts, y=min,yend=max),color="black") + 
	     		geom_rect(data=tmp2.dt,aes(xmin = ts-(1.5/12), xmax = ts+(1.5/12), ymin = q25, ymax = q75),fill="gray70",color="black") +
	     		geom_segment(data=tmp2.dt,aes(x = ts-(1.5/12), xend = ts+(1.5/12), y = q50, yend = q50),,color="black") + 
	     		geom_hline(yintercept = 0) +
	  			geom_point( aes(x=ts, y=std_residual, fill=fishery),size=3,color="black",shape=21) + expand_limits(y=c(-3.5,3.5)) +
	  			geom_smooth( aes(x=ts, y=std_residual),color="black",size=1.5,se=FALSE) +
	  			theme_few(base_size = 20) +
	  			scale_fill_manual("Fishery",values = unname(sapply(rc_pal(length(unique(tmp.dt$fishery))+1)[-1],function(x)substr(x,1,7)))) +
	  			scale_color_manual("Fishery",values = unname(sapply(rc_pal(length(unique(tmp.dt$fishery))+1)[-1],function(x)substr(x,1,7))))
				ggsave(filename=paste0(tail(strsplit(cpue_rmse.model_dir_list[m],"/")[[1]],n=1),".stdres.combined.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

				p = tmp.dt %>%  
				ggplot() +
	     		xlab("Year") + ylab("Effort deviate") + facet_wrap(~fishery,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(aes(x=ts,xend=ts, y=0,yend=edev)) + 
	  			geom_point( aes(x=ts, y=edev, fill=penalty),size=3,color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_smooth( aes(x=ts, y=edev),color="black",size=1.5,se=FALSE) +
	  			theme_few(base_size = 20) +
	  			scale_fill_gradientn("Penalty",colours = turbo_vec(100,start = 0.1, end = 0.8))		 
				ggsave(filename=paste0(tail(strsplit(cpue_rmse.model_dir_list[m],"/")[[1]],n=1),".edev.fishery.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			cpue_met.dt_list[[m]] = tmp_met.dt[,.(model,fishery,mse,rmse,nrmse,mean_penalty,mean_cv,runs)]
	  			rm(list=c("p","tmp.dt","tmp2.dt","tmp3.dt","tmp_met.dt","tmp.frq","pointer","tmp.cateffpen","len.cateffpen","tmp_dt"))
	} 

	cpue_met.dt = rbindlist(cpue_met.dt_list)

#___________________________________________________________________________________________________________
# visualize recruitment partitioning between regions
	
	rec_dt.list = as.list(rep(NA,length(model.dir.list)))
	relrec_dt.list = as.list(rep(NA,length(model.dir.list)))
	regpar_dt.list = as.list(rep(NA,length(model.dir.list)))

	# iterate over conventional models
	for(i in 1:length(model.dir.list))
	{	
		# read
		tmp.rep = read.MFCLRep(paste0(model.dir.list[i],"plot-16.par.rep"))
		tmp.endYear = strsplit(model.dir.list[i],"/")
		rec_dt.list[[i]] = as.data.table(popN(tmp.rep)) %>% .[age=="1",.(year,area,value)] %>% setnames(.,"value","rec") %>% .[,year:=as.numeric(year)] %>% 
						  .[,model:=tail(strsplit(model.dir.list[i],"/")[[1]],n=1)] %>% .[,prop:=rec/sum(rec),by=year] %>% .[,.(model,year,area,rec,prop)]
		tmp.par = readLines(paste0(model.dir.list[i],"16.par"))
		pointer = grep("# relative recruitment",tmp.par) + 2
		relrec_dt.list[[i]] = data.table(model=tail(strsplit(model.dir.list[i],"/")[[1]],n=1),year=(1952:as.numeric(tail(strsplit(tail(strsplit(model.dir.list[i],"/")[[1]],n=1),"_")[[1]],n=1)))[-1],relrec = as.numeric(strsplit(trimws(tmp.par[pointer]),"\\s+")[[1]])) # account for lag between spawners and recruits
		tmp.mat = matrix(scan(paste0(model.dir.list[i],"16.par"),skip=grep("# regional recruitment variation",tmp.par),nlines=length(1952:as.numeric(tail(strsplit(tail(strsplit(model.dir.list[i],"/")[[1]],n=1),"_")[[1]],n=1)))),ncol=2,byrow=TRUE)
		colnames(tmp.mat) = 1:2
		tmp.mat = as.data.table(tmp.mat) %>% .[,year:=1952:as.numeric(tail(strsplit(tail(strsplit(model.dir.list[i],"/")[[1]],n=1),"_")[[1]],n=1))] %>% melt(.,id.vars = "year") %>% setnames(.,c("variable","value"),c("area","reg_rec_dev")) %>% .[,model:=tail(strsplit(model.dir.list[i],"/")[[1]],n=1)]
		rec_dt.list[[i]] = merge(rec_dt.list[[i]],tmp.mat)
		pointer = grep("# region parameters",tmp.par) + 1
		regpar_dt.list[[i]] = data.table(model=tail(strsplit(model.dir.list[i],"/")[[1]],n=1),area=as.character(1:2),avg_prop = as.numeric(strsplit(trimws(tmp.par[pointer]),"\\s+")[[1]]))

		rm(list=c("tmp.rep","tmp.par","pointer","tmp.mat"))
	}

	rec_dt = rbindlist(rec_dt.list)
	relrec_dt = rbindlist(relrec_dt.list)
	regpar_dt = rbindlist(regpar_dt.list)

	# srr
		p = diags4MFCL::plot.srr(rep.list, show.legend=TRUE, palette.func=turbo_pal)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("srr.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# total relative rec
		p = relrec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
		ggplot() + 
	    xlab("Year") + ylab("Total relative recruitment") + 
	    geom_hline(yintercept = c(0,1)) +
	  	geom_line( aes(x=year, y=relrec, color=model,linetype=PDH),size=1.15) + expand_limits(y=c(0)) +
	  	geom_smooth( aes(x=year, y=relrec,color=model),size=1.5,se=FALSE) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_pal(unique(relrec_dt$model)))
	  	ggsave(filename="total_relrec.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# regional rec devs
	  	p = rec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
		ggplot() + facet_wrap(~model) +
	    xlab("Year") + ylab("Regional recruitment deviate") + 
	    geom_hline(yintercept = 0) +
	  	geom_line( aes(x=year, y=reg_rec_dev, color=area,linetype=PDH),size=1.15) + expand_limits(y=c(0)) +
	  	geom_smooth( aes(x=year, y=reg_rec_dev,color=area),size=1.25,se=FALSE,alpha=0.5) +
	  	geom_smooth( aes(x=year, y=reg_rec_dev),color="black",size=1.5,se=FALSE) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8))
	  	ggsave(filename="reg_rec_devs.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


	# total reg rec
		p = rec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
		ggplot() + facet_wrap(~model) +
	    xlab("Year") + ylab("Recruitment (millions)") + 
	    geom_hline(yintercept = 0) +
	  	geom_line( aes(x=year, y=rec/1000000, color=area,linetype=PDH),size=1.15) + expand_limits(y=c(0)) +
	  	geom_smooth( aes(x=year, y=rec/1000000),color="black",size=1.5,se=FALSE) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8))
	  	ggsave(filename="total_rec.reg.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


	# total reg rec ribbon
		p = rec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
		.[area=="2",rec:=-rec] %>% .[area=="2",ymin:=rec] %>% .[area=="2",ymax:=0] %>% .[area=="1",ymin:=0] %>% .[area=="1",ymax:=rec] %>%
		ggplot() + facet_wrap(~model,scales="free_y") +
	    xlab("Year") + ylab("Recruitment (millions)") +
	   	geom_ribbon(aes(x=year, ymin=ymin/1000000,ymax=ymax/1000000, fill=area),alpha=0.5) +
	    geom_hline(yintercept = 0) +
	  	geom_line( aes(x=year, y=rec/1000000, color=area,linetype=PDH),size=1.15) + expand_limits(y=c(0)) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8)) +
	  	scale_fill_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8))
	  	ggsave(filename="total_rec.reg.ribbon.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# mean scaled reg rec
	  	p = rec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
	  		.[,rec:=rec/mean(rec),by=.(model,area)] %>%
		ggplot() + facet_wrap(~model) +
	    xlab("Year") + ylab("Relative recruitment") + 
	    geom_hline(yintercept = c(0,1)) +
	  	geom_line( aes(x=year, y=rec, color=area,linetype=PDH),size=1.15) + expand_limits(y=c(0)) +
	  	geom_smooth( aes(x=year, y=rec,color=area),size=1.25,se=FALSE,alpha=0.5) +
	  	geom_smooth( aes(x=year, y=rec),color="black",size=1.5,se=FALSE) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8))
	  	ggsave(filename="scale_rec.reg.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# prop ribbon rec
	  	tmp_regpar = copy(regpar_dt); tmp_regpar = tmp_regpar %>% .[area=="2",avg_prop:=-avg_prop]
	  	p = rec_dt %>%  merge(.,test.dt[,.(model_name,hessian)],by.x="model",by.y="model_name") %>% .[,PDH:=factor(as.character(ifelse(hessian==0,1,0)),levels=c("1","0"))] %>%
	  		.[area=="2",prop:=-prop] %>% .[area=="2",ymin:=prop] %>% .[area=="2",ymax:=0] %>% .[area=="1",ymin:=0] %>% .[area=="1",ymax:=prop] %>%
		ggplot() + facet_wrap(~model) +
	    xlab("Year") + ylab("Relative recruitment") + 
	    geom_ribbon(aes(x=year, ymin=ymin,ymax=ymax, fill=area),alpha=0.5) +
	    geom_hline(yintercept = c(0)) +
	    geom_hline(data=tmp_regpar,aes(yintercept = avg_prop,color=area),size=1) +
	  	geom_line( aes(x=year, y=prop, color=area,linetype=PDH),size=1.15) + expand_limits(y=c(-1,1)) +
	  	theme_few(base_size = 20) +
		scale_linetype_manual("PDH",values = c("solid","dashed")) +
	  	scale_color_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8)) +
	  	scale_fill_manual("Region",values = turbo_vec(2,start = 0.1, end = 0.8))
	  	ggsave(filename="prop_rec.reg.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

#____________________________________________________________________________________________________________
# calculate Mohn's rho
# for recruitment, biomass, and depletion
	  	rec.rho = bio.rho = dep.rho = rep(NA,(nrow(dep.ts)-1))
	  	for(j in 1:(nrow(dep.ts)-1))
	  	{
	  		dep.rho[j] = (dep.ts[nrow(dep.ts)-(j),ncol(dep.ts)-(j)] - dep.ts[nrow(dep.ts),ncol(dep.ts)-(j)])/dep.ts[nrow(dep.ts),ncol(dep.ts)-(j)]
	  		bio.rho[j] = (bio.ts[nrow(bio.ts)-(j),ncol(bio.ts)-(j)] - bio.ts[nrow(bio.ts),ncol(bio.ts)-(j)])/bio.ts[nrow(bio.ts),ncol(bio.ts)-(j)]
	  		rec.rho[j] = (rec.ts[nrow(rec.ts)-(j),ncol(rec.ts)-(j+2)] - rec.ts[nrow(rec.ts),ncol(rec.ts)-(j+2)])/rec.ts[nrow(rec.ts),ncol(rec.ts)-(j+2)]
	  	}
	  	mean(dep.rho) # 0.08156153
	  	mean(bio.rho) # 0.3439664
	  	mean(rec.rho) # 0.08651056


#____________________________________________________________________________________________________________
# launch likelihood profile
	if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	# define profile intervals
	down.int = seq(from=99,by=-1,length.out=60) # c(99, 95, 91, 87, 83, 79, 75, 71, 67, 63, 59, 55, 51)
	up.int = seq(from=101,by=1,length.out=60) # c(101, 105, 109, 113, 117, 121, 125, 129, 133, 137, 141, 145, 149)
	profile.mgc.flag = -4
	llprof.options = expand.grid(direction=c("up","down"),profile_type=c("biomass","depletion"),category=c("adult","total"),time_period=c("average","last_year"),model=apply(testing.options,1,paste0,collapse="_"),stringsAsFactors=FALSE)
	llprof.options = subset(llprof.options, !(profile_type=="depletion" & time_period %in% c("average","initial")))
	rownames(llprof.options) = 1:nrow(llprof.options)

	# iterate across llprof.options
			ssh_exec_wait(session, command = paste0("mkdir ",dir.launch,jobs.group,"LLprof/"))
			
			for(i in 1:nrow(llprof.options))
			{
				runname = paste0(llprof.options[i,c(5,2,3,4,1)],collapse="/")
				model.run.dir = paste0(dir.model_runs,jobs.group,"LLprof/",runname,"/")

				# get largest .par				
				final.par = paste0(max(as.numeric(sapply(list.files(paste0(dir.model_runs,jobs.group,llprof.options[i,5],"/"))[grep(".par$",list.files(paste0(dir.model_runs,jobs.group,llprof.options[i,5],"/")))],function(x)strsplit(x,"[.]")[[1]][1]))),".par")
				final.rep = read.MFCLRep(paste0(paste0(dir.model_runs,jobs.group,llprof.options[i,5],"/"),"plot-",final.par,".rep"))

				# create new directory for model run
					if (! dir.exists(model.run.dir))dir.create(model.run.dir,recursive=TRUE)
					ssh_exec_wait(session, command = paste0("mkdir -p ",dir.launch,jobs.group,"LLprof/",runname))
				# transfer files
					FileList=c("mfcl.cfg","selblocks.dat","swo.frq",final.par)
					file.copy( paste0(dir.model_runs,jobs.group,llprof.options[i,5],"/",FileList),model.run.dir,overwrite=TRUE)

				# copy condor files
				if(llprof.options$profile_type[i] == "biomass")
				{
					if(llprof.options$direction[i] == "down")
					{
						file.copy(paste0(dir.condor,"LLprof/",c("condor.sub","mfcl.X86_64.LINUX.bat","Lprofile.Bio1.sh")),model.run.dir,overwrite=TRUE)
						# rename file
						file.rename(paste0(model.run.dir,c("Lprofile.Bio1.sh")),paste0(model.run.dir,c("doitall.condor")))
					} else {
						file.copy(paste0(dir.condor,"LLprof/",c("condor.sub","mfcl.X86_64.LINUX.bat","Lprofile.Bio2.sh")),model.run.dir,overwrite=TRUE)
						# rename file
						file.rename(paste0(model.run.dir,c("Lprofile.Bio2.sh")),paste0(model.run.dir,c("doitall.condor")))
					}
				} else {
					if(llprof.options$direction[i] == "down")
					{
						file.copy(paste0(dir.condor,"LLprof/",c("condor.sub","mfcl.X86_64.LINUX.bat","Lprofile.Dep1.sh")),model.run.dir,overwrite=TRUE)
						# rename file
						file.rename(paste0(model.run.dir,c("Lprofile.Dep1.sh")),paste0(model.run.dir,c("doitall.condor")))
					} else {
						file.copy(paste0(dir.condor,"LLprof/",c("condor.sub","mfcl.X86_64.LINUX.bat","Lprofile.Dep2.sh")),model.run.dir,overwrite=TRUE)
						# rename file
						file.rename(paste0(model.run.dir,c("Lprofile.Dep2.sh")),paste0(model.run.dir,c("doitall.condor")))
					}

					# modify .bat file
					tmp.bat = readLines(paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"),warn=FALSE)
					tmp.bat = gsub("avg_bio","relative_depletion",tmp.bat)
					writeLines(tmp.bat,con=paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"))
					rm(list=c("tmp.bat"))
				}


				# modify doitall.condor
					tmp.doitall = readLines(paste0(model.run.dir,"doitall.condor"),warn=FALSE)
					tmp.doitall = gsub("AAA.frq","swo.frq",tmp.doitall)
					tmp.doitall = gsub("BBB.par",final.par,tmp.doitall)
					pointer = grep("for Mult in",tmp.doitall)
					if(llprof.options$direction[i] == "down")
					{
						tmp.doitall[pointer] = paste0("for Mult in ",paste0(down.int,collapse = " "))
					} else {
						tmp.doitall[pointer] = paste0("for Mult in ",paste0(up.int,collapse = " "))
					}
					if(llprof.options$category[i] == "adult")
					{
						tmp.doitall = gsub("2 172 0","2 172 1",tmp.doitall)
					}
					if(llprof.options$time_period[i] == "initial")
					{
						tmp.doitall = gsub("Af173=0",paste0("Af173=",unname(dimensions(final.rep)['years'])),tmp.doitall)
						tmp.doitall = gsub("Af174=0",paste0("Af174=",unname(dimensions(final.rep)['years'])-4),tmp.doitall)
					} else if(llprof.options$time_period[i] == "average"){
						tmp.doitall = gsub("Af173=0",paste0("Af173=",unname(dimensions(final.rep)['years'])),tmp.doitall)
						tmp.doitall = gsub("Af174=0","Af174=1",tmp.doitall)
					} else if(llprof.options$time_period[i] == "recent"){
						tmp.doitall = gsub("Af173=0","Af173=5",tmp.doitall)
						tmp.doitall = gsub("Af174=0","Af174=1",tmp.doitall)
					} else {
						tmp.doitall = gsub("Af173=0","Af173=1",tmp.doitall)
						tmp.doitall = gsub("Af174=0","Af174=1",tmp.doitall)
					}
					tmp.doitall = gsub("1 50 -2",paste0("1 50 ",profile.mgc.flag),tmp.doitall)
					writeLines(tmp.doitall,con=paste0(model.run.dir,"doitall.condor"))
					rm(list=c("tmp.doitall","pointer"))

					if(submit.realm != "SUV")
					{
						if(submit.realm == "NOU")
						{
							tmp.sub = readLines(paste0(model.run.dir,"condor.sub"),warn=FALSE)
							pointer = grep("Requirements",tmp.sub)
							tmp.sub[pointer] = 'Requirements = ((Realm == "NOU") && (OpSys == "LINUX") && (Machine =!= "NOUOFPCALC02.corp.spc.int") && (Machine =!= "NOUOFPCAND01") && (Machine =!= "NOUOFPCAND02") && (Machine =!= "NOUOFPCAND03") && (Machine =!= "NOUOFPCAND04") && (Machine =!= "NOUOFPCAND05") && (Machine =!= "NOUOFPCAND06") && (Machine =!= "NOUOFPCAND07") && (Machine =!= "NOUOFPCAND08") && (Machine =!= "nouofpcand09.corp.spc.int") && (Machine =!= "nouofpcand10") && (Machine =!= "nouofpcand11.corp.spc.int") && (Machine =!= "nouofpcand12"))'
							writeLines(tmp.sub,con=paste0(model.run.dir,"condor.sub"))
							rm(list=c("tmp.sub","pointer"))
						} else {
							tmp.sub = readLines(paste0(model.run.dir,"condor.sub"),warn=FALSE)
							pointer = grep("Requirements",tmp.sub)
							tmp.sub[pointer] = 'Requirements = ((OpSys == "LINUX") && (Machine =!= "NOUOFPCALC02.corp.spc.int") && (Machine =!= "NOUOFPCAND01") && (Machine =!= "NOUOFPCAND02") && (Machine =!= "NOUOFPCAND03") && (Machine =!= "NOUOFPCAND04") && (Machine =!= "NOUOFPCAND05") && (Machine =!= "NOUOFPCAND06") && (Machine =!= "NOUOFPCAND07") && (Machine =!= "NOUOFPCAND08") && (Machine =!= "nouofpcand09.corp.spc.int") && (Machine =!= "nouofpcand10") && (Machine =!= "nouofpcand11.corp.spc.int") && (Machine =!= "nouofpcand12"))'
							writeLines(tmp.sub,con=paste0(model.run.dir,"condor.sub"))
							rm(list=c("tmp.sub","pointer"))
						}
					}

					tmp.bat = readLines(paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"),warn=FALSE)
					if(supress_mfcl_IO)
					{
						tmp.bat = gsub("./doitall.condor","./doitall.condor &>/dev/null",tmp.bat,fixed=TRUE)
					}
					writeLines(tmp.bat,con=paste0(model.run.dir,"mfcl.X86_64.LINUX.bat"))

				# make tar with the necessary files
					TarList=c("doitall.condor","mfcl.cfg","selblocks.dat","swo.frq",final.par)
					shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -czf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
				# send all files to the launch machine
			       	scp_upload(session,files=paste0(model.run.dir,c("Start.tar.gz","condor.sub","mfcl.X86_64.LINUX.bat")),to=paste0(dir.launch,jobs.group,"LLprof/",runname))
				# copy mfcl to launch machine run directory
			       	ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'v2080_dev_20210503/mfclo64 ',paste0(dir.launch,jobs.group,"LLprof/",runname)))
			    # remake tar file to include mfcl
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"LLprof/",runname),'; tar -xzf Start.tar.gz; rm Start.tar.gz'))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"LLprof/",runname),'; tar -czf Start.tar.gz ',paste(c(TarList,"mfclo64"),collapse=' ')))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"LLprof/",runname),'; rm ',paste(c(TarList,"mfclo64"),collapse=' ')))
				# launch condor job
			        ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"LLprof/",runname),'; dos2unix mfcl.X86_64.LINUX.bat; chmod 700 mfcl.X86_64.LINUX.bat; condor_submit condor.sub'))   
	
				# clean-up
		        	rm(list=c("runname","model.run.dir","FileList","TarList","final.par","final.rep"))
			}

			ssh_disconnect(session)

#____________________________________________________________________________________________________________
# pull down finished likelihood profiles, untar, and clean directories
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	dir.list = unname(apply(llprof.options[,c(5,2,3,4,1)],1,function(x)paste0(x,collapse="/")))

	for(i in 1:length(dir.list))
	{	

		model.run.dir = paste0(dir.model_runs,jobs.group,"LLprof/",dir.list[i],"/")
		condor.dir = paste0(dir.launch,jobs.group,"LLprof/",dir.list[i],"/")

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
# make plots of likelihood profile by year of retrospective peel
	retroprofile_dt.list = as.list(rep(NA,nrow(llprof.options)))

	A = proc.time()
	for(i in 1:nrow(llprof.options))
	{
		if(llprof.options$direction[i]=="down")
		{
			tmp.int = down.int
		} else {
			tmp.int = c(100,up.int)
		}

		if(llprof.options$profile_type[i]=="biomass")
		{
			tmp.file = "avg_bio"
		} else {
			tmp.file = "relative_depletion"
		}

		tmp.list = as.list(rep(NA,length(tmp.int)))
		for(j in 1:length(tmp.list))
		{
			tmp.out = readLines(paste0(dir.model_runs,jobs.group,"LLprof/",paste0(llprof.options[i,c(5,2,3,4,1)],collapse="/"),"/","test_plot_output_",tmp.int[j]),warn=FALSE)
			tmp.par = readLines(paste0(dir.model_runs,jobs.group,"LLprof/",paste0(llprof.options[i,c(5,2,3,4,1)],collapse="/"),"/","bio",tmp.int[j],"final.par"),warn=FALSE)
			n.fish = length(strsplit(trimws(tmp.out[grep("# Effort_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
			tmp.dt = as.data.frame(matrix(0,nrow=n.fish,ncol=9))
			colnames(tmp.dt) = c("fishery","total","bh","catchability","catch","effort","length","tag","weight")
			tmp.dt = as.data.table(tmp.dt)

			# fishery invariant
			# BH
				tmp.dt$bh = rep(as.numeric(trimws(tmp.out[grep("# BH_steep contribution",tmp.out,fixed=TRUE)+1]))/n.fish,n.fish)
			# tag
				if(length(grep("# tag release", tmp.out))>0)
				{
					tmp.dt$tag = rep(sum(scan(paste0(tmp.dir,"test_plot_output_",tmp.mult[j]),comment.char="#",skip=grep("# tag release",tmp.out,fixed=TRUE)[1],nlines=grep("# Tag likelihood for pooled tag release by fishery groups",tmp.out,fixed=TRUE) - grep("# tag release",tmp.out,fixed=TRUE)[1]))/n.fish,n.fish)
				}
			# fishery specific
			# effort dev
				tmp.dt$effort = as.numeric(strsplit(trimws(tmp.out[grep("# Effort_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
			# catchability
				tmp.dt$catchability = as.numeric(strsplit(trimws(tmp.out[grep("# catchability_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
			# length
				tmp.dt$length = as.numeric(strsplit(trimws(tmp.out[grep("# total length component of likelihood for each fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
			# weight
				tmp.dt$weight = as.numeric(strsplit(trimws(tmp.out[grep("# total weight component of likelihood for each fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
			# catch
				tmp.dt$catch = sapply(strsplit(tmp.out[grep("total catch components of likelihood for fishery", tmp.out)+1], split="[[:blank:]]+"), function(x)sum(as.numeric(x))) 

			# total
				tmp.dt$total = rep(-as.numeric(tmp.par[grep("# Objective function value",tmp.par)+1])/n.fish,n.fish)
					
				tmp.dt$fishery = 1:n.fish
				tmp.dt$grand_total = sum(tmp.dt$total)
				tmp.dt$grand_total_par = sum(tmp.dt$total)
				tmp.dt$other = rep((sum(tmp.dt$total) - (sum(tmp.dt$bh) + sum(tmp.dt$tag) + sum(tmp.dt$effort) + sum(tmp.dt$catchability) + sum(tmp.dt$length) + sum(tmp.dt$weight) + sum(tmp.dt$catch)))/n.fish,n.fish)

				tmp.dt$grand_total_bh = sum(tmp.dt$bh)
				tmp.dt$grand_total_tag = sum(tmp.dt$tag)
				tmp.dt$grand_total_effort = sum(tmp.dt$effort)
				tmp.dt$grand_total_catchability = sum(tmp.dt$catchability)
				tmp.dt$grand_total_length = sum(tmp.dt$length)
				tmp.dt$grand_total_weight = sum(tmp.dt$weight)
				tmp.dt$grand_total_catch = sum(tmp.dt$catch)
				tmp.dt$grand_total_other = sum(tmp.dt$other)
				tmp.dt$mult = tmp.int[j]
				tmp.dt$profile_type = llprof.options$profile_type[i]
				tmp.dt$category = llprof.options$category[i]
				tmp.dt$time_period = llprof.options$time_period[i]
				tmp.dt$endYear = as.numeric(strsplit(llprof.options$model[i],"_")[[1]][2])
				tmp.dt$avg_value = scan(paste0(dir.model_runs,jobs.group,"LLprof/",paste0(llprof.options[i,c(5,2,3,4,1)],collapse="/"),"/",tmp.file))
				tmp.dt$true_value = (tmp.dt$mult)/100 * tmp.dt$avg_value
				tmp.list[[j]] = tmp.dt
				rm(list=c("tmp.dt","tmp.out","tmp.par","n.fish"))
		}

			retroprofile_dt.list[[i]] = rbindlist(tmp.list)		
			rm(list=c("tmp.int","tmp.list"))
	}
	B = proc.time()
	B-A
	retroprofile_dt = rbindlist(retroprofile_dt.list)
	fwrite(retroprofile_dt,file=paste0(paste0(proj.dir,"SWO/Assessment/Model_Summaries/",jobs.group,"retroprofile_dt.csv")))

		p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "last_year"] %>% unique(.,by=c("profile_type","category","time_period","endYear","grand_total_par")) %>%
								.[,.(endYear,true_value,grand_total_par,grand_total_bh,grand_total_effort,grand_total_length,grand_total_weight,grand_total_other)] %>% 
		 						.[,total_dNLL:=abs(grand_total_par-min(grand_total_par)),by=.(endYear)] %>%
		 						.[,bh_dNLL:=abs(grand_total_bh-min(grand_total_bh)),by=.(endYear)] %>%
		 						.[,effort_dNLL:=abs(grand_total_effort-min(grand_total_effort)),by=.(endYear)] %>%
		 						.[,length_dNLL:=abs(grand_total_length-min(grand_total_length)),by=.(endYear)] %>%
		 						.[,weight_dNLL:=abs(grand_total_weight-min(grand_total_weight)),by=.(endYear)] %>%
		 						.[,other_dNLL:=abs(grand_total_other-min(grand_total_other)),by=.(endYear)] %>%
		 						.[,.(endYear,true_value,total_dNLL,bh_dNLL,effort_dNLL,length_dNLL,weight_dNLL,other_dNLL)] %>%
		 						 .[order(endYear,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL","other_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight","Other"))] %>%
		ggplot() + 
	    xlab("Terminal spawning potential (1000's mt)") + ylab("Change in total likelihood") + 
	    facet_wrap(~variable) +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_terminal.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

		p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "average"] %>% unique(.,by=c("profile_type","category","time_period","endYear","grand_total_par")) %>%
								.[,.(endYear,true_value,grand_total_par,grand_total_bh,grand_total_effort,grand_total_length,grand_total_weight,grand_total_other)] %>% 
		 						.[,total_dNLL:=abs(grand_total_par-min(grand_total_par)),by=.(endYear)] %>%
		 						.[,bh_dNLL:=abs(grand_total_bh-min(grand_total_bh)),by=.(endYear)] %>%
		 						.[,effort_dNLL:=abs(grand_total_effort-min(grand_total_effort)),by=.(endYear)] %>%
		 						.[,length_dNLL:=abs(grand_total_length-min(grand_total_length)),by=.(endYear)] %>%
		 						.[,weight_dNLL:=abs(grand_total_weight-min(grand_total_weight)),by=.(endYear)] %>%
		 						.[,other_dNLL:=abs(grand_total_other-min(grand_total_other)),by=.(endYear)] %>%
		 						.[,.(endYear,true_value,total_dNLL,bh_dNLL,effort_dNLL,length_dNLL,weight_dNLL,other_dNLL)] %>%
		 						 .[order(endYear,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL","other_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight","Other"))] %>%
		ggplot() + 
	    xlab("Average spawning potential (1000's mt)") + ylab("Change in total likelihood") + 
	    facet_wrap(~variable) +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_average.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)


		p = retroprofile_dt %>% .[profile_type == "depletion" & category == "adult" & time_period == "last_year"] %>% unique(.,by=c("profile_type","category","time_period","endYear","grand_total_par")) %>%
								.[,.(endYear,true_value,grand_total_par,grand_total_bh,grand_total_effort,grand_total_length,grand_total_weight,grand_total_other)] %>% 
		 						.[,total_dNLL:=abs(grand_total_par-min(grand_total_par)),by=.(endYear)] %>%
		 						.[,bh_dNLL:=abs(grand_total_bh-min(grand_total_bh)),by=.(endYear)] %>%
		 						.[,effort_dNLL:=abs(grand_total_effort-min(grand_total_effort)),by=.(endYear)] %>%
		 						.[,length_dNLL:=abs(grand_total_length-min(grand_total_length)),by=.(endYear)] %>%
		 						.[,weight_dNLL:=abs(grand_total_weight-min(grand_total_weight)),by=.(endYear)] %>%
		 						.[,other_dNLL:=abs(grand_total_other-min(grand_total_other)),by=.(endYear)] %>%
		 						.[,.(endYear,true_value,total_dNLL,bh_dNLL,effort_dNLL,length_dNLL,weight_dNLL,other_dNLL)] %>%
		 						 .[order(endYear,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL","other_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight","Other"))] %>%
		ggplot() + 
	    xlab("Depletion of spawning biomass relative to static initial spawning biomass") + ylab("Change in total likelihood") + 
	    facet_wrap(~variable) +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.dep_adult_terminal.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "last_year"] %>% 
								.[,.(endYear,fishery,true_value,total)] %>% 
		 						.[,total_dNLL:=abs(total-min(total)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,total_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
		 						 .[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Terminal spawning potential (1000's mt)") + ylab("Change in total likelihood (Total)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_terminal.fsh_total.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "last_year" & fishery %in% 14:18] %>% 
								.[,.(endYear,fishery,true_value,effort)] %>% 
		 						.[,effort_dNLL:=abs(effort-min(effort)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,effort_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Terminal spawning potential (1000's mt)") + ylab("Change in total likelihood (CPUE)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_terminal.fsh_cpue.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "last_year" & length!=0] %>% 
								.[,.(endYear,fishery,true_value,length)] %>% 
		 						.[,length_dNLL:=abs(length-min(length)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,length_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Terminal spawning potential (1000's mt)") + ylab("Change in total likelihood (Length)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_terminal.fsh_length.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "last_year" & weight!=0] %>% 
								.[,.(endYear,fishery,true_value,weight)] %>% 
		 						.[,weight_dNLL:=abs(weight-min(weight)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,weight_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Terminal spawning potential (1000's mt)") + ylab("Change in total likelihood (Weight)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_terminal.fsh_weight.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)



	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "average"] %>% 
								.[,.(endYear,fishery,true_value,total)] %>% 
		 						.[,total_dNLL:=abs(total-min(total)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,total_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
		 						 .[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Average spawning potential (1000's mt)") + ylab("Change in total likelihood (Total)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_average.fsh_total.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "average" & fishery %in% 14:18] %>% 
								.[,.(endYear,fishery,true_value,effort)] %>% 
		 						.[,effort_dNLL:=abs(effort-min(effort)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,effort_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Average spawning potential (1000's mt)") + ylab("Change in total likelihood (CPUE)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_average.fsh_cpue.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "average" & length!=0] %>% 
								.[,.(endYear,fishery,true_value,length)] %>% 
		 						.[,length_dNLL:=abs(length-min(length)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,length_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Average spawning potential (1000's mt)") + ylab("Change in total likelihood (Length)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_average.fsh_length.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = retroprofile_dt %>% .[profile_type == "biomass" & category == "adult" & time_period == "average" & weight!=0] %>% 
								.[,.(endYear,fishery,true_value,weight)] %>% 
		 						.[,weight_dNLL:=abs(weight-min(weight)),by=.(endYear,fishery)] %>%
		 						.[,.(endYear,true_value,fishery,weight_dNLL)] %>%
		 						 .[order(endYear,fishery,true_value)] %>% .[,endYear:=as.character(endYear)] %>% melt(.,id.vars=c("endYear","true_value","fishery")) %>% .[,variable:=factor(variable,levels=c("total_dNLL","bh_dNLL","effort_dNLL","length_dNLL","weight_dNLL"),labels=c("Total","BH-SRR","CPUE","Length","Weight"))] %>%
								.[,fishery:=factor(as.character(fishery),levels=as.character(1:length(fishery_names)),labels=fishery_names)] %>%
		ggplot() + 
	    xlab("Average spawning potential (1000's mt)") + ylab("Change in total likelihood (Weight)") + 
	    facet_wrap(~fishery,scale="free_y") +
	    geom_line(aes(x=true_value/1000,y=value,group=endYear,color=endYear),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_color_manual("Terminal year",values = turbo_vec(length(unique(retroprofile_dt$endYear)),start = 0.1, end = 0.8))
	  	ggsave(filename="retro_likelihood.bio_adult_average.fsh_weight.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

