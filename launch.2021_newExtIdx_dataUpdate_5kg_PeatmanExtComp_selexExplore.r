

# Nicholas Ducharme-Barth
# 07/07/2021
# 2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp_selexExplore
# Use "2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp/v2080_6_1_none_0_0_0_0_0_1_18" as a jumping off point to explore selectivity changes which produce a better fit
# Selex block for fishery 12
# Group EU fishery; possible selex block
# remove non-decreasing for fishery 3
# remove zero selex at age 1 for fishery 2

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
	dir.stem = paste0(proj.dir,"SWO/Assessment/Model_Runs/2017_DiagCase_2080_test_idxFsh/v2080_6_1_1_0_1_0/")
	dir.frq = paste0(proj.dir,"SWO/Assessment/Data_Prep/frq_file/")
	frq.prefix = "2021_frq_5kg_newExtractComp"
	dir.launch = "/home/nicholasd/"
	dir.mfcl.launch = "/home/nicholasd/mfcl/"
	submit.realm = "NOU"
	supress_mfcl_IO = TRUE

#________________________________________________________________________________________________________________________________________________________________________________________________________
# create array of testing runs
	testing.options = expand.grid(model="v2080",pf50=c(6),no_JP=c(1),TW_cpue=c("none"),JP_share_regsel = c(0),TW_share_regsel=c(0),shift_growth=c(0),start1994=c(0),increase_weighting_dwfn_comp=c(0),est_SDB=c(1),ff3=c(12,14,16,18),group_f5_f11=c(0,1),block_f11=c(0,2011),block_f12=c(0,2010),rm_nodec_f3=c(0,1),rm_0slex_f2=c(0,1),rm_0slex_f8=c(0,1),stringsAsFactors=FALSE)
	testing.options = subset(testing.options,!(TW_cpue == "none"&TW_share_regsel==1))
	testing.options = unique(testing.options)
	rownames(testing.options) = 1:nrow(testing.options)
	testing.options = subset(testing.options,!(no_JP == 1&JP_share_regsel==1))
	testing.options = unique(testing.options)
	rownames(testing.options) = 1:nrow(testing.options)
	testing.options = subset(testing.options,!(no_JP == 1&TW_cpue == "none"&increase_weighting_dwfn_comp == 1))
	testing.options = unique(testing.options)
	rownames(testing.options) = 1:nrow(testing.options)

#________________________________________________________________________________________________________________________________________________________________________________________________________
# define biology - use code from test_2sex_v2 to get the average biology, this will be called 'combined'
	load(paste0(proj.dir,"SWO/Assessment/Data_Prep/model_ensemble/ensemble_output/post.samp.RData"))

		# define growth options
			tmp_a.vec = 0:19
			tmp_laa_mean = colMeans(t(apply(post.samp,1,function(x)(x[['Linf_F']]*(1-exp(-x[['k_F']]*(tmp_a.vec-x[['t0']])))+x[['Linf_M']]*(1-exp(-x[['k_M']]*(tmp_a.vec-x[['t0']]))))/2)))
			L1_new = tmp_laa_mean[1]
			L2_new = tmp_laa_mean[20]

#____________________________________________________________________________________________________________
# launch condor runs
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}
		jobs.group = "2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp_selexExplore/"

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
				# modify doitall, convergence tolerance
					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
					pointer = grep("# modify convergence criteria from this phase",tmp.doitall,fixed=TRUE)
					tmp.doitall[pointer] = paste0("  1 50 -",testing.options$pf50[i],"        # modify convergence criteria from this phase")
					writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
					rm(list=c("tmp.doitall","pointer"))

				# modify frq, TW_cpue
					if(testing.options$TW_cpue[i] %in% c("none"))
					{
						if(testing.options$TW_cpue[i] == "none")
						{
							idx_fsh = 14:22
							jp_idx_fsh = 17:22
							tw_idx_fsh = 0

							tmp.frqit = readfrq(paste0(dir.frq,"input_frq/",frq.prefix,".swo.frq"))
							# update header information
								n_fisheries(tmp.frqit) = max(idx_fsh)
								region_fish(tmp.frqit) = FLQuant(c(c(aperm(region_fish(tmp.frqit),c(3,1,2,4,5,6)))[1:n_fisheries(tmp.frqit)]),dim=c(1,1,n_fisheries(tmp.frqit),1,1,1))
							  	data_flags(tmp.frqit) = data_flags(tmp.frqit)[,1:n_fisheries(tmp.frqit)]
							# update cateffpen
							  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[fishery <= n_fisheries(tmp.frqit)]
								cateffpen(tmp.frqit) = as.data.frame(catch.dt)
								lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
								valid_fshinst = apply(cateffpen(tmp.frqit)[,c("year","month","fishery")],1,paste0,collapse="_")
							# update lf
								ln.dt = as.data.table(lnfrq(tmp.frqit)) %>% .[fishery <= n_fisheries(tmp.frqit)] %>% .[,id:=paste0(year,"_",month,"_",fishery)] %>% data.table::setcolorder(., neworder=c("id","year","month","week","fishery",as.character(seq(from=lf_range(tmp.frqit)['LFFirst'],by=lf_range(tmp.frqit)['LFWidth'],length.out=lf_range(tmp.frqit)['LFIntervals'])))) %>%
										.[id %in% valid_fshinst]
								ln.dt = ln.dt[, !"id"]
								lnfrq(tmp.frqit) = as.data.frame(ln.dt)

							writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
	  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
	  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
	  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
	  						rm(list=c("catch.dt","ln.dt","valid_fshinst"))
						}
					} else {
						idx_fsh = 14:26
						jp_idx_fsh = 17:22
						tw_idx_fsh = 23:26
						tmp.frqit = readfrq(paste0(dir.frq,"input_frq/",frq.prefix,".swo.frq"))
						writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
					}

					if(testing.options$no_JP[i] == 1)
					{
						if(length(tw_idx_fsh) > 1)
						{
							idx_fsh = c(14:16,seq(from=17,by=1,length.out=length(tw_idx_fsh)))
							old_tw_idx_fsh = tw_idx_fsh
							tw_idx_fsh = seq(from=17,by=1,length.out=length(tw_idx_fsh))
							tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
								# update header information
									n_fisheries(tmp.frqit) = max(idx_fsh)
									region_fish(tmp.frqit) = FLQuant(c(c(aperm(region_fish(tmp.frqit),c(3,1,2,4,5,6)))[1:n_fisheries(tmp.frqit)]),dim=c(1,1,n_fisheries(tmp.frqit),1,1,1))
								  	data_flags(tmp.frqit) = data_flags(tmp.frqit)[,1:n_fisheries(tmp.frqit)]
								# update cateffpen
								  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[!(fishery %in% jp_idx_fsh)] %>% .[fishery>16,fishery:=tw_idx_fsh[match(fishery,old_tw_idx_fsh)]]
									cateffpen(tmp.frqit) = as.data.frame(catch.dt)
									lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
								# update lf
									ln.dt = as.data.table(lnfrq(tmp.frqit))  %>% .[!(fishery %in% jp_idx_fsh)] %>% .[fishery>16,fishery:=tw_idx_fsh[match(fishery,old_tw_idx_fsh)]]
									lnfrq(tmp.frqit) = as.data.frame(ln.dt)

								writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
		  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
		  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
		  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
		  						rm(list=c("catch.dt","ln.dt"))

							jp_idx_fsh = 0
						} else {
							idx_fsh = 14:16
							jp_idx_fsh = 0
							tw_idx_fsh = 0

							tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
								# update header information
									n_fisheries(tmp.frqit) = max(idx_fsh)
									region_fish(tmp.frqit) = FLQuant(c(c(aperm(region_fish(tmp.frqit),c(3,1,2,4,5,6)))[1:n_fisheries(tmp.frqit)]),dim=c(1,1,n_fisheries(tmp.frqit),1,1,1))
								  	data_flags(tmp.frqit) = data_flags(tmp.frqit)[,1:n_fisheries(tmp.frqit)]
								# update cateffpen
								  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[fishery %in% 1:n_fisheries(tmp.frqit)] 
									cateffpen(tmp.frqit) = as.data.frame(catch.dt)
									lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
								# update lf
									ln.dt = as.data.table(lnfrq(tmp.frqit))  %>% .[fishery %in% 1:n_fisheries(tmp.frqit)] 
									lnfrq(tmp.frqit) = as.data.frame(ln.dt)

								writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
		  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
		  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
		  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
		  						rm(list=c("catch.dt","ln.dt"))
						}
					}

					rm(list=c("tmp.frqit","tmp.frq"))

				# modify frq, start year 1994
				if(testing.options$start1994[i]==1)
				{
					if(length(jp_idx_fsh)>1)
					{
						jp_idx_fsh = 17:18
						if(length(tw_idx_fsh)>1)
						{
							tw_idx_fsh = seq(from=19,by=1,length.out=length(tw_idx_fsh))
							idx_fsh = 14:max(tw_idx_fsh)
						} else {
							idx_fsh = 14:max(jp_idx_fsh)
						}

						tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
								# update header information
									n_fisheries(tmp.frqit) = max(idx_fsh)
									region_fish(tmp.frqit) = FLQuant(c(c(aperm(region_fish(tmp.frqit),c(3,1,2,4,5,6)))[1:n_fisheries(tmp.frqit)]),dim=c(1,1,n_fisheries(tmp.frqit),1,1,1))
								  	data_flags(tmp.frqit) = data_flags(tmp.frqit)[,1:n_fisheries(tmp.frqit)]
								  	range(tmp.frqit)['minyear'] = 1952
								# update cateffpen
								  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[!(fishery%in% 17:20)] %>% .[!(fishery > 13 & year < 1994)] %>% .[fishery>16,fishery:=fishery-4] 
									cateffpen(tmp.frqit) = as.data.frame(catch.dt)
									lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
								# update lf
									ln.dt = as.data.table(lnfrq(tmp.frqit))  %>% .[!(fishery%in% 17:20)] %>% .[!(fishery > 13 & year < 1994)] %>%.[fishery>16,fishery:=fishery-4] 
									lnfrq(tmp.frqit) = as.data.frame(ln.dt)
								# update wf
									wt.dt = as.data.table(wtfrq(tmp.frqit))  %>% .[!(fishery%in% 17:20)] %>% .[!(fishery > 13 & year < 1994)] %>%.[fishery>16,fishery:=fishery-4] 
									wtfrq(tmp.frqit) = as.data.frame(wt.dt)

								writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
		  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
		  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
		  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
		  						rm(list=c("catch.dt","ln.dt","wt.dt"))
					} else {
						if(length(tw_idx_fsh)>1)
						{
							tw_idx_fsh = seq(from=17,by=1,length.out=length(tw_idx_fsh))
							idx_fsh = 14:max(tw_idx_fsh)
						} 

							tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
								# update header information
									n_fisheries(tmp.frqit) = max(idx_fsh)
									region_fish(tmp.frqit) = FLQuant(c(c(aperm(region_fish(tmp.frqit),c(3,1,2,4,5,6)))[1:n_fisheries(tmp.frqit)]),dim=c(1,1,n_fisheries(tmp.frqit),1,1,1))
								  	data_flags(tmp.frqit) = data_flags(tmp.frqit)[,1:n_fisheries(tmp.frqit)]
								  	range(tmp.frqit)['minyear'] = 1952
								# update cateffpen
								  	catch.dt = as.data.table(cateffpen(tmp.frqit)) %>% .[!(fishery > 13 & year < 1994)]  
									cateffpen(tmp.frqit) = as.data.frame(catch.dt)
									lf_range(tmp.frqit)['Datasets'] = nrow(cateffpen(tmp.frqit))
								# update lf
									ln.dt = as.data.table(lnfrq(tmp.frqit))  %>% .[!(fishery > 13 & year < 1994)] 
									lnfrq(tmp.frqit) = as.data.frame(ln.dt)
								# update wf
									wt.dt = as.data.table(wtfrq(tmp.frqit))  %>% .[!(fishery > 13 & year < 1994)] 
									wtfrq(tmp.frqit) = as.data.frame(wt.dt)

								writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
		  						tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
		  						tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
		  						writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
		  						rm(list=c("catch.dt","ln.dt","wt.dt"))
					}
				}

				# modify doitall, set-up index structure based on idx_fsh
					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
					# Phase 1, cubic spline selex
						pointer = grep("-9 61 5       # except for fishery 9 which has 5 nodes",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",idx_fsh," 61 5"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 1, non-decreasing selex (default is non-decreasing)
						pointer = grep("-9 16 1",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",idx_fsh," 16 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 1, re-order fisheries with selectivity block (AU index is now 14)
						tmp.doitall = gsub("-15 71 1","-14 71 1",tmp.doitall)
					# Phase 1, set selex to zero for JP & TW idx
						pointer = grep("-9 75 6",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:3))]
						# if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==0)
						# {
						# 	tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",jp_idx_fsh," 75 1"),paste0("    -",tw_idx_fsh," 75 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						# }
						# if(testing.options$TW_cpue[i]=="none"&testing.options$no_JP[i]==0){
						# 	tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",jp_idx_fsh," 75 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						# }
						# if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==1){
						# 	tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",tw_idx_fsh," 75 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						# }
					# Phase 1, grouping of fisheries with shared selex (default is JP all together and TW all together)
						pointer = grep("-13 24 12",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						if(testing.options$TW_cpue[i]=="full" & testing.options$no_JP[i]==0)
						{
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",14:16," 24 ",13:15),paste0("    -",jp_idx_fsh," 24 16"),paste0("    -",tw_idx_fsh," 24 17"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						} else if (testing.options$TW_cpue[i]=="none" & testing.options$no_JP[i]==0){
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",14:16," 24 ",13:15),paste0("    -",jp_idx_fsh," 24 16"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						} else if (testing.options$TW_cpue[i]=="full" & testing.options$no_JP[i]==1){
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",14:16," 24 ",13:15),paste0("    -",tw_idx_fsh," 24 16"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						} else {
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",14:16," 24 ",13:15),tmp.doitall[(pointer+1):length(tmp.doitall)])
						}
					# Phase 1, catchability grouping start
						pointer = grep("-13 60 13",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",idx_fsh," 60 ",idx_fsh),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 1, catchability grouping general
						pointer = grep("-13 29 13",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",idx_fsh," 29 ",idx_fsh),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 2, increase weight on catch likelihood
						pointer = grep("-999 45 100000    # Increase weight on catch likelihood",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("   -",idx_fsh," 45 10000"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 3, increase weight on catch likelihood
						pointer = grep("2 110 50        # set penalty for recruitment deviates",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("   -",idx_fsh," 45 50000"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 4, increase weight on catch likelihood
						pointer = grep("2 68 0          # do not estimate movement",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("   -",idx_fsh," 45 100000"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					# Phase 5, seasonal q
						pointer = grep("-999 27 0       # estimate seasonal catchability for all fisheries",tmp.doitall,fixed=TRUE)
						# remove old index fleet structure
						tmp.doitall = tmp.doitall[-c(pointer+(1:5))]
						tmp.doitall = c(tmp.doitall[1:pointer],paste0("   -",idx_fsh," 27 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
					writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
					rm(list=c("tmp.doitall","pointer"))

				# modify doitall, JP share selex by region
					if(testing.options$JP_share_regsel[i] == 1)
					{
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						if(testing.options$TW_cpue[i] == "none")
						{
							jp_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-17 24",tmp.doitall)])," ")[[1]][3])
							for(j in 1:length(jp_idx_fsh))
							{
								if(rev(jp_idx_fsh %% 2)[j] == 1)
								{
									tmp.doitall[grep(paste0("-",jp_idx_fsh[j]," 24"),tmp.doitall)] = paste0("    -",jp_idx_fsh[j]," 24 ",jp_selex_group+1)
								}
							}
						} else {
							jp_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-17 24",tmp.doitall)])," ")[[1]][3])
							for(j in 1:length(jp_idx_fsh))
							{
								if(rev(jp_idx_fsh %% 2)[j] == 1)
								{
									tmp.doitall[grep(paste0("-",jp_idx_fsh[j]," 24"),tmp.doitall)] = paste0("    -",jp_idx_fsh[j]," 24 ",jp_selex_group+1)
								}
							}
							tw_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-22 24",tmp.doitall)])," ")[[1]][3]) +1
							for(j in 1:length(tw_idx_fsh))
							{
								tmp.doitall[grep(paste0("-",tw_idx_fsh[j]," 24"),tmp.doitall)] = paste0("    -",tw_idx_fsh[j]," 24 ",tw_selex_group)
							}
						}
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
					}

				# modify doitall, TW share selex by region
					if(testing.options$TW_share_regsel[i] == 1)
					{
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						if(testing.options$no_JP[i] == 1)
						{
							tw_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-16 24",tmp.doitall)])," ")[[1]][3]) +1
						} else {
							if(testing.options$start1994[i]==1)
							{
								tw_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-18 24",tmp.doitall)])," ")[[1]][3]) +1
							} else {
								tw_selex_group = as.numeric(strsplit(trimws(tmp.doitall[grep("-22 24",tmp.doitall)])," ")[[1]][3]) +1
							}
						}
						for(j in 1:length(tw_idx_fsh))
						{
								if(rev(tw_idx_fsh %% 2)[j] == 1)
								{
									tmp.doitall[grep(paste0("-",tw_idx_fsh[j]," 24"),tmp.doitall)] = paste0("    -",tw_idx_fsh[j]," 24 ",tw_selex_group+1)
								}
						}
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
					}

				# modify doitall, start year 1994
					if(testing.options$start1994[i]==1)
					{
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						tmp.doitall = gsub("2 199 64       # B-H SRR calculation begins in first model period; 1952",paste0("2 199 ",length(1994:2019),"       # B-H SRR calculation begins in first model period; ",1994),tmp.doitall,fixed=TRUE)
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
					} else {
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						tmp.doitall = gsub("2 199 64       # B-H SRR calculation begins in first model period; 1952",paste0("2 199 ",length(1952:2019),"       # B-H SRR calculation begins in first model period; 1952"),tmp.doitall,fixed=TRUE)
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
					}

				# modify doitall & selblocks.dat
						tmp.selblocks = readLines(paste0(model.run.dir,"selblocks.dat"),warn=FALSE)
						writeLines(tmp.selblocks[-1],con=paste0(model.run.dir,"selblocks.dat"))
						rm(list=c("tmp.selblocks"))
						
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("-14 71 1",tmp.doitall,fixed=TRUE)
						tmp.doitall = tmp.doitall[-pointer]
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall","pointer"))

				# modify doitall, af95
					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
					tmp.doitall = gsub("2 95 5         # initial age structure based on Z for 1st 5 calendar years","2 95 2         # initial age structure based on Z for 1st 2 calendar years",tmp.doitall,fixed=TRUE)
					writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
					rm(list=c("tmp.doitall"))
				
				# modify ini, change L1 and L2
					if(testing.options$shift_growth[i] == 1)
					{
						tmp.ini=readLines(paste0(model.run.dir,"swo.ini"))
						pointer = grep("# ML1",tmp.ini,fixed=TRUE)+1
						tmp.ini[pointer] = paste0(L1_new," 50 150")
						pointer = grep("# ML2",tmp.ini,fixed=TRUE)+1
						tmp.ini[pointer] = paste0(L2_new," 130 300")
						writeLines(tmp.ini,con=paste0(model.run.dir,"swo.ini"))
						rm(list=c("tmp.ini","pointer"))
					}

				# modify doitall, increase weighting on composition data
					if(testing.options$increase_weighting_dwfn_comp[i] == 1)
					{
						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("-9 49 10      # except for fishery 9",tmp.doitall,fixed=TRUE)
						if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==0)
						{
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",jp_idx_fsh," 49 1"),paste0("    -",tw_idx_fsh," 49 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						}
						if(testing.options$TW_cpue[i]=="none"&testing.options$no_JP[i]==0){
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",jp_idx_fsh," 49 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						}
						if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==1){
							tmp.doitall = c(tmp.doitall[1:pointer],paste0("    -",tw_idx_fsh," 49 1"),tmp.doitall[(pointer+1):length(tmp.doitall)])
						} 
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall","pointer"))
					}

				# # modify doitall, length based selectivity
				# 	if(testing.options$length_selex[i]==1)
				# 	{
				# 		tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
				# 		pointer = grep("-999 26 2       # sets age-dependent selectivity option",tmp.doitall,fixed=TRUE)
				# 		tmp.doitall[pointer] = "   -999 26 3       # length based"
				# 		writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
				# 		rm(list=c("tmp.doitall","pointer"))
				# 	}

				# modify ini & doitall, estimate SDB
					if(testing.options$est_SDB[i] == 1)
					{
						tmp.ini=readLines(paste0(model.run.dir,"swo.ini"))
						pointer = grep("# Generic SD of length at age",tmp.ini,fixed=TRUE)+1
						tmp.ini[pointer] = "25 5 60"
						pointer = grep("# Length-dependent SD",tmp.ini,fixed=TRUE)+1
						tmp.ini[pointer] = "0 -1.5 1.5"
						writeLines(tmp.ini,con=paste0(model.run.dir,"swo.ini"))
						rm(list=c("tmp.ini","pointer"))

						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("1 32 6",tmp.doitall,fixed=TRUE)
						tmp.doitall = c(tmp.doitall[1:pointer],"      1 16 1          # estimate scalar of length dependent SD",tmp.doitall[(pointer+1):length(tmp.doitall)])
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall","pointer"))
					}

				# # modify doitall, estimate sel_dev_coffs for index fisheries
				# 	if(testing.options$sel_dev_coffs[i] == 1)
				# 	{
				# 		tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
				# 		pointer = grep("# Selectivity Section",tmp.doitall,fixed=TRUE) + 1
				# 		if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==0)
				# 		{
				# 			tmp.doitall = c(tmp.doitall[1:pointer],"   1 323 1",paste0("    -",jp_idx_fsh," 19 4"),paste0("    -",tw_idx_fsh," 19 4"),tmp.doitall[(pointer+1):length(tmp.doitall)])
				# 		}
				# 		if(testing.options$TW_cpue[i]=="none"&testing.options$no_JP[i]==0){
				# 			tmp.doitall = c(tmp.doitall[1:pointer],"   1 323 1",paste0("    -",jp_idx_fsh," 19 4"),tmp.doitall[(pointer+1):length(tmp.doitall)])
				# 		}
				# 		if(testing.options$TW_cpue[i]=="full"&testing.options$no_JP[i]==1){
				# 			tmp.doitall = c(tmp.doitall[1:pointer],"   1 323 1",paste0("    -",tw_idx_fsh," 19 4"),tmp.doitall[(pointer+1):length(tmp.doitall)])
				# 		} 
				# 		writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
				# 		rm(list=c("tmp.doitall","pointer"))
				# 	}

				# modify doitall, ff3
				tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
				tmp.doitall = gsub("-999 3 18",paste0("-999 3 ",testing.options$ff3[i]),tmp.doitall,fixed=TRUE)
				writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
				rm(list=c("tmp.doitall"))

				# modify doitall, remove grouping of fisheries 7 & 8 in new extraction fisheries structure
				tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
				# remove restriction for zero selex after age given by ff3
				pointer = grep("-7 16 2",tmp.doitall,fixed=TRUE)
				tmp.doitall = tmp.doitall[-pointer]
				# remove restriction for selex to be zero at early ages
				pointer = grep("-7 75 2",tmp.doitall,fixed=TRUE)
				tmp.doitall = tmp.doitall[-pointer]
				# regroup selectivities
				for(f in c(8:13,idx_fsh))
				{
					pointer = grep(paste0("-",f," 24 "),tmp.doitall,fixed=TRUE)
					tmp.doitall[pointer] = paste0("    -",f," 24 ",as.numeric(strsplit(trimws(tmp.doitall[pointer])," ")[[1]][3])+1)
				}
				writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
				rm(list=c("tmp.doitall"))

				# modify frq, fishery 9 length comp
					tmp.frqit = readfrq(paste0(model.run.dir,"swo.frq"))
					tmp.frqit2 = readfrq(paste0(dir.frq,"input_frq/2021_frq_5kg.swo.frq"))
					valid_fshinst = apply(cateffpen(tmp.frqit)[,c("year","month","fishery")],1,paste0,collapse="_")
					# update lf
						ln.dt2 = as.data.table(lnfrq(tmp.frqit2)) %>% .[fishery == 9]
						ln.dt = as.data.table(lnfrq(tmp.frqit)) %>% .[fishery!=9] %>% rbind(.,ln.dt2) %>% .[order(fishery,year,month,week)] %>% .[,id:=paste0(year,"_",month,"_",fishery)] %>% data.table::setcolorder(., neworder=c("id","year","month","week","fishery",as.character(seq(from=lf_range(tmp.frqit)['LFFirst'],by=lf_range(tmp.frqit)['LFWidth'],length.out=lf_range(tmp.frqit)['LFIntervals'])))) %>%
								.[id %in% valid_fshinst]
						ln.dt = ln.dt[, !"id"]
						lnfrq(tmp.frqit) = as.data.frame(ln.dt)

					writefrq(tmp.frqit,paste0(model.run.dir,"swo.frq"))
	  				tmp.frq = readLines(paste0(model.run.dir,"swo.frq"))
	  				tmp.frq[grep("#  Region     Fisheries    diffusion    tag groups",tmp.frq)+1] = paste0(c(n_regions(tmp.frqit), n_fisheries(tmp.frqit), as.numeric(generic_diffusion(tmp.frqit)), n_tag_groups(tmp.frqit),range(tmp.frqit)['minyear'], "1", as.numeric(frq_age_len(tmp.frqit)), n_recs_yr(tmp.frqit), rec_month(tmp.frqit), frq_version(tmp.frqit)), collapse='    ')
	  				writeLines(tmp.frq,paste0(model.run.dir,"swo.frq"))
	  				rm(list=c("ln.dt","valid_fshinst"))

	  			# modify doitall, group_f5_f11
	  				if(testing.options$group_f5_f11[i]==1)
	  				{
	  					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						tmp.doitall = gsub("-11 24 11","-11 24 5",tmp.doitall,fixed=TRUE)
						# regroup selectivities
						for(f in c(12:13,idx_fsh))
						{
							pointer = grep(paste0("-",f," 24 "),tmp.doitall,fixed=TRUE)
							tmp.doitall[pointer] = paste0("    -",f," 24 ",as.numeric(strsplit(trimws(tmp.doitall[pointer])," ")[[1]][3])-1)
						}
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall","pointer"))
	  				}

	  			# modify doitall & selblocks.dat, block_f11
	  				if(testing.options$block_f11[i]!=0)
	  				{
	  					if(testing.options$group_f5_f11[i]==1)
	  					{
	  						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
							pointer = grep("-4 71 1",tmp.doitall,fixed=TRUE)
							tmp.doitall = c(tmp.doitall[1:pointer],"     -5 71 1","    -11 71 1",tmp.doitall[(pointer+1):length(tmp.doitall)])
							writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
							rm(list=c("tmp.doitall","pointer"))

							tmp.selblocks = readLines(paste0(model.run.dir,"selblocks.dat"),warn=FALSE)
							writeLines(c(tmp.selblocks,rep(testing.options$block_f11[i],2)),con=paste0(model.run.dir,"selblocks.dat"))
							rm(list=c("tmp.selblocks"))
	  					} else {
	  						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
							pointer = grep("-4 71 1",tmp.doitall,fixed=TRUE)
							tmp.doitall = c(tmp.doitall[1:pointer],"    -11 71 1",tmp.doitall[(pointer+1):length(tmp.doitall)])
							writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
							rm(list=c("tmp.doitall","pointer"))

							tmp.selblocks = readLines(paste0(model.run.dir,"selblocks.dat"),warn=FALSE)
							writeLines(c(tmp.selblocks,rep(testing.options$block_f11[i],1)),con=paste0(model.run.dir,"selblocks.dat"))
							rm(list=c("tmp.selblocks"))
	  					}
	  				}

	  			# modify doitall & selblocks.dat, block_f12
	  				if(testing.options$block_f12[i]!=0)
	  				{

	  						tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
							pointer = grep("# Number of age classes starting from zero where selectivity is 0",tmp.doitall,fixed=TRUE) -1
							tmp.doitall = c(tmp.doitall[1:pointer],"    -12 71 1",tmp.doitall[(pointer+1):length(tmp.doitall)])
							writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
							rm(list=c("tmp.doitall","pointer"))

							tmp.selblocks = readLines(paste0(model.run.dir,"selblocks.dat"),warn=FALSE)
							writeLines(c(tmp.selblocks,rep(testing.options$block_f12[i],1)),con=paste0(model.run.dir,"selblocks.dat"))
							rm(list=c("tmp.selblocks"))
	  				}

				# modify doitall, rm_nodec_f3
	  				if(testing.options$rm_nodec_f3[i]==1)
	  				{
	  					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("-3 16 1",tmp.doitall,fixed=TRUE)
						tmp.doitall = tmp.doitall[-pointer]
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
	  				}

				# modify doitall, rm_0slex_f2
	  				if(testing.options$rm_0slex_f2[i]==1)
	  				{
	  					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("-2 75 1",tmp.doitall,fixed=TRUE)
						tmp.doitall = tmp.doitall[-pointer]
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
	  				}
				# modify doitall, rm_0slex_f8
	  				if(testing.options$rm_0slex_f8[i]==1)
	  				{
	  					tmp.doitall = readLines(paste0(model.run.dir,"condor_doitall.swo"),warn=FALSE)
						pointer = grep("-8 75 2",tmp.doitall,fixed=TRUE)
						tmp.doitall = tmp.doitall[-pointer]
						writeLines(tmp.doitall,con=paste0(model.run.dir,"condor_doitall.swo"))
						rm(list=c("tmp.doitall"))
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
						   .[,pf50:=as.numeric(NA)] %>%
						   .[,no_JP:=as.numeric(NA)] %>%
						   .[,TW_cpue	:=as.character(NA)] %>%
						   .[,JP_share_regsel:=as.numeric(NA)] %>%
						   .[,TW_share_regsel:=as.numeric(NA)] %>%
						   .[,shift_growth:=as.numeric(NA)] %>%
						   .[,start1994:=as.numeric(NA)] %>%
						   .[,increase_weighting_dwfn_comp:=as.numeric(NA)] %>%
						   .[,est_SDB:=as.numeric(NA)] %>%
						   .[,ff3:=as.numeric(NA)] %>%
						   .[,group_f5_f11:=as.numeric(NA)] %>%
						   .[,block_f11:=as.numeric(NA)] %>%
						   .[,block_f12:=as.numeric(NA)] %>%
						   .[,rm_nodec_f3:=as.numeric(NA)] %>%
						   .[,rm_0slex_f2:=as.numeric(NA)] %>%
						   .[,rm_0slex_f8:=as.numeric(NA)] %>%
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
		test.dt$version[i] = strsplit(dir.list[i],"_")[[1]][1]
		test.dt$pf50[i] = -as.numeric(strsplit(dir.list[i],"_")[[1]][2])
		test.dt$no_JP[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][3])
		test.dt$TW_cpue[i] = as.character(strsplit(dir.list[i],"_")[[1]][4])
		test.dt$JP_share_regsel[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][5])
		test.dt$TW_share_regsel[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][6])
		test.dt$shift_growth[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][7])
		test.dt$start1994[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][8])
		test.dt$increase_weighting_dwfn_comp[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][9])
		test.dt$est_SDB[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][10])
		test.dt$ff3[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][11])
		test.dt$group_f5_f11[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][12])
		test.dt$block_f11[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][13])
		test.dt$block_f12[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][14])
		test.dt$rm_nodec_f3[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][15])
		test.dt$rm_0slex_f2[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][16])
		test.dt$rm_0slex_f8[i] = as.numeric(strsplit(dir.list[i],"_")[[1]][17])
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
	
	# min aic is -103603.5
	# previous best model v2080_6_1_none_0_0_0_0_0_1_18_0_0_0_0_0_0 had aic of (-103595.1)
	test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&aic <= (-103603.5)+2][order(aic)]

	# stepwise selection
	# ff3
		ff3test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&group_f5_f11==0&block_f11==0&block_f12==0&rm_nodec_f3==0&rm_0slex_f2==0&rm_0slex_f8==0][order(aic)]$model_name
		# no model has an aic 2 points greater than the baseline model, so retain all models to next round
	# group_f5_f11
		test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&block_f11==0&block_f12==0&rm_nodec_f3==0&rm_0slex_f2==0&rm_0slex_f8==0][order(aic)]
		# 3 models have an aic 2 points greater than the baseline model, so retain those 3 models to next round
	# block_f11
		test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f12==0&rm_nodec_f3==0&rm_0slex_f2==0&rm_0slex_f8==0][order(aic)]
		# no models have an aic 2 points greater than the baseline model, so retain 3 original models to next round
	# block_f12
		test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0&rm_nodec_f3==0&rm_0slex_f2==0&rm_0slex_f8==0][order(aic)]
		# the spread of aic for all models is less than 2 so retain all models for next round
	# rm_nodec_f3
		test.dt[hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0&block_f12 %in% c(0,2010)&rm_0slex_f2==0&rm_0slex_f8==0][order(aic)]
		# eliminate 1 model for not being within 2 aic points of the best model v2080_6_1_none_0_0_0_0_0_1_14_1_0_0_1_0_0
	# rm_0slex_f2
		test.dt[model_name!="v2080_6_1_none_0_0_0_0_0_1_14_1_0_0_1_0_0"&hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0&block_f12 %in% c(0,2010)&rm_nodec_f3%in%c(0,1)&rm_0slex_f8==0][order(aic)]
		# eliminate all models 2 aic points more than lowest aic model
	# rm_0slex_f8
		test.dt[aic <= -103601.8 + 2 &model_name!="v2080_6_1_none_0_0_0_0_0_1_14_1_0_0_1_0_0"&hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0&block_f12 %in% c(0,2010)&rm_nodec_f3%in%c(0,1)&rm_0slex_f2%in%c(0,1)][order(aic)]
		# eliminate all models 2 aic points more than lowest aic model
	# best models
		test.dt[aic <= -103603.5 + 2 &model_name!="v2080_6_1_none_0_0_0_0_0_1_14_1_0_0_1_0_0"&hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0][order(aic)]
		best_models = test.dt[aic <= -103603.5 + 2 &model_name!="v2080_6_1_none_0_0_0_0_0_1_14_1_0_0_1_0_0"&hessian==0&bound_SDA==0&bound_SDB==0&bound_edev==0&ff3 %in% c(14,16,18) &group_f5_f11==1&block_f11==0][order(aic)]$model_name
		# best model: v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1; all things being equal according to aic, most parsimonious model with the best likelihood and run time of less than 3 hours
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

	tmp.dep = as.data.table(dep.ts) %>% .[,model:=rownames(dep.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","depletion")) %>% na.omit(.)
	add.dep = tmp.dep[model=="v2080_6_1_none_0_0_0_0_0_1_18_0_0_0_0_0_0"] %>% .[,year:=as.numeric(as.character(year))]


	p = tmp.dep  %>% .[,.(model,year,depletion)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,runtime_min)],by.x="model",by.y="model_name") %>% 
				.[,runtime_min:=ifelse(runtime_min<180,"0","1")] %>% .[model %in% best_models]
	p = p %>%			
				ggplot() +
	     		xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=depletion, color=model,linetype=runtime_min),size=2) + 
	     		geom_point(data=add.dep, aes(x=year, y=depletion),color="gray50",size=1.5) +
	  			expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(unique(p$model)))
	ggsave(filename="model.dep.lines.subset.best.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	tmp.bio = as.data.table(bio.ts) %>% .[,model:=rownames(bio.ts)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("year","bio")) %>% .[,bio:=bio/1000] %>% na.omit(.)
	add.bio = tmp.bio[model=="v2080_6_1_none_0_0_0_0_0_1_18_0_0_0_0_0_0"] %>% .[,year:=as.numeric(as.character(year))]


	p = tmp.bio  %>% .[,.(model,year,bio)] %>% .[,year:=as.numeric(as.character(year))] %>%  merge(.,test.dt[,.(model_name,runtime_min)],by.x="model",by.y="model_name") %>% 
				.[,runtime_min:=ifelse(runtime_min<180,"0","1")] %>% .[model %in% best_models]
	p = p %>%			
				ggplot() +
	     		xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=bio, color=model,linetype=runtime_min),size=2) + 
	     		geom_point(data=add.bio, aes(x=year, y=bio),color="gray50",size=1.5) +
	  			expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(unique(p$model)))
	ggsave(filename="model.bio.lines.subset.best.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	

	tmp.time = as.data.table(runtime_minutes.mat) %>% .[,model:=rownames(runtime_minutes.mat)] %>% melt(.,id.vars=c("model")) %>% setnames(.,c("variable","value"),c("phase","time")) 

	p = tmp.time  %>% .[,.(model,phase,time)] %>% .[,phase:=as.numeric(as.character(phase))] %>%  merge(.,test.dt[,.(model_name,runtime_min)],by.x="model",by.y="model_name") %>% 
				.[,runtime_min:=ifelse(runtime_min<180,"0","1")] %>% .[model %in% best_models] %>%
				ggplot() +
	     		xlab("Estimation phase") + ylab("Elapsed time (minutes)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=phase, y=time, color=model,linetype=runtime_min),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model", values=turbo_pal(best_models))
	ggsave(filename="model.etime.lines.subset.best.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
#___________________________________________________________________________________________________________
# plot selectivity...
	rep.list= lapply(paste0(paste0(dir.model_runs,jobs.group,c(best_models),"/"),"plot-",16,".par.rep"),read.MFCLRep)
	names(rep.list) = c(best_models)

	agg.years = TRUE
	agg.regions=TRUE
	biomass.type = "SSB"
	biomass.units=1000
	recdist.year_range=NULL
	yaxis.free=FALSE
	LRP=0.20
	TRP=NULL

	# fsh.lab = c(1:3,"4a","4b",5:14,"15a","15b",16:18)
	
	# selectivity - age
		p = diags4MFCL::plot.selectivity(rep.list,names(rep.list),sel.basis="AGE",palette.func=turbo_pal)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("sel.age.best.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# selectivity - length
		p = diags4MFCL::plot.selectivity(rep.list,names(rep.list),sel.basis="LENGTH",palette.func=turbo_pal)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("sel.length.best.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

#___________________________________________________________________________________________________________
# plot fit to length & weight
	fishery_names = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1","05_EU_1","06_Other_1","07_DW_2N","08_DW_2C","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU")

	source(paste0(proj.dir,"Utilities/parse.fit_file.r"))

	p = plot_fit_file(parse_fit_file(paste0(dir.model_runs,jobs.group,"v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","/"),type="length","swo"),type="length",fishery_names=fishery_names,save_dir=paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),save_name_stem="v2080_6_1_none_0_0_0_0_0_1_12")
	p = plot_fit_file(parse_fit_file(paste0(dir.model_runs,jobs.group,"v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","/"),type="weight","swo"),type="weight",fishery_names=fishery_names,save_dir=paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),save_name_stem="v2080_6_1_none_0_0_0_0_0_1_12")


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

	cpue_rmse.model_dir_list = paste0(dir.model_runs,jobs.group,c("v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1"),"/")
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
	
	rec_dt.list = as.list(rep(NA,length(cpue_rmse.model_dir_list)))
	relrec_dt.list = as.list(rep(NA,length(cpue_rmse.model_dir_list)))
	regpar_dt.list = as.list(rep(NA,length(cpue_rmse.model_dir_list)))

	# iterate over conventional models
	for(i in 1:length(cpue_rmse.model_dir_list))
	{	
		# read
		tmp.rep = read.MFCLRep(paste0(cpue_rmse.model_dir_list[i],"plot-16.par.rep"))
		rec_dt.list[[i]] = as.data.table(popN(tmp.rep)) %>% .[age=="1",.(year,area,value)] %>% setnames(.,"value","rec") %>% .[,year:=as.numeric(year)] %>% 
						  .[,model:=tail(strsplit(cpue_rmse.model_dir_list[i],"/")[[1]],n=1)] %>% .[,prop:=rec/sum(rec),by=year] %>% .[,.(model,year,area,rec,prop)]
		tmp.par = readLines(paste0(cpue_rmse.model_dir_list[i],"16.par"))
		pointer = grep("# relative recruitment",tmp.par) + 2
		relrec_dt.list[[i]] = data.table(model=tail(strsplit(cpue_rmse.model_dir_list[i],"/")[[1]],n=1),year=(1952:2015)[-1],relrec = as.numeric(strsplit(trimws(tmp.par[pointer]),"\\s+")[[1]])) # account for lag between spawners and recruits
		tmp.mat = matrix(scan(paste0(cpue_rmse.model_dir_list[i],"16.par"),skip=grep("# regional recruitment variation",tmp.par),nlines=64),ncol=2,byrow=TRUE)
		colnames(tmp.mat) = 1:2
		tmp.mat = as.data.table(tmp.mat) %>% .[,year:=1952:2015] %>% melt(.,id.vars = "year") %>% setnames(.,c("variable","value"),c("area","reg_rec_dev")) %>% .[,model:=tail(strsplit(cpue_rmse.model_dir_list[i],"/")[[1]],n=1)]
		rec_dt.list[[i]] = merge(rec_dt.list[[i]],tmp.mat)
		pointer = grep("# region parameters",tmp.par) + 1
		regpar_dt.list[[i]] = data.table(model=tail(strsplit(cpue_rmse.model_dir_list[i],"/")[[1]],n=1),area=as.character(1:2),avg_prop = as.numeric(strsplit(trimws(tmp.par[pointer]),"\\s+")[[1]]))

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
# launch likelihood profile
		if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	dir.list = paste0("v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","/")
	model.dir.list = paste0(dir.model_runs,jobs.group,dir.list)

	# get largest .par
	final.par = paste0(max(as.numeric(sapply(list.files(model.dir.list)[grep(".par$",list.files(model.dir.list))],function(x)strsplit(x,"[.]")[[1]][1]))),".par")
	final.rep = read.MFCLRep(paste0(model.dir.list,"plot-",final.par,".rep"))

	# define profile intervals
	down.int = seq(from=99,by=-1,length.out=60) # c(99, 95, 91, 87, 83, 79, 75, 71, 67, 63, 59, 55, 51)
	up.int = seq(from=101,by=1,length.out=60) # c(101, 105, 109, 113, 117, 121, 125, 129, 133, 137, 141, 145, 149)
	profile.mgc.flag = -4
	llprof.options = expand.grid(direction=c("up","down"),profile_type=c("biomass","depletion"),category=c("adult","total"),time_period=c("initial","average","recent","last_year"),stringsAsFactors=FALSE)
	llprof.options = subset(llprof.options, !(profile_type=="depletion" & time_period %in% c("average","initial")))
	rownames(llprof.options) = 1:nrow(llprof.options)

	# iterate across llprof.options
			ssh_exec_wait(session, command = paste0("mkdir ",dir.launch,jobs.group,"LLprof/"))
			
			for(i in 1:nrow(llprof.options))
			{
				runname = paste0(llprof.options[i,c(2,3,4,1)],collapse="/")
				model.run.dir = paste0(dir.model_runs,jobs.group,"LLprof/",runname,"/")

				# create new directory for model run
					if (! dir.exists(model.run.dir))dir.create(model.run.dir,recursive=TRUE)
					ssh_exec_wait(session, command = paste0("mkdir -p ",dir.launch,jobs.group,"LLprof/",runname))
				# transfer files
					FileList=c("mfcl.cfg","selblocks.dat","swo.frq",final.par)
					file.copy( paste0(model.dir.list,FileList),model.run.dir,overwrite=TRUE)

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
					if(supress_mfcl_IO)
					{
						tmp.bat = gsub("./doitall.condor","./doitall.condor &>/dev/null",tmp.bat,fixed=TRUE)
					}
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
		        	rm(list=c("runname","model.run.dir","FileList","TarList"))
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

	dir.list = unname(apply(llprof.options[,c(2,3,4,1)],1,function(x)paste0(x,collapse="/")))

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
# plot likelihood profile
	source(paste0(proj.dir,"Utilities/plot.likelihood.profiles.r"))

	plot.dir = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group)
	bio.units = 1000
	fishery_names = c("01_DW_1N","02_DW_1C","03_DW_1S","04_AU_1","05_EU_1","06_Other_1","07_DW_2N","08_DW_2C","09_DW_2S","10_NZ_2","11_EU_2","12_Other_2N","13_Other_2C","14_idx_AU","15_idx_NZ","16_idx_EU")

	for(i in 1:nrow(llprof.options))
	{
		llprof.dir.stem = paste0(proj.dir,"SWO/Assessment/Model_Runs/",jobs.group,"LLprof/",paste0(llprof.options[i,c(2,3,4)],collapse="/"),"/")
		if(llprof.options$profile_type[i] == "biomass")
		{
			ymax.total=40
			ymax.fsh=15
		} else {
			ymax.total=15
			ymax.fsh=5
		}
		plot.llprof(llprof.dir.stem,plot.dir,bio.units,fishery_names,ymax.total,ymax.fsh,profile_type=llprof.options$profile_type[i],category=llprof.options$category[i],time_period=llprof.options$time_period[i])
	}

#___________________________________________________________________________________________________________
# bivariate likelihood profile on M and steepness
	if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	dir.list = paste0("v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","/")
	model.dir.list = paste0(dir.model_runs,jobs.group,dir.list)

	# get largest .par
	final.par = paste0(max(as.numeric(sapply(list.files(model.dir.list)[grep(".par$",list.files(model.dir.list))],function(x)strsplit(x,"[.]")[[1]][1]))),".par")

	# define profile intervals
	m.int = c(0.25,0.5,0.75,1,1.25,1.5,1.75) 
	h.int = c(0.6,0.7,0.8,0.9,0.99)
	llprof.options = expand.grid(m_mult=m.int,h=h.int)

	# iterate across llprof.options
			ssh_exec_wait(session, command = paste0("mkdir ",dir.launch,jobs.group,"M_h_prof/"))
			
			for(i in 1:nrow(llprof.options))
			{
				runname = paste0(llprof.options[i,],collapse="_")
				model.run.dir = paste0(dir.model_runs,jobs.group,"M_h_prof/",runname,"/")

				# create new directory for model run
					if (! dir.exists(model.run.dir))dir.create(model.run.dir,recursive=TRUE)
					ssh_exec_wait(session, command = paste0("mkdir -p ",dir.launch,jobs.group,"M_h_prof/",runname))
				# transfer files
					FileList=c("mfcl.cfg","selblocks.dat","swo.frq",final.par)
					file.copy( paste0(model.dir.list,FileList),model.run.dir,overwrite=TRUE)
					file.rename(paste0(model.run.dir,final.par),paste0(model.run.dir,c("start.par")))

				# copy condor files
					file.copy(paste0(dir.condor,"max_iter/",c("condor.sub","mfcl.max_iter.bat","max_iter.sh")),model.run.dir,overwrite=TRUE)


				# modify start.par
					tmp.par = readLines(paste0(model.run.dir,"start.par"),warn=FALSE)
					pointer = grep("# natural mortality coefficient",tmp.par) + 2
					tmp.par[pointer] = format(as.numeric(tmp.par[pointer]) * llprof.options$m_mult[i],scientific=TRUE,digits=15)
					pointer = grep("# Seasonal growth parameters",tmp.par) + 1
					tmp.vec = sapply(tmp.par[pointer],function(x)strsplit(trimws(x),"\\s+")[[1]])
					tmp.vec[29] = format(llprof.options$h[i],scientific=TRUE,digits=15,nsmall=15)
					tmp.par[pointer] = paste0(" ",tmp.vec,collapse=" ")
					writeLines(tmp.par,con=paste0(model.run.dir,"start.par"))
					rm(list=c("tmp.par","pointer"))

				# modify .bat file
					tmp.bat = readLines(paste0(model.run.dir,"mfcl.max_iter.bat"),warn=FALSE)
					if(supress_mfcl_IO)
					{
						tmp.bat = gsub("./max_iter.sh","./max_iter.sh &>/dev/null",tmp.bat,fixed=TRUE)
					}
					writeLines(tmp.bat,con=paste0(model.run.dir,"mfcl.max_iter.bat"))
					rm(list=c("tmp.bat"))

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


				# make tar with the necessary files
					TarList=c("max_iter.sh","mfcl.cfg","selblocks.dat","swo.frq","start.par")
					shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -czf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
				# send all files to the launch machine
			       	scp_upload(session,files=paste0(model.run.dir,c("Start.tar.gz","condor.sub","mfcl.max_iter.bat")),to=paste0(dir.launch,jobs.group,"M_h_prof/",runname))
				# copy mfcl to launch machine run directory
			       	ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'v2080_dev_20210503/mfclo64 ',paste0(dir.launch,jobs.group,"M_h_prof/",runname)))
			    # remake tar file to include mfcl
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"M_h_prof/",runname),'; tar -xzf Start.tar.gz; rm Start.tar.gz'))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"M_h_prof/",runname),'; tar -czf Start.tar.gz ',paste(c(TarList,"mfclo64"),collapse=' ')))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"M_h_prof/",runname),'; rm ',paste(c(TarList,"mfclo64"),collapse=' ')))
				# launch condor job
			        ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"M_h_prof/",runname),'; dos2unix mfcl.max_iter.bat; chmod 700 mfcl.max_iter.bat; condor_submit condor.sub'))   
	
				# clean-up
		        	rm(list=c("runname","model.run.dir","FileList","TarList"))
			}

			ssh_disconnect(session)
	
#___________________________________________________________________________________________________________
# pull down finished likelihood profiles, untar, and clean directories
	if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}


	for(i in 1:nrow(llprof.options))
	{	

		runname = paste0(llprof.options[i,],collapse="_")
		model.run.dir = paste0(dir.model_runs,jobs.group,"M_h_prof/",runname,"/")
		condor.dir = paste0(dir.launch,jobs.group,"M_h_prof/",runname,"/")

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

#___________________________________________________________________________________________________________
# plot bivariate profile	
	
	bivariate_dt.list = as.list(rep(NA,nrow(llprof.options)))

	for(i in 1:nrow(llprof.options))
	{
		tmp.out = readLines(paste0(dir.model_runs,jobs.group,"M_h_prof/",paste0(llprof.options[i,],collapse="_"),"/","test_plot_output"),warn=FALSE)
		tmp.par = readLines(paste0(dir.model_runs,jobs.group,"M_h_prof/",paste0(llprof.options[i,],collapse="_"),"/","max_iter.par"),warn=FALSE)
		tmp.neigen = scan(paste0(dir.model_runs,jobs.group,"M_h_prof/",paste0(llprof.options[i,],collapse="_"),"/","neigenvalues"),nlines = 1)[1]
		n.fish = length(strsplit(trimws(tmp.out[grep("# Effort_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
		tmp.df = as.data.frame(matrix(0,nrow=n.fish,ncol=12))
		colnames(tmp.df) = c("m","h","fishery","total","bh","catchability","catch","effort","length","tag","weight","other")

		# fishery invariant
		# BH
			tmp.df$bh = rep(as.numeric(trimws(tmp.out[grep("# BH_steep contribution",tmp.out,fixed=TRUE)+1]))/n.fish,n.fish)
		# tag
			if(length(grep("# tag release", tmp.out))>0)
			{
				tmp.df$tag = rep(sum(scan(paste0(tmp.dir,"test_plot_output_",tmp.mult[j]),comment.char="#",skip=grep("# tag release",tmp.out,fixed=TRUE)[1],nlines=grep("# Tag likelihood for pooled tag release by fishery groups",tmp.out,fixed=TRUE) - grep("# tag release",tmp.out,fixed=TRUE)[1]))/n.fish,n.fish)
			}
		# fishery specific
		# effort dev
			tmp.df$effort = as.numeric(strsplit(trimws(tmp.out[grep("# Effort_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
		# catchability
			tmp.df$catchability = as.numeric(strsplit(trimws(tmp.out[grep("# catchability_dev_penalty_by_fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
		# length
			tmp.df$length = as.numeric(strsplit(trimws(tmp.out[grep("# total length component of likelihood for each fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
		# weight
			tmp.df$weight = as.numeric(strsplit(trimws(tmp.out[grep("# total weight component of likelihood for each fishery",tmp.out,fixed=TRUE)+1])," ")[[1]])
		# catch
			tmp.df$catch = sapply(strsplit(tmp.out[grep("total catch components of likelihood for fishery", tmp.out)+1], split="[[:blank:]]+"), function(x)sum(as.numeric(x))) 

		# total
			tmp.df$total = rep(-as.numeric(tmp.par[grep("# Objective function value",tmp.par)+1])/n.fish,n.fish)
			tmp.df$other = rep((sum(tmp.df$total) - (sum(tmp.df$bh) + sum(tmp.df$tag) + sum(tmp.df$effort) + sum(tmp.df$catchability) + sum(tmp.df$length) + sum(tmp.df$weight) + sum(tmp.df$catch)))/n.fish,n.fish)	
			tmp.df$m= as.numeric(tmp.par[grep("# natural mortality coefficient",tmp.par)+2])
			tmp.df$h= llprof.options$h[i]
			tmp.df$fishery = 1:n.fish
			tmp.df$PDH = ifelse(tmp.neigen>0,"0","1")
			bivariate_dt.list[[i]] = tmp.df		
			rm(list=c("tmp.out","n.fish","tmp.df","tmp.par","tmp.neigen"))
	}

	bivariate_dt = rbindlist(bivariate_dt.list)

		p = bivariate_dt %>% .[,.(PDH=unique(PDH),total=sum(total)),by=.(m,h)] %>% .[,delta_NLL:=abs(total-min(total))] %>% 
		ggplot() + 
	    xlab("Natural mortality (M)") + ylab("Steepness (h)") + 
	    geom_tile(aes(x=m,y=h,fill=delta_NLL)) +
	    geom_point(aes(x=m,y=h,color=PDH),size=2) + 
	  	theme_few(base_size = 20) +
	  	scale_fill_gradientn("Change in likelihood",colours = turbo_vec(100,start = 0.1, end = 0.8)) +
	  	scale_color_manual("PDH",values=c("red","black"))
	  	ggsave(filename="bivariate.m_h.total_likelihood.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	
	  	p = bivariate_dt %>% .[,lapply(.SD,sum),by=.(m,h),.SDcols = c("total","bh","catchability","catch","effort","length","tag","weight","other")] %>%
	  					 merge(.,unique(bivariate_dt[,.(m,h)]),by=c("m","h")) %>%
	  					 .[,total:=abs(total-min(total)),by=m] %>% .[,bh:=abs(bh-min(bh)),by=m] %>% .[,catchability:=abs(catchability-min(catchability)),by=m] %>%
						 .[,catch:=abs(catch-min(catch)),by=m] %>% .[,effort:=abs(effort-min(effort)),by=m] %>% .[,length:=abs(length-min(length)),by=m] %>%
						 .[,tag:=abs(tag-min(tag)),by=m] %>% .[,weight:=abs(weight-min(weight)),by=m] %>% .[,other:=abs(other-min(other)),by=m] %>% melt(.,id.vars=c("m","h")) %>% .[order(m,h),] %>% .[variable %in% c("total","bh","effort","length","weight","other")] %>%
						 .[,variable:=factor(variable,levels=c("total","bh","effort","length","weight","other"),labels=c("Total","BH-SRR","Effort","Length","Weight","Other penalties"))] %>%
						 ggplot() + 
					    xlab("Steepness (h)") + ylab("Change in likelihood") +
					    geom_hline(yintercept=0) + facet_wrap(~variable,scales="free_y") + 
					    geom_line(aes(x=h,y=value,color=m,group=m),size=1) +
					  	theme_few(base_size = 20) +
					  	scale_color_gradientn("Natural mortality (M)",colours = turbo_vec(100,start = 0.1, end = 0.8))
				ggsave(filename="bivariate.h_profile.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	  	p = bivariate_dt %>% .[,lapply(.SD,sum),by=.(m,h),.SDcols = c("total","bh","catchability","catch","effort","length","tag","weight","other")] %>%
	  					 merge(.,unique(bivariate_dt[,.(m,h)]),by=c("m","h")) %>%
	  					 .[,total:=abs(total-min(total)),by=h] %>% .[,bh:=abs(bh-min(bh)),by=h] %>% .[,catchability:=abs(catchability-min(catchability)),by=h] %>%
						 .[,catch:=abs(catch-min(catch)),by=h] %>% .[,effort:=abs(effort-min(effort)),by=h] %>% .[,length:=abs(length-min(length)),by=h] %>%
						 .[,tag:=abs(tag-min(tag)),by=h] %>% .[,weight:=abs(weight-min(weight)),by=h] %>% .[,other:=abs(other-min(other)),by=h] %>% melt(.,id.vars=c("m","h")) %>% .[order(h,m),] %>% .[variable %in% c("total","bh","effort","length","weight","other")] %>%
						 .[,variable:=factor(variable,levels=c("total","bh","effort","length","weight","other"),labels=c("Total","BH-SRR","Effort","Length","Weight","Other penalties"))] %>%
						 ggplot() + 
					    xlab("Natural mortality (M)") + ylab("Change in likelihood") +
					    geom_hline(yintercept=0) + facet_wrap(~variable,scales="free_y") + 
					    geom_line(aes(x=m,y=value,color=h,group=h),size=1) +
					  	theme_few(base_size = 20) +
					  	scale_color_gradientn("Steepness (h)",colours = turbo_vec(100,start = 0.1, end = 0.8))
				ggsave(filename="bivariate.m_profile.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)

	rep_dt.list = as.list(rep(NA,nrow(llprof.options)))
	for(i in 1:length(rep_dt.list))
	{
		tmp.rep = read.MFCLRep(paste0(dir.model_runs,jobs.group,"M_h_prof/",paste0(llprof.options[i,],collapse="_"),"/","plot-max_iter.par.rep"))
		tmp.dt = data.table(model=paste0(llprof.options[i,],collapse="_"),sb=as.vector(SB(rep=tmp.rep,mean_nyears=1,lag_nyears=0)),
		sb_dep=as.vector(SBSBF0(rep=tmp.rep,sb_mean_nyears=1, sb_lag_nyears=0, sbf0_mean_nyears=1, sbf0_lag_nyears=0)),
		m=llprof.options$m[i],h=llprof.options$h[i])
		tmp.dt$year = seq(from=1952,by=1,length.out=nrow(tmp.dt))
		tmp.dt$PDH = ifelse(scan(paste0(dir.model_runs,jobs.group,"M_h_prof/",paste0(llprof.options[i,],collapse="_"),"/","neigenvalues"),nlines = 1)[1]>0,"0","1")
		rep_dt.list[[i]] = tmp.dt[,.(model,m,h,PDH,year,sb,sb_dep)]
		rm("tmp.dt","tmp.rep")
	}
	rep_dt = rbindlist(rep_dt.list)


		  	p = rep_dt %>%
		  		ggplot() +
				xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=sb_dep, group=model,color=as.numeric(m),linetype=PDH),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_linetype_manual("PDH",values = c("dashed","solid")) +
	  			scale_colour_gradientn("Natural mortality (M) multiplier", colors=rev(turbo_vec(100,start = 0.1, end = 0.8)))
				ggsave(filename="bivariate.m.dep.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)
	  	 	p = rep_dt %>%
		  		ggplot() +
				xlab("Year") + ylab(expression("SB"/"SB"["F=0"])) + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=sb_dep, group=model,color=as.numeric(h),linetype=PDH),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_linetype_manual("PDH",values = c("dashed","solid")) +
	  			scale_colour_gradientn("Steepness (h)", colors=rev(turbo_vec(100,start = 0.1, end = 0.8)))
				ggsave(filename="bivariate.h.dep.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)
	
		  	p = rep_dt %>%
		  		ggplot() +
				xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=sb/1000, group=model,color=as.numeric(m),linetype=PDH),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_linetype_manual("PDH",values = c("dashed","solid")) +
	  			scale_colour_gradientn("Natural mortality (M) multiplier", colors=rev(turbo_vec(100,start = 0.1, end = 0.8)))
				ggsave(filename="bivariate.m.bio.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)
	  	 	p = rep_dt %>%
		  		ggplot() +
				xlab("Year") + ylab("Spawning potential (1000's mt)") + geom_hline(yintercept = 0) +
	  			geom_line( aes(x=year, y=sb/1000, group=model,color=as.numeric(h),linetype=PDH),size=2) + expand_limits(y=0) +
	  			theme_few(base_size = 20) +
	  			scale_linetype_manual("PDH",values = c("dashed","solid")) +
	  			scale_colour_gradientn("Steepness (h)", colors=rev(turbo_vec(100,start = 0.1, end = 0.8)))
				ggsave(filename="bivariate.h.bio.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  	 		scale = 1.25, width = 16, height = 9, units = c("in"),
	  	 		dpi = 300, limitsize = TRUE)


#___________________________________________________________________________________________________________
# ASPM diagnostic
	if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}

	supress_mfcl_IO = FALSE	
	dir.list = paste0("v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","/")
	model.dir.list = paste0(dir.model_runs,jobs.group,dir.list)

	# get largest .par
	final.par = paste0(max(as.numeric(sapply(list.files(model.dir.list)[grep(".par$",list.files(model.dir.list))],function(x)strsplit(x,"[.]")[[1]][1]))),".par")
	aspm.options = expand.grid(type=c("basic","rdev"))

	# iterate across aspm.options
			ssh_exec_wait(session, command = paste0("mkdir ",dir.launch,jobs.group,"aspm/"))
			
			for(i in 1:nrow(aspm.options))
			{
				runname = aspm.options$type[i]
				model.run.dir = paste0(dir.model_runs,jobs.group,"aspm/",runname,"/")

				# create new directory for model run
					if (! dir.exists(model.run.dir))dir.create(model.run.dir,recursive=TRUE)
					ssh_exec_wait(session, command = paste0("mkdir -p ",dir.launch,jobs.group,"aspm/",runname))
				# transfer files
					FileList=c("mfcl.cfg","selblocks.dat","swo.frq",final.par)
					file.copy( paste0(model.dir.list,FileList),model.run.dir,overwrite=TRUE)
					file.rename(paste0(model.run.dir,final.par),paste0(model.run.dir,c("start.par")))

				# copy condor files
					file.copy(paste0(dir.condor,"max_iter/",c("condor.sub","mfcl.max_iter.bat","max_iter.sh")),model.run.dir,overwrite=TRUE)


				# modify doitall
					tmp.sh = readLines(paste0(model.run.dir,"max_iter.sh"),warn=FALSE)
					pointer = grep("$MFCL $FRQ $PAR tmp.par -switch 7 1 1 10000 1 50 -8 1 186 0 1 187 0 1 188 0 1 189 0 1 190 0",tmp.sh,fixed=TRUE)
					if(aspm.options$type[i] == "basic")
					{
						tmp.sh[pointer] = "$MFCL $FRQ $PAR tmp.par -switch 17 1 1 10000 1 50 -8 1 186 0 1 187 0 1 188 0 1 189 0 1 190 0 -999 48 0 -999 49 100000 -999 50 100000 1 15 0 1 16 0 2 30 0 2 70 0 2 71 0 -100000 1 0 -100000 2 0"
					} else {
						tmp.sh[pointer] = "$MFCL $FRQ $PAR tmp.par -switch 12 1 1 10000 1 50 -8 1 186 0 1 187 0 1 188 0 1 189 0 1 190 0 -999 48 0 -999 49 100000 -999 50 100000 1 15 0 1 16 0"
					}
					writeLines(tmp.sh,con=paste0(model.run.dir,"max_iter.sh"))
					rm(list=c("tmp.sh","pointer"))

				# modify start.par
					if(aspm.options$type[i] == "basic")
					{
						tmp.par = readLines(paste0(model.run.dir,"start.par"),warn=FALSE)
						pointer = grep("# relative recruitment",tmp.par,fixed=TRUE) + 2
						tmp.par[pointer] = paste0(rep(1,length(strsplit(trimws(tmp.par[pointer])," ")[[1]])),collapse=" ")
						tmp.par[grep("# regional recruitment variation",tmp.par) + (1:68)] = " 0 0"
						writeLines(tmp.par,con=paste0(model.run.dir,"start.par"))
						rm(list=c("tmp.par","pointer"))
					} 

				# modify .bat file
					tmp.bat = readLines(paste0(model.run.dir,"mfcl.max_iter.bat"),warn=FALSE)
					if(supress_mfcl_IO)
					{
						tmp.bat = gsub("./max_iter.sh","./max_iter.sh &>/dev/null",tmp.bat,fixed=TRUE)
					}
					writeLines(tmp.bat,con=paste0(model.run.dir,"mfcl.max_iter.bat"))
					rm(list=c("tmp.bat"))
					
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


				# make tar with the necessary files
					TarList=c("max_iter.sh","mfcl.cfg","selblocks.dat","swo.frq","start.par")
					shell(paste0("cd ",model.run.dir,"& C:/cygwin64/bin/tar.exe -czf Start.tar.gz ",paste(TarList,collapse=' ')),translate=TRUE)
				# send all files to the launch machine
			       	scp_upload(session,files=paste0(model.run.dir,c("Start.tar.gz","condor.sub","mfcl.max_iter.bat")),to=paste0(dir.launch,jobs.group,"aspm/",runname))
				# copy mfcl to launch machine run directory
			       	ssh_exec_wait(session,command=paste0('cp ',dir.mfcl.launch,'v2080_dev_20210503/mfclo64 ',paste0(dir.launch,jobs.group,"aspm/",runname)))
			    # remake tar file to include mfcl
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"aspm/",runname),'; tar -xzf Start.tar.gz; rm Start.tar.gz'))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"aspm/",runname),'; tar -czf Start.tar.gz ',paste(c(TarList,"mfclo64"),collapse=' ')))
			       	ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"aspm/",runname),'; rm ',paste(c(TarList,"mfclo64"),collapse=' ')))
				# launch condor job
			        ssh_exec_wait(session,command=paste0('cd ',paste0(dir.launch,jobs.group,"aspm/",runname),'; dos2unix mfcl.max_iter.bat; chmod 700 mfcl.max_iter.bat; condor_submit condor.sub'))   
	
				# clean-up
		        	rm(list=c("runname","model.run.dir","FileList","TarList"))
			}

			ssh_disconnect(session)
	
#___________________________________________________________________________________________________________
# pull down finished likelihood profiles, untar, and clean directories
	if(submit.realm == "SUV")
		{
			session = ssh_connect("nicholasd@SUVOFPSUBMIT")
		} else {
			session = ssh_connect("nicholasd@NOUOFPCALC02")
		}


	for(i in 1:nrow(aspm.options))
	{	

		runname = paste0(aspm.options[i,],collapse="_")
		model.run.dir = paste0(dir.model_runs,jobs.group,"aspm/",runname,"/")
		condor.dir = paste0(dir.launch,jobs.group,"aspm/",runname,"/")

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

#___________________________________________________________________________________________________________
# plot aspm results
	library(diags4MFCL)
	library(frqit)

	model.path = paste0(dir.model_runs,jobs.group,c("v2080_6_1_none_0_0_0_0_0_1_16_1_0_0_1_0_1","aspm/rdev"),"/")
	rep.list= lapply(paste0(model.path,"plot-",c(16,"max_iter"),".par.rep"),read.MFCLRep)
	names(rep.list) = c("Full model","ASPM-R")
	model.names = names(rep.list)

	agg.years = TRUE
	agg.regions=TRUE
	biomass.type = "SSB"
	biomass.units=1000
	recdist.year_range=NULL
	yaxis.free=FALSE
	LRP=0.20
	TRP=NULL
	source(paste0(proj.dir,"Utilities/turbo.r"))

# depletion - regional
		p = plot.depletion(rep.list,model.names,agg.years,agg.regions=FALSE,biomass.type=biomass.type,LRP=LRP,TRP=TRP)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("aspm.dep.reg.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
				
	# depletion
		p = plot.depletion(rep.list,model.names,agg.years,agg.regions=TRUE,biomass.type=biomass.type,LRP=LRP,TRP=TRP)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("aspm.dep.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
				
	# biomass - regional
		p = plot.biomass(rep.list,model.names,agg.years,agg.regions=FALSE,biomass.type=biomass.type,biomass.units=biomass.units)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("aspm.bio.reg.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)
				
	# biomass
		p = plot.biomass(rep.list,model.names,agg.years,agg.regions=TRUE,biomass.type=biomass.type,biomass.units=biomass.units)
		p = p + theme_few(base_size = 20) 
		ggsave(filename=paste0("aspm.bio.png"), plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),scale = 1, width = 16, height = 9, units = c("in"),dpi = 300, limitsize = TRUE)

	# fit2cpue
			tmp_dt = cbind(cateffpen(readfrq(paste0(model.path[1],"swo.frq"))),unlist(effort_dev_coffs(read.MFCLPar(paste0(model.path[1],"16.par")))))
			colnames(tmp_dt)[8] = "edev"
			tmp_dt = as.data.table(tmp_dt) %>% .[,model:="Full model"]
			tmp_dt2 = cbind(cateffpen(readfrq(paste0(model.path[2],"swo.frq"))),unlist(effort_dev_coffs(read.MFCLPar(paste0(model.path[2],"max_iter.par")))))
			colnames(tmp_dt2)[8] = "edev"
			tmp_dt2 = as.data.table(tmp_dt2) %>% .[,model:="ASPM-R"]
			tmp_dt = rbind(tmp_dt,tmp_dt2)


			tmp.dt = tmp_dt %>% .[catch==-1,catch:=NA] %>% .[effort==-1,effort:=NA] %>% .[!is.na(effort) & !is.na(catch)] %>%
					 .[,mean_effort := mean(effort),by=.(fishery,model)] %>% .[,norm_effort:=effort/mean(effort),by=.(fishery,model)] %>% .[,cpue_obs:=catch/norm_effort] %>% .[,cpue_pred := catch/(norm_effort*exp(edev))] %>%
					 .[,input_cv := 1/(sqrt(2*penalty))] %>% .[,l_se:=exp(log(cpue_obs)-input_cv)] %>% .[,u_se:=exp(log(cpue_obs)+input_cv)] %>% .[,ts:=year+(month-1)/12] %>%
					 .[,fishery:=factor(as.character(fishery),levels=as.character(sort(unique(fishery))),labels=fishery_names)] %>% 
					 .[,max_u_se:=max(u_se,cpue_pred),by=fishery] %>% .[fishery %in% c("14_idx_AU","15_idx_NZ","16_idx_EU")]
		
			# metrics
				tmp_met.dt = tmp.dt %>% .[,.(mse=mean((cpue_obs - cpue_pred)^2),rmse=sqrt(mean((cpue_obs - cpue_pred)^2)),nrmse=sqrt(mean((cpue_obs - cpue_pred)^2))/mean(cpue_obs),max_u_se=max(max_u_se)),by=.(fishery,model)] %>%
				 			            .[,x:=min(tmp.dt$ts)] %>% .[,x2:=max(tmp.dt$ts)]
				tmp_met.dt1 = tmp_met.dt %>% .[model == "Full model"]
				tmp_met.dt2 = tmp_met.dt %>% .[model == "ASPM-R"]

				obs.dt = tmp.dt %>% .[,.(fishery,ts,cpue_obs,l_se,u_se,penalty)] %>% unique(.)
				
			# plots 
				p = tmp.dt %>%  
				ggplot() +
	     		xlab("Year") + ylab("Standardized CPUE") + facet_wrap(~fishery,scales="free_y") +
	     		geom_hline(yintercept = 0) +
	     		geom_segment(data=obs.dt,aes(x=ts,xend=ts, y=l_se,yend=u_se),color="gray65",size=0.75) + 
	  			geom_point(data=obs.dt, aes(x=ts, y=cpue_obs),size=3,fill="gray90",color="black",shape=21) + expand_limits(y=c(0)) +
	  			geom_line( aes(x=ts, y=cpue_pred,color=model),size=1.25) +
	  			theme_few(base_size = 20) +
	  			scale_colour_manual("Model",values = c("#2979ff","black")) +
	  			geom_rect(data=tmp_met.dt1,aes(xmin = x, xmax = 0.4*(x2-x)+x, ymin = max_u_se*0.85, ymax = max_u_se*1.05),fill="white",alpha=0.85) +
				geom_text(data=tmp_met.dt2,aes(x=x,y=max_u_se,label=paste0("CV (RMSE): ",round(nrmse,digits=3))),color="#2979ff",hjust="inward") +
	  			geom_text(data=tmp_met.dt1,aes(x=x,y=max_u_se*0.95,label=paste0("CV (RMSE): ",round(nrmse,digits=3))),color="black",hjust="inward")
				ggsave(filename="aspm.fit2cpue.fishery.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# recruitment
		rec_dt.list = as.list(rep(NA,length(rep.list)))

		# iterate over conventional models
		for(i in 1:length(rep.list))
		{	
			# read
			tmp.rep = rep.list[[i]]
			rec_dt.list[[i]] = as.data.table(popN(tmp.rep)) %>% .[age=="1",.(year,area,value)] %>% setnames(.,"value","rec") %>% .[,year:=as.numeric(year)] %>% 
							  .[,model:=names(rep.list[i])] %>% .[,prop:=rec/sum(rec),by=year] %>% .[,.(model,year,area,rec,prop)]
			
			rm(list=c("tmp.rep"))
		}

		rec_dt = rbindlist(rec_dt.list)
		# total reg rec
			p = rec_dt %>% 
			ggplot() + facet_wrap(~area) +
		    xlab("Year") + ylab("Recruitment (millions)") + 
		    geom_hline(yintercept = 0) +
		  	geom_line( aes(x=year, y=rec/1000000, color=model),size=1.15) + expand_limits(y=c(0)) +
		  	theme_few(base_size = 20) +
	  		scale_colour_manual("Model",values = c("#2979ff","black"))
		  	ggsave(filename="aspm.total_rec.reg.png", plot = p, device = "png", path = paste0(proj.dir,"SWO/Assessment/Model_Plots/",jobs.group),
		  			scale = 1.25, width = 16, height = 9, units = c("in"),
		  			dpi = 300, limitsize = TRUE)

		  	