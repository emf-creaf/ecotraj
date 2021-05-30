# .representativeTrajectory2<-function(lsd, clcore, eps) {
#   D_core = as.matrix(lsd$Dseg)[clcore,clcore, drop=FALSE]
#   nameseg = row.names(D_core)
#   Dini_core = as.matrix(lsd$Dini)[clcore,clcore, drop=FALSE]
#   Dfin_core = as.matrix(lsd$Dfin)[clcore,clcore, drop=FALSE]
#   Dfinini_core = lsd$Dfinini[clcore,clcore, drop=FALSE]
#   Dinifin_core = lsd$Dinifin[clcore,clcore, drop=FALSE]
#   ncore = sum(clcore)
#   projmatIni<-matrix(NA, nrow=ncore, ncore)
#   projmatFin<-matrix(NA, nrow=ncore, ncore)
#   neighborhood<-rep(0,ncore)
#   for(i in 1:ncore) {
#     for(j in 1:ncore) {
#       if(D_core[i,j]<eps) neighborhood[i]=neighborhood[i]+1
#       pini = .distanceToSegment(Dfinini_core[i,i], Dini_core[i,j], Dfinini_core[i,j])
#       projmatIni[i,j] = pini[1]/Dfinini_core[i,i]
#       pfin = .distanceToSegment(Dfinini_core[i,i], Dinifin_core[i,j], Dfin_core[i,j])
#       projmatFin[i,j] = pfin[1]/Dfinini_core[i,i]
#     }
#   }
#   print(neighborhood)
#   
#   #Initial medoid segment
#   dists = numeric(0)
#   trajvec = character(0)
#   inivec = logical(0)
#   med = which.max(neighborhood)[1]
#   print(med)
#   ini = TRUE
#   
#   piniTomed = (projmatIni[med,] > 0) & (projmatIni[med,]<=1)
#   pfinTomed = (projmatFin[med,] > 0) & (projmatFin[med,]<=1)
#   
#   while(sum(piniTomed)+ sum(pfinTomed)>0) {
#     maxNeigh = max(neighborhood[piniTomed],neighborhood[pfinTomed])
#     indini = which(piniTomed)[which(neighborhood[piniTomed] == maxNeigh)]
#     indfin = which(pfinTomed)[which(neighborhood[pfinTomed] == maxNeigh)]
#     if(ini) {
#       dindini = Dini_core[med,indini]
#       dindfin = Dinifin_core[med,indfin]
#     } else {
#       dindini = Dfinini_core[med,indini]
#       dindfin = Dfin_core[med,indfin]
#     }
#     mind = which.min(c(dindini,dindfin))[1]
#     if(mind<=length(indini)) {
#       med = indini[mind]
#       ini = TRUE
#     } else {
#       med = indfin[mind-length(indini)]
#       ini = FALSE
#     }
#     dists = c(dists, mind)
#     trajvec = c(trajvec,nameseg[med])
#     inivec = c(inivec, TRUE)
#     print(rbind(trajvec,inivec,dists))
#     piniTomed = (projmatIni[med,] > 0) & (projmatIni[med,]<=1)
#     pfinTomed = (projmatFin[med,] > 0) & (projmatFin[med,]<=1)
#   }
#   
#   return(data.frame(segments = trajvec, ini = inivec, distances = dists))
# }
# .representativeTrajectory<-function(lsd, clcore, eps){
#   
#   D_core = as.matrix(lsd$Dseg)[clcore,clcore, drop=FALSE]
#   Dini_core = as.matrix(lsd$Dini)[clcore,clcore, drop=FALSE]
#   Dfin_core = as.matrix(lsd$Dfin)[clcore,clcore, drop=FALSE]
#   Dfinini_core = lsd$Dfinini[clcore,clcore, drop=FALSE]
#   Dinifin_core = lsd$Dinifin[clcore,clcore, drop=FALSE]
#   ncore = sum(clcore)
#   ## Determine medoid
#   med = which.min(rowSums(D_core))
#   # cat("1",med,"\n")
#   trajvec = med
#   dists = Dinifin_core[med,med]
#   if(ncore>1) {
#     which.mod <-function(d, eps) {
#       w = which(d<eps)
#       if(length(w)==0) w = which.min(d)
#       return(w)
#     }
#     ## Add segments after medoid  
#     # print(Dfinini_core)
#     # print(eps)
#     nb = apply(Dfinini_core,1,which.mod, eps)
#     if(length(nb)>0) {
#       nbm = nb[[med]]
#       nbm = nbm[!(nbm %in% trajvec)]
#       # print(nbm)
#       if(length(nbm)>0) {
#         selnb = rep(FALSE,length(nbm))
#         for(i in 1:length(nbm)) {
#           # cat(nbm[i],Dinifin_core[med,med], Dini_core[med,nbm[i]], Dfinini_core[med,nbm[i]],"\n")
#           dfin = .distanceToSegmentC(Dinifin_core[med,med], Dini_core[med,nbm[i]], Dfinini_core[med,nbm[i]])
#           selnb[i] = (dfin[1] == Dinifin_core[med,med])
#         }
#         nbm = nbm[selnb]
#       }
#     } else {
#       nbm = numeric(0)
#     }
#     
#     while(length(nbm)>0) {
#       # print(trajvec)
#       medp = med
#       # cat("2",med,"\n")
#       med = nbm[which.min(rowSums(D_core[nbm,nbm, drop=FALSE]))]
#       trajvec = c(trajvec,med)
#       dists = c(dists,Dfinini_core[medp,med] ,Dinifin_core[med,med])
#       nbm = nb[[med]]
#       nbm = nbm[!(nbm %in% trajvec)]
#       # print(nbm)
#       if(length(nbm)>0) {
#         selnb = rep(FALSE,length(nbm))
#         for(i in 1:length(nbm)) {
#           # cat(nbm[i],Dinifin_core[med,med], Dini_core[med,nbm[i]], Dfinini_core[med,nbm[i]],"\n")
#           dfin = .distanceToSegmentC(Dinifin_core[med,med], Dini_core[med,nbm[i]], Dfinini_core[med,nbm[i]])
#           selnb[i] = (dfin[1] == Dinifin_core[med,med])
#         }
#         nbm = nbm[selnb]
#       }
#     }
#     
#     ## Add segments before medoid
#     nb = apply(Dinifin_core,1,which.mod, eps)
#     if(length(nb)>0) {
#       med  =trajvec[1]
#       nbm = nb[[med]]
#       nbm = nbm[!(nbm %in% trajvec)]
#       # print(nbm)
#       if(length(nbm)>0) {
#         selnb = rep(FALSE,length(nbm))
#         for(i in 1:length(nbm)) {
#           # cat(nbm[i],Dinifin_core[med,med], Dinifin_core[med,nbm[i]], Dfin_core[med,nbm[i]],"\n")
#           dfin = .distanceToSegmentC(Dinifin_core[med,med], Dinifin_core[med,nbm[i]], Dfin_core[med,nbm[i]])
#           selnb[i] = (dfin[1] == 0)
#         }
#         nbm = nbm[selnb]
#       }    
#     } else {
#       nbm = numeric(0)
#     }
#     
#     while(length(nbm)>0) {
#       medp = med
#       med = nbm[which.min(rowSums(D_core[nbm,nbm, drop=FALSE]))]
#       trajvec = c(med,trajvec)
#       dists = c(Dinifin_core[med,med],Dfinini_core[medp,med],dists)
#       
#       # cat("3",med,"\n")
#       nbm = nb[[med]]
#       nbm = nbm[!(nbm %in% trajvec)]
#       # print(nbm)
#       if(length(nbm)>0) {
#         selnb = rep(FALSE,length(nbm))
#         for(i in 1:length(nbm)) {
#           # cat(nbm[i],Dinifin_core[med,med], Dinifin_core[med,nbm[i]], Dfin_core[med,nbm[i]],"\n")
#           dfin = .distanceToSegmentC(Dinifin_core[med,med], Dinifin_core[med,nbm[i]], Dfin_core[med,nbm[i]])
#           selnb[i] = (dfin[1] == 0)
#         }
#         nbm = nbm[selnb]
#       }
#     }
#   }
#   return(list(segments = names(trajvec), distances = dists))
# }
# 
# .representativeTrajectories<-function(lsd, cluster, core, eps) {
#   clids <-unique(cluster)
#   clids<-clids[!is.na(clids)]
#   res<-vector("list", length(clids))
#   for(i in 1:length(clids)) {
#     clid = clids[i]
#     clcore = (cluster==clid) & (core)
#     clcore[is.na(clcore)] = FALSE
#     if(sum(clcore)>0) {
#       res[[i]]<-.representativeTrajectory(lsd,clcore,eps)
#     }
#   }
#   return(res)
# }
# 
# .trajectoryClustering<-function(d, nt, eps, minlns, distance.type = "directed-segment") {
#   lsd = segmentDistances(d,nt,distance.type = distance.type)
#   dmat = as.matrix(lsd$Dseg)
#   nb = apply(dmat<=eps,1,which)
#   nseg = length(nb)
#   cluster = rep(NA,nseg)
#   core = rep(FALSE,nseg)
#   clid = 1
#   queue = numeric(0)
#   for(s in 1:nseg) {
#     # cat("segment",s,"\n")
#     if(is.na(cluster[s])) {
#       ns = nb[[s]]
#       # cat("neighbors",ns,"\n")
#       if(length(ns)>=minlns) {
#         core[s] = TRUE
#         for(i in 1:length(ns)) {
#           cluster[ns[i]] = clid
#           if(ns[i]!=s) queue = c(queue, ns[i])
#         }  
#         # cat("queue",queue,"\n")
#         #Expand cluster
#         while(length(queue)>0) {
#           m = queue[1]
#           # cat("m",m,"\n")
#           nsm = nb[[m]]
#           # cat("m neighbors",nsm,"\n")
#           if(length(nsm)>=minlns) {
#             core[m] = TRUE
#             for(i in 1:length(nsm)) {
#               unclassified = is.na(cluster[nsm[i]])
#               if(unclassified || cluster[nsm[i]]==-1) cluster[nsm[i]] = clid
#               if(unclassified) queue = c(queue, nsm[i])
#             }
#           }
#           queue = queue[-1]
#         }
#         clid = clid + 1
#       } else {
#         cluster[s] = -1 #Noise
#       }
#     }
#   }
#   cluster[cluster==-1] = NA
#   
#   res<-list(pointDistances = d,
#             params = list(nt = nt, eps=eps, minlns=minlns, distance.type = distance.type),
#             segmentDistances=lsd,
#             cluster=cluster,
#             core=core,
#             representativeTrajectories = .representativeTrajectories(lsd,cluster, core, eps))
#   class(res)<-c("list","trajclust")
#   return(res)
# }
# 
# .plot.trajclust<-function(object, cmdscale.add=TRUE, draw.representatives = TRUE, exclude.unclassified=FALSE, ...) {
#   
#   D_points = object$pointDistances
#   D_segments = object$segmentDistances$Dseg
#   nt = object$params$nt
#   cluster = object$cluster
#   core = object$core
#   nseg = length(cluster)
#   nobj=nseg/(nt-1)
#   ninis = numeric(nseg)
#   nfins = numeric(nseg)
#   
#   for(i in 1:nseg) {
#     ninis[i] = i
#     nfins[i] = i+nobj
#   }
#   # print(cbind(ninis,nfins))
#   if(exclude.unclassified) {
#     rem = is.na(cluster)
#     # print(which(rem))
#     D_segments = as.dist(as.matrix(D_segments)[!rem,!rem])
#     ids = which(rem)
#     ids = unique(c(ids, ids+nobj))
#     # print(ids)
#     # print(dim(as.matrix(D_points)))
#     D_points = as.dist(as.matrix(D_points)[-ids,-ids])
#     cluster = cluster[!rem]
#     core = core[!rem]
#     for(i in 1:nseg) {
#       ninis[i] = ninis[i] - sum(ids<ninis[i])
#       nfins[i] = nfins[i] - sum(ids<nfins[i])
#     }
#     nseg = length(cluster)
#     ninis = ninis[!rem]
#     nfins = nfins[!rem]
#     # print(cbind(ninis,nfins))
#     # print(dim(as.matrix(D_points)))
#   }
#   
#   clids = unique(cluster)
#   clids = clids[!is.na(clids)]
#   clids = sort(clids)
#   core.colors = rainbow(length(clids))
#   noncore.colors = rainbow(length(clids), alpha=0.5)
#   cmd<-cmdscale(D_points,eig=TRUE, k=2, add=cmdscale.add)
#   x<-cmd$points[,1]
#   y<-cmd$points[,2]
#   plot(x,y, type="n", asp=1, xlab=paste0("PCoA 1 (", round(100*cmd$eig[1]/sum(cmd$eig)),"%)"), 
#        ylab=paste0("PCoA 2 (", round(100*cmd$eig[2]/sum(cmd$eig)),"%)"),...)
#   
#   for(i in 1:nseg) {
#     nini = ninis[i]
#     nfin = nfins[i]
#     unclassified = TRUE
#     for(c in 1:length(clids)) {
#       if(unclassified) {
#         cl = clids[c]
#         clnoncore =  cluster[i]==cl && !core[i] 
#         clcore =  cluster[i]==cl && core[i] 
#         if(is.na(clnoncore)) clnoncore = FALSE
#         if(is.na(clcore)) clcore = FALSE
#         if(clnoncore) arrows(x[nini],y[nini],x[nfin],y[nfin], col=noncore.colors[c], length=0.1, lwd=2)
#         else if(clcore) arrows(x[nini],y[nini],x[nfin],y[nfin], col=core.colors[c], length=0.1, lwd=2)
#         unclassified = !(clnoncore || clcore)
#       }
#     }
#     if(unclassified) arrows(x[nini],y[nini],x[nfin],y[nfin], col="black", length=0.05, lwd=1)
#     # text(x[ni93], y[ni93], labels = plotcodes[which(sel1)[i]], pos=3, cex=0.7)
#   }
#   if(draw.representatives) {
#     for(i in 1:length(clids)) {
#       trajnames = object$representativeTrajectories[[i]]$segments
#       # print(trajnames)
#       xl = numeric(0)
#       yl = numeric(0)
#       for(j in 1:length(trajnames)) {
#         k = which(row.names(as.matrix(D_segments))==trajnames[j])
#         xl = c(xl,x[ninis[k]], x[nfins[k]])
#         yl = c(yl,y[ninis[k]], y[nfins[k]])
#       }
#       lines(xl, yl, col="black", lwd=2)
#     }
#   }
#   if(length(clids)>1) {
#     legend("bottomright", legend = paste("cluster",clids), bty="n", col=core.colors, lty=1, lwd=2)
#   }
# }