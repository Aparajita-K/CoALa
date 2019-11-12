SimilarityFromFile<-function(W)
{   
     n=dim(W)[1]
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     AN=DH%*%W%*%DH
     I=diag(n)
     L=I+AN
     return(list(L=L,W=AN))
}
GraphLaplacian<-function(Dmat,rbfsg=1000,mod=1)
{
     n=dim(Dmat)[1]
     W=matrix(0,n,n)
     for(i in 1:n)
     {
        for(j in 1:i)
        {
             W[i,j]=W[j,i]=exp(-sum((Dmat[i,]-Dmat[j,])^2)/(2*rbfsg*rbfsg))
        }
     }
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     AN=DH%*%W%*%DH
     I=diag(n)
     L=I+AN
     return(list(L=L,W=AN))
}

gramschmidt <- function(x) {
  x <- as.matrix(x)
  n <- ncol(x)
  m <- nrow(x)
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)
   
  for (j in 1:n) {
    v = x[,j] 
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] 
        v <- v - r[i,j] * q[,i] 
      }      
    }
    r[j,j] <- sqrt(sum(v^2))
    q[,j] <- v / r[j,j]
  }
  qrcomp <- list(Q=q, R=r)
  return(qrcomp)
}
rowNormalize<-function(Data)
{
	n=dim(Data)[1]
	for(i in 1:n)
	{	
		if(sum(Data[i,]^2)>0)
	##To perform row normalization on the eigenvectors of the joint subspace, please uncomment the next line. 
	#Although row normalization is not recommended in the paper, so it is kept commented in the next line.
			Data[i,]=Data[i,]#/sqrt(sum(Data[i,]^2))
	}
	return(Data)
}


ConvexLaplacian<-function(Data,K,rank=NULL,alpha=NULL,beta=NULL,simFromFile=FALSE,mod=NULL)
{
	cat("\014")  
    M=length(Data)
    if(is.null(rank))
    {
		startRk=ceiling(K/M)
		maxRk=50 	#Set the maximum rank r to vary upto
    }
	else
		startRk=maxRk=rank
    step=1
    nkmeans=30
    rkSeq=seq(startRk,maxRk,step)
    SilOpt=-1
    VrOpt=NULL
    PirOpt=NULL
	show=5
	nstart=100
	ShiLm=list()
	ShiLmperp=list()
	Datas=list()
	Rel=vector(mode="numeric",length=M)
	Alpha=vector(mode="numeric",length=M)
	order=vector(mode="numeric",length=M)
	Dist=vector(mode="numeric",length=M)
	sigma=vector(mode="numeric",length=M)
	n=dim(Data[[1]])[1]
    for(rkFrac in rkSeq)
    {	 
    	
        nsc=rkFrac
        evind=seq(1,nsc)
        for(m in 1:M)
        {
        	tstartm=proc.time()
	        if(simFromFile==TRUE)
	           Lp=SimilarityFromFile(Data[[m]])
	        else
	        {
    	        Datas[[m]]=scale(Data[[m]],center=T,scale=F)
	            sgFrac=0.5
	            Dist[m]=max(as.numeric(dist(Datas[[m]])))
	            sigma[m]=sgFrac*Dist[m]
                Lp=GraphLaplacian(Datas[[m]],rbfsg=sigma[m],mod=m)
            }
            Lm=Lp$L
            evi=eigen(Lm)
            evi$values[ abs(evi$values)<1e-10 ] <- 0
            evi$values<-round(evi$values,digits=10)
            evind=which(evi$values!=2)[1:nsc]
            ln=length(evi$values)
            Um=evi$vectors[,evind,drop=FALSE]
            Dm=evi$values[evind]
            Dmat=evi$vectors[,2,drop=FALSE]
            km=kmeans(Dmat,2,iter.max=100,nstart=nstart)
			Rel[m]=(abs(sum(evi$values[2:2])))*silhouette(Dmat,km$cluster)
	        ShiLm[[m]]=list(U=Um,D=Dm,L=Lm)
        }
        order=sort(Rel,index.return=TRUE,decreasing=TRUE)$ix
        if(is.null(beta))
        	beta=1.25
        alphaW=rep(0,M)    
        if(is.null(alpha)==FALSE)
        	alphaW=alpha
        else
        {
        	for(m in 1:M)
	        	alphaW[order[m]]=Rel[order[m]]*(1/(beta^m)) 
        }
	    Alpha=alphaW/sum(alphaW)
	    if(is.null(mod))
        	mod=seq(1,M)
        if(rkFrac==rkSeq[1])
        {
        	for(m in 1:M)
        	{
	            cat("\nModality=",mod[m],":")
		        cat("\n\t\tRelevance=",Rel[m])
		        wh=which(order==m)
		        cat("\n\t\tRelevance Based Position in Ordering=",wh)
		        cat("\n\t\tAlpha=",Alpha[m],"\n\n")   
        	}  
        	cat("Optimizing Rank r:")      	            
        }
        Gamma<-list()
        Gamma[[1]]=as.matrix(ShiLm[[1]]$U)
        Basis=as.matrix(Gamma[[1]])
        for(m in 2:M)
        {
	        Sm=as.matrix(t(Basis)%*%ShiLm[[m]]$U)
	        Pm=as.matrix(Basis%*%Sm)
	        Qm=ShiLm[[m]]$U-Pm
	        Gm=gramschmidt(Qm)$Q
	        Gamma[[m]]=Gm
	        Basis=cbind(Basis,Gm)
        }
        HmList<-list()
        for(m in 1:M)
        {
	        HmTemp=matrix(0,M*nsc,M*nsc)
	        for(p in 1:M)
	        {
	            for(q in 1:M)
	            {
	                if(nsc==1)
	                    Blockpq=t(Gamma[[p]])%*%ShiLm[[m]]$U*ShiLm[[m]]$D%*%t(ShiLm[[m]]$U)%*%Gamma[[q]]
	                else
	                    Blockpq=t(Gamma[[p]])%*%ShiLm[[m]]$U%*%diag(ShiLm[[m]]$D)%*%t(ShiLm[[m]]$U)%*%Gamma[[q]]
	                Blockpq[ abs(Blockpq)<1e-10 ] <- 0
		            rs=(p-1)*nsc+1
		            rl=(p-1)*nsc+nsc
			        cs=(q-1)*nsc+1
			        cl=(q-1)*nsc+nsc
			        HmTemp[rs:rl,cs:cl]=Blockpq
		        }
	        }
	        HmList[[m]]=HmTemp
        }     
            tkind=seq(1,K)
            H=matrix(0,M*nsc,M*nsc)
            for(m in 1:M)
                H=H+Alpha[m]*HmList[[m]]
            eigH=eigen(H)
            R=eigH$vectors[,tkind,drop=FALSE]
            Pik=eigH$values[tkind]
            Vk=Basis%*%R
            Lkstar=Vk
            for(i in 1:ncol(Lkstar))
            {
            	nrm=sqrt(sum(Lkstar[,i]^2))
            	Lkstar[,i]=Lkstar[,i]/nrm
            }
            Lkstar=rowNormalize(Lkstar)
            km=kmeans(Lkstar,K,iter.max=100,nstart=nstart)$cluster
            RankSil=silhouette(Lkstar,km)
            cat("\nAt rank r=",nsc,"  Silhoutte index=",RankSil)
            if(RankSil>SilOpt)
            {
            	VrOpt=Lkstar
            	PirOpt=Pik
            	SilOpt=RankSil
            	rStar=nsc
            }
    }
    cat("\n\nOptimal rank r*=",rStar)
    cat("\nEigenvectors of joint Laplacian Lr* written to file: LrStar_Eigenvectors.txt")
    write.table(VrOpt,row.names=FALSE,col.names=FALSE,quote=FALSE,file="LrStar_Eigenvectors.txt")
    cat("\nEigenvalues of joint Laplacian Lr* written to file: LrStar_Eigenvalues.txt")
    cat("Eigenvalues Lr*=",PirOpt,file="LrStar_Eigenvalues.txt")
	return(list(U=VrOpt,D=PirOpt,Sil=SilOpt))
}



#Internal Indices
#Silhouette
sumsq<-function(v1,v2)
{
	eqd=sum((v1-v2)^2)
	return(eqd)
}
proximity<-function(Dmat)
{
	n=dim(Dmat)[1]
	proxim=matrix(0,n,n)
	for(i in 1:n)
	{
		for(j in 1:(i-1))
		{
			proxim[i,j]=sqrt(sumsq(Dmat[i,],Dmat[j,]))
			proxim[j,i]=proxim[i,j]
		}
	}
	return(proxim)
}
silhouette<-function(Dmat,kmclust)
{
    proxim=proximity(Dmat)
    K=max(kmclust)
	cl=seq(1:K)
	n=dim(Dmat)[1]
	s<-vector(mode="numeric",length=n)
	bvec<-vector(mode="numeric",length=K-1)
	for(i in 1:n)
	{
		own=kmclust[i]
		temp=proxim[i,]
		restown=setdiff(which(kmclust==own),i)
		if(length(restown)>0)
		{
			owndist=temp[restown]
			ai=mean(owndist)
		}
		else
			ai=0
		oc=setdiff(cl,own)
		for(j in 1:length(oc))
		{
			oth=oc[j]
			otherclust=which(kmclust==oth)
			othdist=temp[otherclust]
			bvec[j]=mean(othdist)
		}
		bi=min(bvec)
		s[i]=(bi-ai)/max(ai,bi)
	}
	clustsil=rep(0,K)
	for(i in 1:K)
	{
		temp=which(kmclust==i)
		clsil=s[temp]
		clustsil[i]=mean(clsil)
	}
    sil=mean(s)
	return(sil)
}
