#model based on Bailey et al. 2013 "Can Cooperative Management of Tuna Fisheries in the Western Pacific Solve the Growth Overfishing Problem?" Strategic Behavior and the Environment 3:31-66

#updated and modified by Ciara Willis 2019


library(ggplot2)
library(reshape)
library(gridExtra)
library(viridis)

#Model output----
output2=read.csv('output2.csv')

#Input parameters----
#species order: skipjack, yellowfin, bigeye
#gears order: longline, purse seine, handline, small scale surface

set.seed(5)

nsp = 3	#number of species
ro = c(100000000000,10000000000,5000000000) #Unfished age-1 recruits
h = c(0.90,0.75,0.75) #steepness of stock recruitment
kap = 4*h/(1-h) #Recruitment compensation ratio
a = c(0.0000086388,0.00002152,0.000019729)
b = c(3.2174,2.9369,3.0274)
linf = c(106,175,180) #max length
winf = a*linf^b #max weight
wmat = c(1.56,19,23.47) #weight at maturity
vbk = c(0.3105,0.392,0.31525) #von Bert metabolic coefficient
m = 1.5*vbk #Natural mortality rate
A1 = c(5,6,7)#max age
abar = 1.65/m #Age of 50% maturity
gbar = 0.1*abar #Std in maturity


ngr = 4 #number of fishing gears

#LL
ahatLL = c(40, 3, 3) #Age at 50% vulnerability. Big # for skipjack bc no catch
ghat = 0.25 #Std in vulnerbablity

#PS
ahatPS = c(1, 40, 40)

lh1p = c(30, 60, 40)
lh2p = c(80, 160, 120)
sd1p = c(5, 10, 10)
sd2p = c(1, 10, 10)

#HL
ahatHL = c(40, 4, 4) #Age at 50% vulnerability
ghat = 0.25 #Std in vulnerbablity

#SSS
lh1 = c(15, 20, 20)
lh2 = c(60, 60, 60)
sd1 = c(1, 5, 5)
sd2 = c(5, 20, 10)



q = list(c(0, 7.142857E-06*.5, 0, 7.142857E-08*.5), #catchability
         c(1.050000E-10*2, 3.230769E-06*.5, 1.615437E-08*.75, 1.615437E-08*.5),
         c(1.050000E-10*.75, 1.153846E-06*.5, 5.769981E-09*.75, 5.769981E-09*.5))


cost = list(c(0, 22000, 0, 100), #cost per unit effort
            c(1, 0, 100, 0),
            c(0, 0, 0, 0))


p=list(c(0, 1600, 0, 1300), #landed value (price per t)
       c(8000, 2300, 7700, 2000),
       c(11000, 2600, 10700, 2300))

subs = list(c(0.2105273, 0.2139164, 0.02659722, 0.02659722), #subsidy intensity
            c(0.2882288, 0.1688027, 0.0260168, 0.0260168),  
            c(0.2744555, 0.1924171, 0.02730034, 0.02730034))

msy= c(1891600, 586400, 159020) #max sust yield from Brouwer 2018

theta=list(ro=ro, kap=kap, winf=winf, vbk=vbk, m=m, A1=A1,
           abar=abar, gbar=gbar, ngr=ngr, q=q,
           ahatLL=ahatLL, ahatPS=ahatPS, ahatHL=ahatHL, ghat=ghat, sd1p=sd1p, sd2p=sd2p, lh1p=lh1p, lh2p=lh2p, sd1=sd1, sd2=sd2, lh1=lh1, lh2=lh2,
           p=p, cost=cost, subs=subs, msy=msy)


#Find optimal efforts by gear ----
equil.age <- function(theta,eff=rep(0,ngr))
{
  with(theta, {
    fn.equil <- function(ii) #ii is species index
    { 
      age = 0:A1[ii]
      la = linf[ii]*(1-exp(-vbk[ii]*(age+0.5))) #length at age
      # plot(la)
      wa = winf[ii]*(1-exp(-vbk[ii]*(age+0.5)))^3
      ma = m[ii]*(linf[ii]/la) #mortality at age
      fec = wa-wmat[ii]; fec[fec<0]=0 #fecundity at age
      
      A=length(age)
      lx=vector() #survivorship
      lx[1]=1; for(i in 2:A)
        lx[i] = lx[i-1]*exp(-ma[i-1])
      lx[A]=lx[A]/(1-exp(-ma[A]))
      phie = sum(lx*fec)	#Unfished eggs per recruit
      
      #fished conditions
      ##selectivities:
      valLL = sapply(age,plogis,ahatLL[ii],ghat) #This is logistic selectivity for longline
      valPS = sapply(age,plogis,ahatPS[ii],ghat) #logistic selectivity for PS
      vadPS = (1/(1+exp(-(1/sd1p[ii])*(la-lh1p[ii]))))*(1/(1+exp((1/sd2p[ii])*(la-lh2p[ii])))) #dome selectivity for PS
      valHL = sapply(age,plogis,ahatHL[ii],ghat) #logistic selectivity for HL
      vadSSS = (1/(1+exp(-(1/sd1[ii])*(la-lh1[ii]))))*(1/(1+exp((1/sd2[ii])*(la-lh2[ii])))) 	# dome selectivity for SSS
      
      if (ii==1) {
        va = rbind(valLL, valPS, valHL, vadSSS) # skipjack
      }
      
      if (ii==2 | ii==3) {
        va = rbind(valLL, vadPS, valHL, vadSSS) #yellowfin, bigeye
      }
      
      rownames(va)=c("LL", "PS", "HL", "SSS")
      
      fe = q[[ii]]*eff #catchability * effort
      fa = fe*va  #age & gear specific fishing mortality
      
      za = ma + colSums(fa)	# Total mortality by age
      qa = va/za*(1-exp(-za))	# Mortality by gear type
      
      lz=vector() #survivorship
      lz[1]=1; for(i in 2:A)
        lz[i]=lz[i-1]*exp(-za[i-1])
      lz[A]=lz[A]/(1-exp(-za[A]))
      # plot(lz)
      
      phif = sum(lz*fec)		#Fished eggs per recruit
      phiq = rowSums(lz*wa*qa)	#vulnerable biomass per recruit for each gear
      
      re = max(0,ro[ii]*(kap[ii]-phie/phif)/(kap[ii]-1))
      sb = re*phif/1000		#Spawning biomass
      vb = re*phiq/1000 #vulnerable biomass per gear
      ye = fe*re*phiq/1000		#equilibrium yield for each gear type.
      yet = sum(ye)
      
      #limit model outputs so all species' catch is below 80% MSY
      if (yet > (msy[ii]*0.8)) {
        ye = rep(NA,3)
        
      }
      
      revs = ye*p[[ii]] #revenue
      costs = cost[[ii]]*eff #costs
      PR = revs-costs #private rent
      SR = revs*(1-subs[[ii]])-costs #social rent
      
      
      results=(c(sb=sb,vb=vb,ye=ye,revs=revs,costs=costs,PR=PR,SR=SR))
      
      return(results)
    }
    xx=sapply(1:nsp,fn.equil, simplify=FALSE)
    LLpi=xx[[1]][18]+xx[[2]][18]+xx[[3]][18] #LL PR
    PSpi=xx[[1]][19]+xx[[2]][19]+xx[[3]][19]
    HLpi=xx[[1]][20]+xx[[2]][20]+xx[[3]][20]
    SSSpi=xx[[1]][21]+xx[[2]][21]+xx[[3]][21]
    nc.PR.vector=c(LLpi,PSpi, HLpi, SSSpi)
    nc.PR.sum=LLpi+PSpi+HLpi+SSSpi
    
    LLsri=xx[[1]][22]+xx[[2]][22]+xx[[3]][22] #LL SR
    PSsri=xx[[1]][23]+xx[[2]][23]+xx[[3]][23]
    HLsri=xx[[1]][24]+xx[[2]][24]+xx[[3]][24]
    SSSsri=xx[[1]][25]+xx[[2]][25]+xx[[3]][25]
    nc.SR.sum=LLsri+PSsri+HLsri+SSSsri
    
    LLyi=xx[[1]][6]+xx[[2]][6]+xx[[3]][6] #LL yield
    PSyi=xx[[1]][7]+xx[[2]][7]+xx[[3]][7]
    HLyi=xx[[1]][8]+xx[[2]][8]+xx[[3]][8]
    SSSyi=xx[[1]][9]+xx[[2]][9]+xx[[3]][9]
    nc.ye.sum=LLyi+PSyi+HLyi+SSSyi
    
    psy=c(nc.PR.sum, nc.SR.sum, nc.ye.sum)
    
    
    ifelse(test = nc.PR.vector[1] < 0 | nc.PR.vector[2] < 0 | nc.PR.vector[3] < 0 | nc.PR.vector[4] < 0 | is.na(nc.PR.vector), return(rep(NA,3)),
           return(psy)
    )
    
    
    
  })
}



#Run model----
VL=20 # number of effort levels per gear

#input effort levels per gear
fe1=seq(1,6e8, length=VL) #LL
fe2=seq(1,4.5e4,length=VL) #PS
fe3=seq(1,1.9e7,length=VL) #HL
fe4=seq(1,1e7, length=VL) #SSS


eff=expand.grid(fe1, fe2, fe3, fe4) #all possible combinations of gear efforts
names(eff)=c("LL", "PS", "HL", "SSS")

pi.fn=function(eff) equil.age(theta, eff)
PROF=apply(eff,1,pi.fn) #apply a fn (pi.fn) to matrix (eff) over rows of the matrix (1)


PROF=as.data.frame(t(PROF))
names(PROF)=c("PR", "SR", "Yield")

output=data.frame(eff,PROF)
output2=na.omit(output)

#Index optimal scenarios----
# plot(output$PR)
# plot(output$SR)
# plot(output$Yield)

pr.optf=output2[which(output2$PR==max(output2$PR)),]
sr.optf=output2[which(output2$SR==max(output2$SR)),]
ye.optf=output2[which(output2$Yield==max(output2$Yield)),]

pr.optf
sr.optf
ye.optf

#Optimal scenarios detailed results----
eff=as.numeric(pr.optf[2:5])
# eff=sr.optf
# eff=ye.optf

xx=sapply(1:nsp,fn.equil, simplify=FALSE)

xx

#landed values
lv=c(xx[[1]][10]+xx[[2]][10]+xx[[3]][10],
     xx[[1]][11]+xx[[2]][11]+xx[[3]][11],
     xx[[1]][12]+xx[[2]][12]+xx[[3]][12],
     xx[[1]][13]+xx[[2]][13]+xx[[3]][13])
lv
sum(lv)
lv.actual=c(1.46e9-0.34e9, 3.4e9) #LL total - albacore
as.integer((lv[1:2]/lv.actual)*100) #percent of current catch for LL and PS

#fishing costs
costs.out=c(xx[[1]][14]+xx[[2]][14]+xx[[3]][14],
            xx[[1]][15]+xx[[2]][15]+xx[[3]][15],
            xx[[1]][16]+xx[[2]][16]+xx[[3]][16],
            xx[[1]][17]+xx[[2]][17]+xx[[3]][17])
names(costs.out)=c("LL", "PS", "HL", "SSS")
costs.out

#Private rent
LLpi=xx[[1]][18]+xx[[2]][18]+xx[[3]][18]
PSpi=xx[[1]][19]+xx[[2]][19]+xx[[3]][19]
HLpi=xx[[1]][20]+xx[[2]][20]+xx[[3]][20]
SSSpi=xx[[1]][21]+xx[[2]][21]+xx[[3]][21]
nc.PR.vector=c(LLpi,PSpi,HLpi, SSSpi)
nc.PR.sum=LLpi+PSpi+HLpi+SSSpi
nc.PR.vector
nc.PR.sum

#Social rent
LLsri=xx[[1]][22]+xx[[2]][22]+xx[[3]][22]
PSsri=xx[[1]][23]+xx[[2]][23]+xx[[3]][23]
HLsri=xx[[1]][24]+xx[[2]][24]+xx[[3]][24]
SSSsri=xx[[1]][25]+xx[[2]][25]+xx[[3]][25]
nc.SR.vector=c(LLsri,PSsri,HLsri, SSSsri)
nc.SR.sum=LLsri+PSsri+HLsri+SSSsri
nc.SR.vector
nc.SR.sum


#Yields
LLy2=c(xx[[1]][6], xx[[2]][6], xx[[3]][6])
PSy2=c(xx[[1]][7],xx[[2]][7],xx[[3]][7])
HLy2=c(xx[[1]][8],xx[[2]][8],xx[[3]][8])
SSSy2=c(xx[[1]][9],xx[[2]][9],xx[[3]][9])
ye2=rbind(LLy2, PSy2, HLy2, SSSy2)
colnames(ye2)[1:3]=c("SJ", "YF", "BE")

ye2=as.data.frame(ye2)
ye2$gear=c("LL", "PS", "HL", "SSS")
ye2

ye2m=melt(ye2, id.vars = "gear")
names(ye2m)=c("gear", "species", "catch")
ye2m

#catch by gear
yeg = c(sum(ye2[1,1:3]), sum(ye2[2,1:3]), sum(ye2[3,1:3]), sum(ye2[4,1:3]))
as.integer(yeg)
yegactual = c(1.45e5, 1.65e6, 514000, 514000) #SSF estimate from Gillet for WCPO
as.integer((yeg/yegactual)*100)

#catch by species
yes = colSums(ye2[,1:3])
as.integer(yes)
yesactual = c(1.6e6, 550000, 150000)
as.integer((yes/yesactual)*100) #as percent current catch
as.integer((yes/msy)*100) #as percent of MSY

#jobs
jobs.vector = as.integer(c(eff[1]*0.0001125, eff[2]*0.1875, eff[3]*2, eff[4]*2))
jobs.sum = as.integer(sum(jobs.vector))
jobs.vector
jobs.sum



#Figure 1 vulnerability----
txts=10
txtl=15

va4 = data.frame(LL=vector(),PS=vector(),HL=vector(),SSS=vector(),sp=vector(),age=vector(),size=vector())

for (ii in 1:nsp) {
  age = 0:A1[ii]
  la = linf[ii]*(1-exp(-vbk[ii]*(age+0.5)))
  wa = winf[ii]*(1-exp(-vbk[ii]*(age+0.5)))^3
  ma = m[ii]*(linf[ii]/la)	
  fec = wa-wmat[ii]; fec[fec<0]=0
  
  A=length(age) 
  lx=vector() 
  lx[1]=1; for(i in 2:A)
    lx[i] = lx[i-1]*exp(-ma[i-1])
  lx[A]=lx[A]/(1-exp(-ma[A]))
  phie = sum(lx*fec)
  
  valLL = sapply(age,plogis,ahatLL[ii],ghat)
  valPS = sapply(age,plogis,ahatPS[ii],ghat)
  vadPS = (1/(1+exp(-(1/sd1p[ii])*(la-lh1p[ii]))))*(1/(1+exp((1/sd2p[ii])*(la-lh2p[ii]))))
  valHL = sapply(age,plogis,ahatHL[ii],ghat) #log for HL
  vadSSS = (1/(1+exp(-(1/sd1[ii])*(la-lh1[ii]))))*(1/(1+exp((1/sd2[ii])*(la-lh2[ii])))) 
  
  if (ii==1) {
    va = rbind(valLL, valPS, valHL, vadSSS)
  }
  
  if (ii==2 | ii==3) {
    va = rbind(valLL, vadPS, valHL, vadSSS)
  }
  
  rownames(va)=c("LL", "PS", "HL", "SSS")
  
  vat=t(va)
  vat=as.data.frame(vat)
  vat$sp=as.factor(ifelse(ii==1, "Skipjack", ifelse(ii==2, "Yellowfin", "Bigeye")))
  vat$age=age
  vat$size=la
  
  va4=rbind(va4, vat)
  
}


va4m=melt(va4, id.vars = c("age", "size", "sp"))

ggplot(va4m, aes(x=age, y=value))+
  facet_grid(.~sp, scales = "free")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA))+ #publication
  theme(axis.text.x = element_text(colour = "black", size = txts),
        axis.title.x = element_text(size=txtl), 
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title.y = element_text(size=txtl),
        legend.text = element_text(size=txtl),
        legend.title = element_blank(), 
        legend.key = element_rect(colour="transparent", fill="white"),
        legend.position = "top")+
  theme(strip.text = element_text(size=txtl))+
  geom_line(aes(colour=variable, linetype=variable), size=1.5)+
  scale_color_manual(values = c("orangered2", "gold2", "deepskyblue", "midnightblue"))+
  labs(x="Age", y="Vulnerability")
# ggsave("Fig 1 Vulnerability.png", height = 4, width = 8.5, unit="in")
# ggsave("Figure_1.pdf", height = 4*18/8.5, width = 18, unit="cm")
# ggsave("Fig 1.png", height = 4*18/8.5, width = 18, unit="cm") #poster

#Figure 2 rents and catch by scenario-----
txts=10
txtl=15

eff=list(c(as.numeric(pr.optf[2:5])),
         c(as.numeric(sr.optf[2:5])),
         c(as.numeric(ye.optf[2:5])))

rents3 = data.frame(scen=vector(), type = vector(), variable = vector(), value = vector())


yield3 = data.frame(scen=vector(), gear = vector(), variable = vector(), value = vector())

for (j in 1:length(eff)) {
  for (ii in 1:nsp) {
    
    
    fn.equil <- function(ii)
    {
      age = 0:A1[ii]
      la = linf[ii]*(1-exp(-vbk[ii]*(age+0.5)))
      wa = winf[ii]*(1-exp(-vbk[ii]*(age+0.5)))^3
      ma = m[ii]*(linf[ii]/la)
      fec = wa-wmat[ii]; fec[fec<0]=0 
      
      A=length(age)
      lx=vector()
      lx[1]=1; for(i in 2:A)
        lx[i] = lx[i-1]*exp(-ma[i-1])
      lx[A]=lx[A]/(1-exp(-ma[A]))
      phie = sum(lx*fec)

      valLL = sapply(age,plogis,ahatLL[ii],ghat)
      valPS = sapply(age,plogis,ahatPS[ii],ghat)
      vadPS = (1/(1+exp(-(1/sd1p[ii])*(la-lh1p[ii]))))*(1/(1+exp((1/sd2p[ii])*(la-lh2p[ii]))))
      valHL = sapply(age,plogis,ahatHL[ii],ghat)
      vadSSS = (1/(1+exp(-(1/sd1[ii])*(la-lh1[ii]))))*(1/(1+exp((1/sd2[ii])*(la-lh2[ii]))))
      
      if (ii==1) {
        va = rbind(valLL, valPS, valHL, vadSSS)
      }
      
      if (ii==2 | ii==3) {
        va = rbind(valLL, vadPS, valHL, vadSSS)
      }
      
      rownames(va)=c("LL", "PS", "HL", "SSS")
      
      fe = q[[ii]]*eff[[j]]
      fa = fe*va
      
      za = ma + colSums(fa)
      qa = va/za*(1-exp(-za))
      
      lz=vector()
      lz[1]=1; for(i in 2:A)
        lz[i]=lz[i-1]*exp(-za[i-1])
      lz[A]=lz[A]/(1-exp(-za[A]))
      
      phif = sum(lz*fec)
      phiq = rowSums(lz*wa*qa)
      
      re = max(0,ro[ii]*(kap[ii]-phie/phif)/(kap[ii]-1))
      sb = re*phif/1000
      vb = re*phiq/1000
      ye = fe*re*phiq/1000
      yet = sum(ye)
      
      revs = ye*p[[ii]]
      costs = cost[[ii]]*eff[[j]]
      PR = revs-costs
      SR = revs*(1-subs[[ii]])-costs
      
      results=(c(sb=sb,vb=vb,ye=ye,revs=revs,costs=costs,PR=PR,SR=SR))
      
      return(results)
    }
    xx=sapply(1:nsp,fn.equil, simplify=FALSE)
    
    LLpi=xx[[1]][18]+xx[[2]][18]+xx[[3]][18]
    PSpi=xx[[1]][19]+xx[[2]][19]+xx[[3]][19]
    HLpi=xx[[1]][20]+xx[[2]][20]+xx[[3]][20]
    SSSpi=xx[[1]][21]+xx[[2]][21]+xx[[3]][21]
    nc.PR.vector=c(LLpi,PSpi,HLpi, SSSpi)
    
    
    LLsri=xx[[1]][22]+xx[[2]][22]+xx[[3]][22]
    PSsri=xx[[1]][23]+xx[[2]][23]+xx[[3]][23]
    HLsri=xx[[1]][24]+xx[[2]][24]+xx[[3]][24]
    SSSsri=xx[[1]][25]+xx[[2]][25]+xx[[3]][25]
    nc.SR.vector=c(LLsri,PSsri,HLsri, SSSsri)
    
    rents = rbind(nc.PR.vector, nc.SR.vector)
    rents = as.data.frame(rents)
    rents$type = c("PR", "SR")
    names(rents)[1:4]=c("LL", "PS", "HL", "SSS")
    
    rents.m = melt(rents, id.vars = "type")
    
    rents2 = data.frame(scen=as.factor(ifelse(j==1, "Max PR", ifelse(j==2, "Max SR", "Max Yield"))), rents.m)

    rents3=rbind(rents3, rents2)
  
    
    LLy2=c(xx[[1]][6], xx[[2]][6], xx[[3]][6])
    PSy2=c(xx[[1]][7],xx[[2]][7],xx[[3]][7])
    HLy2=c(xx[[1]][8],xx[[2]][8],xx[[3]][8])
    SSSy2=c(xx[[1]][9],xx[[2]][9],xx[[3]][9])
    nc.ye2.vector=rbind(LLy2, PSy2, HLy2, SSSy2)
    colnames(nc.ye2.vector)[1:3]=c("SJ", "YF", "BE")
    
    nc.ye2.vector=as.data.frame(nc.ye2.vector)
    nc.ye2.vector$gear=c("LL", "PS", "HL", "SSS")
    
    yield1=melt(nc.ye2.vector, id.vars = "gear")
    yield2 = data.frame(scen=as.factor(ifelse(j==1, "Max PR", ifelse(j==2, "Max SR", "Max Yield"))), yield1)
    yield3=rbind(yield3, yield2)
    
  }
  
}

yield3$gear=factor(yield3$gear, levels = c("LL", "PS", "HL", "SSS"))


rents=ggplot(rents3, aes(x=variable, y=value))+
  facet_grid(.~scen, scales = "free")+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = txts), 
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title.y = element_text(size=txtl),
        legend.title = element_blank(),
        legend.position = c(.72,.83),
        legend.background = element_blank(),
        legend.text = element_text(size=txtl))+
  theme(strip.text = element_text(size=txtl))+
  geom_col(aes(fill=type), position = "dodge")+
  scale_fill_manual(values = c("red3", "orange"))+
  labs(y="Rent (USD)")
rents

scaleFUN <- function(x) sprintf("%.e", x)

yields = ggplot(yield3, aes(x=gear, y=value))+
  facet_grid(.~scen, scales = "fixed")+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = txts), 
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title.y = element_text(size=txtl),
        legend.title = element_blank(),
        legend.position = c(.72,.76), #pdf
        # legend.position = c(.72,.785), #png
        legend.background = element_blank(),
        legend.text = element_text(size=txtl))+
  theme(strip.text = element_text(size=txtl))+
  geom_col(aes(fill=variable))+
  scale_fill_manual(values = c("dodgerblue3", "gold", "green3"))+
  scale_y_continuous(labels = scaleFUN)+
  labs(y="Catch (t)")
# yields

# grid.arrange(rents, yields, ncol=1)
g2=arrangeGrob(rents, yields, ncol=1)
# ggsave("Fig 2 Rent and Yield GS.png", g2, height = 11/2, width = 8.5, units = "in")
# ggsave("Figure_2.pdf", g2, height = (11/2)*18/8.5, width = 18, units = "cm")
# ggsave("Fig 2.png", g2, height = (11/2)*18/8.5, width = 18, units = "cm") #poster

#Figure 3 catch by age and gear by scenario----
txts=10
txtl=15

eff=list(c(as.numeric(pr.optf[2:5])),
         c(as.numeric(sr.optf[2:5])),
         c(as.numeric(ye.optf[2:5])))

ye4 = data.frame(scen=vector(), sp = vector(), age = vector(), gear = vector(), yield= vector())


for (j in 1:length(eff)) {
  for (ii in 1:nsp) {
    age = 0:A1[ii]
    la = linf[ii]*(1-exp(-vbk[ii]*(age+0.5)))
    wa = winf[ii]*(1-exp(-vbk[ii]*(age+0.5)))^3
    ma = m[ii]*(linf[ii]/la)
    fec = wa-wmat[ii]; fec[fec<0]=0
    
    A=length(age)
    lx=vector()
    lx[1]=1; for(i in 2:A)
      lx[i] = lx[i-1]*exp(-ma[i-1])
    lx[A]=lx[A]/(1-exp(-ma[A]))
    phie = sum(lx*fec)
    phie2 = lx*fec
    
    valLL = sapply(age,plogis,ahatLL[ii],ghat)
    valPS = sapply(age,plogis,ahatPS[ii],ghat)
    vadPS = (1/(1+exp(-(1/sd1p[ii])*(la-lh1p[ii]))))*(1/(1+exp((1/sd2p[ii])*(la-lh2p[ii]))))
    valHL = sapply(age,plogis,ahatHL[ii],ghat)
    vadSSS = (1/(1+exp(-(1/sd1[ii])*(la-lh1[ii]))))*(1/(1+exp((1/sd2[ii])*(la-lh2[ii])))) 
    
    if (ii==1) {
      va = rbind(valLL, valPS, valHL, vadSSS) 
    }
    
    if (ii==2 | ii==3) {
      va = rbind(valLL, vadPS, valHL, vadSSS) 
    }
    
    rownames(va)=c("LL", "PS", "HL", "SSS")
    
    fe = q[[ii]]*eff[[j]]
    fa = fe*va
    
    za = ma + colSums(fa)
    qa = va/za*(1-exp(-za))
    
    lz=vector()
    lz[1]=1; for(i in 2:A)
      lz[i]=lz[i-1]*exp(-za[i-1])
    lz[A]=lz[A]/(1-exp(-za[A]))
    
    phif = sum(lz*fec)
    phif2 = lz*fec
    phiq = rowSums(lz*wa*qa)	
    phiq2 = lz*wa*qa 
    
    re = max(0,ro[ii]*(kap[ii]-phie/phif)/(kap[ii]-1))
    
    sb = re*phif/1000	
    sb2 = re*phif2/1000
    vb = re*phiq/1000 
    ye = fe*re*phiq/1000
    ye2 = fe*re*phiq2/1000
    yet = sum(ye)
    
    
    ye3 = data.frame(scen=as.factor(ifelse(j==1, "Max PR", ifelse(j==2, "Max SR", "Max Yield"))),
                     sp=as.factor(ifelse(ii==1, "Skipjack", ifelse(ii==2, "Yellowfin", "Bigeye"))),
                     age = c(rep(0:(length(age)-1), 4)),
                     gear=c(rep("LL",length(age)), rep("PS", length(age)), rep("HL",length(age)), rep("SSS", length(age))),
                     yield=c(ye2[1,], ye2[2,], ye2[3,], ye2[4,]))
    
    ye3$gear=factor(ye3$gear, levels = c("LL", "PS", "HL", "SSS"))
    
    ye4 = rbind(ye4, ye3)
  }
  
}


ggplot(ye4, aes(x=age, y=yield))+
  facet_grid(scen~sp, scales = "free_x")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA))+ #publication
  theme(axis.text.x = element_text(colour = "black", size = txts),
        axis.title.x = element_text(size=txtl), 
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title.y = element_text(size=txtl),
        legend.text = element_text(size=txtl), legend.title = element_blank(),
        legend.position = c(.735,.917), #pdf
        legend.background = element_blank()
  )+
  theme(strip.text = element_text(size=txtl))+
  scale_y_continuous(labels = scales::scientific)+
  scale_fill_manual(values = c("orangered2", "gold2", "deepskyblue", "midnightblue"))+
  geom_col(aes(fill=gear))+
  labs(x="Age", y="Yield (t)")
# ggsave("Figure_3.pdf", height = 18, width = 18, units = "cm")
# ggsave("Fig 3.png", height = 18, width = 18, units = "cm") #poster



#Figure 4 PSY histogram----
scaleFUN1 <- function(x) sprintf("%.1e", x)

output2m=melt(output2[,6:8])

ggplot(output2m, aes(fill=variable))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(colour = "black", size = txts, angle = 15, vjust = .8),
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title.y = element_text(size=txtl),
        axis.title.x = element_blank(), #THIS PLOT ONLY
        legend.title = element_text(size=txtl),
        legend.background = element_blank(),
        legend.text = element_text(size=txts),
        legend.position = "none")+
  theme(strip.text = element_text(size=txtl))+
  geom_histogram(aes(x=value), colour="white")+
  scale_fill_manual(values=rep("#440154FF",3))+
  scale_x_continuous(label=scaleFUN1)+
  facet_grid(.~variable, scales = "free_x")+
  labs(y="Frequency")
# ggsave("Figure_4.pdf", height = ((11/2)*18/8.5)/2, width = 18, units = "cm")
# ggsave("Fig 4.png", height = ((11/2)*18/8.5)/2, width = 18, units = "cm")

#percent scenarios above 80% max
100*dim(output2[which(output2$PR>max(output2$PR)*0.8),])[1]/dim(output2)[1] #26%
100*dim(output2[which(output2$PR>max(output2$PR)*0.9),])[1]/dim(output2)[1] #4%

100*dim(output2[which(output2$SR>max(output2$SR)*0.8),])[1]/dim(output2)[1] #29%
100*dim(output2[which(output2$SR>max(output2$SR)*0.9),])[1]/dim(output2)[1] #5%

100*dim(output2[which(output2$Yield>max(output2$Yield)*0.8),])[1]/dim(output2)[1] #40%
100*dim(output2[which(output2$Yield>max(output2$Yield)*0.9),])[1]/dim(output2)[1] #13%


#Figure 5 2D PSY histograms----
cor(output2$PR, output2$SR, method = "spearman")
cor(output2$PR, output2$Yield, method = "spearman")
cor(output2$SR, output2$Yield, method = "spearman")

psh = ggplot(output2, aes(x=PR, y=SR))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black", size = txts-1),
        axis.title = element_text(size=txtl),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=txts-1),
        legend.position = c(.105,.71))+ #pdf
  geom_hex()+
  geom_abline(intercept = c(0,0), slope = 1, size = 1, colour="gray60", linetype = 2)+
  geom_smooth(method = "lm", size = 1, colour="gray70", formula = y ~ 0 + x)+
  annotate("text", x=3.3e9, y=0.1e9, label = "SR = 0.80*PR", size=4)+
  scale_fill_viridis()
# psh

pslm=lm(SR~0+PR, data = output2)
pslm
summary(pslm)


# hist((output2$SR/output2$PR))
# min(output2$SR/output2$PR) #66%
# max(output2$SR/output2$PR) #96%
# mean(output2$SR/output2$PR) #81


scaleFUN1 <- function(x) sprintf("%.1e", x)
yph = ggplot(output2, aes(x=Yield, y=PR))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        # axis.text.x = element_text(colour = "black", size = txts, angle = 45, vjust = .7),
        axis.text = element_text(colour = "black", size = txts-1),
        axis.title = element_text(size=txtl),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=txts-1),
        legend.position = c(.105,.71))+ #pdf
  geom_hex()+
  scale_fill_viridis()+
  scale_x_continuous(labels = scaleFUN1)
# yph

ysh = ggplot(output2, aes(x=Yield, y=SR))+
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        # axis.text.x = element_text(colour = "black", size = txts, angle = 30, vjust = .7),
        axis.text = element_text(colour = "black", size = txts-1),
        axis.title = element_text(size=txtl),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=txts-1),
        # legend.position = c(.1,.74))+ #png
        legend.position = c(.105,.71))+ #pdf
  geom_hex()+
  scale_fill_viridis()+
  scale_x_continuous(labels = scaleFUN1)
# ysh


# grid.arrange(psh, yph, ysh, ncol=1)
g=arrangeGrob(psh, yph, ysh, ncol=1)
# ggsave("Figure_5.pdf", g, height = 9*8/3.5, width = 8, units = "cm")
# ggsave("Fig 5.png", g, height = 9*8/3.5, width = 8, units = "cm")


#Figure 6 3D PSY----
# library(plotly)
# plot_ly(output2, x = ~PR, y = ~SR, z = ~Yield, color = ~Yield) %>% add_markers()


#Figure S1 plot PSY over all model results----
#NB that fig S1 uses a reduced version of the model with VL=10 to make it easier to view these plots

output2.10=read.csv('output2.10.csv')

txtl=12
txts=7


scaleFUN <- function(x) sprintf("%.e", x)
ggplot(output2.10, aes(x=LL, y=PS, fill=SR))+ #change variable here
  theme(axis.text.x = element_text(colour = "black", size = txts, angle = 90), 
        axis.text.y = element_text(colour = "black", size = txts),
        axis.title = element_text(size=txtl),
        legend.title = element_text(size=txtl),
        legend.background = element_blank(),
        legend.text = element_text(size=txts))+
  theme(strip.text = element_text(size=txts))+
  theme(panel.background = element_rect(fill="gray90", colour="gray80"),
        panel.grid.major = element_line(linetype = "solid", colour="white", size=0.25))+
    geom_tile()+
  scale_fill_viridis(labels = scaleFUN)+ #PR and SR
  # scale_fill_viridis(labels = scaleFUN, breaks=c(1e6, 2e6))+ #Yield only
  facet_grid(as.integer(HL)~SSS)+
  
#circle top value for PR
  # geom_point(size=4, shape=21, colour="#440154FF", fill="transparent", stroke=1,
  #            aes(x=output2.10[which(output2.10$PR==max(output2.10$PR)),]$LL, y=output2.10[which(output2.10$PR==max(output2.10$PR)),]$PS),
  #            data= subset(output2.10, HL == output2.10[which(output2.10$PR==max(output2.10$PR)),]$HL & SSS == output2.10[which(output2.10$PR==max(output2.10$PR)),]$SSS) )

#circle top value for SR 
  geom_point(size=4, shape=21, colour="#440154FF", fill="transparent", stroke=1,
             aes(x=output2.10[which(output2.10$SR==max(output2.10$SR)),]$LL, y=output2.10[which(output2.10$SR==max(output2.10$SR)),]$PS),
             data= subset(output2.10, HL == output2.10[which(output2.10$SR==max(output2.10$SR)),]$HL & SSS == output2.10[which(output2.10$SR==max(output2.10$SR)),]$SSS) )
  
#circle top value for Yield
  # geom_point(size=4, shape=21, colour="#440154FF", fill="transparent", stroke=1,
  #            aes(x=output2.10[which(output2.10$Yield==max(output2.10$Yield)),]$LL, y=output2.10[which(output2.10$Yield==max(output2.10$Yield)),]$PS),
  #            data= subset(output2.10, HL == output2.10[which(output2.10$Yield==max(output2.10$Yield)),]$HL & SSS == output2.10[which(output2.10$Yield==max(output2.10$Yield)),]$SSS) )

# ggsave("Figure_S1.pdf", height = 15, width = 18, units = "cm") #PR
# ggsave("Figure_S2.pdf", height = 15, width = 18, units = "cm") #SR
# ggsave("Figure_S3.pdf", height = 15, width = 18, units = "cm") #Yield
