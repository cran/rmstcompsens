
#' @name diftimermst
#' @aliases diftimermst
#' @title Comparing restricted mean survival time when survival curves have divergent tails
#' @description Performs two-sample comparisons using the restricted mean survival time (RMST) when survival curves end at different time points between groups.
#' This package implements a sensitivity approach that allows the threshold  timepoint tau to be specified after the longest survival time in the shorter survival group.
#' Two kinds of between-group contrast estimators (the difference in RMST and the ratio of RMST) are computed.

#' @usage diftimermst(data, tau, alpha = 0.05)
#' @param data The data set must contains three variables named time, arm and event.\cr
#' time: The follow-up time for right censored.\cr
#' arm: The group indicator for comparison. 1=treat group, 0=control group.\cr
#' event: The event indicator, 1=event, 0=right censored.\cr
#' We assume that control group KM is shoter than treat group.
#' @param tau The truncation time point for the RMST calculation. If tau is bigger than the largest observed time in shorter KM, then conduct a sensitivity approach which we propose.
#' @param alpha The default is 0.05. 1-alpha confidence intervals are reported.
#' @return list of RMST comparing results.
#' @return 1. Point estimates
#' @return \item{RMST arm=1}{RMST results in treat group.}
#' @return \item{RMST arm=0 Method1}{RMST results in control group. Method1:Deal with the censored case as event case.}
#' @return \item{RMST arm=0 Method2}{RMST results in control group. Method2:Extend the observation time up to tau and deal with event.}
#' 2. Comparison between two group
#' @return \item{RMST difference Method1 (arm=1)-(arm=0)}{Results of the RMST difference in Method1.}
#' @return \item{RMST difference Method2 (arm=1)-(arm=0)}{Results of the RMST difference in Method2.}
#' @return \item{RMST ratio Method1 (arm=1)/(arm=0)}{Results of the RMST ratio in Method1.}
#' @return \item{RMST ratio Method2 (arm=1)/(arm=0)}{Results of the RMST ratio in Method2.}
#' The values below are generated when several patients are censored at the largest survival time in the shorter-survival group.
#' @return 1. Point estimates
#' @return \item{RMST arm=1}{RMST results in treat group.}
#' @return \item{RMST (difference p value=max)}{RMST results in control group. The pattern that maximize the p value when comparing by difference. Usually, the result of Method2 is calculated.}
#' @return \item{RMST (difference p value=min)}{RMST results in control group. The pattern that minimize the p value when comparing by difference. Usually, the result of Method1 is calculated.}
#' @return \item{RMST (ratio p value=max)}{RMST results in control group. The pattern that maximize the p value when comparing by ratio.}
#' @return \item{RMST (ratio p value=min)}{RMST results in control group. The pattern that minimize the p value when comparing by ratio.}
#' 2. Comparison between two group
#' @return \item{RMST difference (p value=max)}{Results of the RMST difference when the p value maximum.}
#' @return \item{RMST difference (p value=min)}{Results of the RMST difference when the p value minimum.}
#' @return \item{RMST ratio (p value=max)}{Results of the RMST ratio when the p value maximum.}
#' @return \item{RMST ratio (p value=max)}{Results of the RMST ratio when the p value minimum.}
#' @references Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D,
#' Takeuchi M, Uyama Y, Zhao L, Skali H, Solomon S, Jacobus S, Hughes M,
#' Packer M, Wei LJ. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis. Journal of clinical Oncology 2014, 32, 2380-2385.
#'
#' Hajime Uno, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui and James Bell (2020). survRM2: Comparing Restricted
#' Mean Survival Time. R package version 1.0-3. https://CRAN.R-project.org/package=survRM2
#'
#' Ueno K, Morita S. Sensitivity Analysis for Restricted Mean Survival Time When Survival Curves Have Divergent Tails. Ther Innov Regul Sci (2023).
#'
#' @author Kentaro Ueno
#' @examples
##' #--- sample data ---#
#' time  <- c(0.7,1.6,3.1,4.5,7.6,11,13.5,18.6,22.7,26.5,0.4,2.2,2.9,3.8,5.2,8.6,9.8,10.1,13.3,14.9)
#' event <- c(0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,1,0)
#' arm   <- c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
#' data <- data.frame(time,event,arm)
#' #--- analysis ---#
#' library(rmstcompsens)
#' a = diftimermst(data,24)
#' print(a)

#'
#' @import  survival dplyr
#' @importFrom stats pnorm qnorm time
#' @export

#########################################
#2022.1.10
#compare RMST at different time point
#assume control group KM is shoter than treat group
#short KM group: arm=0
#long KM group: arm=1
#########################################



diftimermst=function(data, tau, alpha=0.05){
  #-- data
  #-- tau (calculate RMST until tau)
  #-- alpha
  arm = NULL
  event = NULL
  flag = NULL
  Var1 = NULL
  Var2 = NULL
  Var3 = NULL
  #--- contents ---
  #1     case 1:number of last observation is one
  #1.1     calculate RMST (arm=1:longer survival group
  #1.2     calculate RMST (arm=0:shorter survival group and last observation is event)
  #1.3     calculate RMST (arm=0:shorter survival group and last observation is censor)
  #1.3.1     Convert censor to event
  #1.3.2     Convert survival time to tau
  #1.4     calculate difference and ratio
  #1.4.1     contrast (RMST difference)
  #1.4.2     contrast (RMST ratio)
  #1.5     organize for output
  #2     case 2:number of last observation are more than one
  #2.1.1     Convert censor to event
  #2.1.2     Convert survival time to tau
  #2.2     calculate RMST (arm=1:longer survival group)
  #2.3     calculate RMST (arm=0:shorter survival group and last observations are censor)
  #2.4     calculate difference and ratio
  #2.4.1     contrast (RMST difference)
  #2.4.2     contrast (RMST ratio)
  #2.5     organize for output
  #2.6     pick up max and min p value result
  #2.6.1     result of rmst difference
  #2.6.2     result of rmst ratio


  data=arrange(data,arm,time)
  data.short = subset(data,data$arm==0) #make arm=0(shorter KM) subset
  data.long = subset(data,data$arm==1) #make arm=1(longer KM) subset ##
  time.max = table(data.short$time==max(data.short$time)) #count number of max time in control group
  time.max.short = max(data.short$time) #max time in control group ##
  time.max.long = max(data.long$time) #max time in treat group ##
  flag.time = if(tau>time.max.long){flag.time=1}else if(tau<=time.max.short){flag.time=2}else{flag.time=3} ##

  #====================================================================
  #---0 case 0:tau is shorter than max time in control group
  #====================================================================
  if(flag.time==1){stop(paste("The truncation time, tau, needs to be shorter than or equal to max time"))} ##
  else if(flag.time==2){
    #==================================
    #  0.1 calculate RMST (arm=1:longer survival group)
    #==================================
    data.long = subset(data,arm==1)

    ft= survfit(Surv(time, event)~1,data.long)
    idx=ft$time<=tau

    wk.time=sort(c(ft$time[idx],tau))
    wk.surv=ft$surv[idx]
    wk.n.risk =ft$n.risk[idx]
    wk.n.event=ft$n.event[idx]

    time.diff <- diff(c(0, wk.time))
    areas <- time.diff * c(1, wk.surv)
    rmst = sum(areas)

    wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                     wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
    wk.var =c(wk.var,0)
    rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
    rmst.se  = sqrt(rmst.var)

    #--- output ---
    Zlong=list()
    Zlong$rmst = rmst
    Zlong$rmst.var = rmst.var

    #==================================
    #  0.2 calculate RMST (arm=0:shorter survival group)
    #==================================
    data.short = subset(data,arm==0)

    ft= survfit(Surv(time, event)~1,data.short)
    idx=ft$time<=tau

    wk.time=sort(c(ft$time[idx],tau))
    wk.surv=ft$surv[idx]
    wk.n.risk =ft$n.risk[idx]
    wk.n.event=ft$n.event[idx]

    time.diff <- diff(c(0, wk.time))
    areas <- time.diff * c(1, wk.surv)
    rmst = sum(areas)

    wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                     wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
    wk.var =c(wk.var,0)
    rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
    rmst.se  = sqrt(rmst.var)

    #--- output ---
    Zshort=list()
    Zshort$rmst = rmst
    Zshort$rmst.var = rmst.var

    #--- calculate difference and ratio ---
    rmst.diff     = Zlong$rmst-Zshort$rmst
    rmst.diff.se  = sqrt(Zlong$rmst.var + Zshort$rmst.var)
    rmst.diff.low = rmst.diff - qnorm(1-alpha/2)*rmst.diff.se
    rmst.diff.upp = rmst.diff + qnorm(1-alpha/2)*rmst.diff.se
    rmst.diff.pval   = pnorm(-abs(rmst.diff)/rmst.diff.se)*2
    rmst.diff.result = c(rmst.diff, rmst.diff.low, rmst.diff.upp, rmst.diff.pval)

    rmst.log.ratio     = log(Zlong$rmst) - log(Zshort$rmst)
    rmst.log.ratio.se  = sqrt(Zlong$rmst.var/Zlong$rmst/Zlong$rmst + Zshort$rmst.var/Zshort$rmst/Zshort$rmst)
    rmst.log.ratio.low = rmst.log.ratio - qnorm(1-alpha/2)*rmst.log.ratio.se
    rmst.log.ratio.upp = rmst.log.ratio + qnorm(1-alpha/2)*rmst.log.ratio.se
    rmst.log.ratio.pval   = pnorm(-abs(rmst.log.ratio)/rmst.log.ratio.se)*2
    rmst.ratio.result     = c(exp(rmst.log.ratio), exp(rmst.log.ratio.low), exp(rmst.log.ratio.upp),rmst.log.ratio.pval)

    #--- organize for output ---
    point.est = rbind(Zlong$rmst, Zshort$rmst)
    point.est.var = rbind(Zlong$rmst.var, Zshort$rmst.var)
    point.est.lci = rbind(Zlong$rmst-qnorm(1-alpha/2)*sqrt(Zlong$rmst.var), Zshort$rmst-qnorm(1-alpha/2)*sqrt(Zshort$rmst.var))
    point.est.uci = rbind(Zlong$rmst+qnorm(1-alpha/2)*sqrt(Zlong$rmst.var), Zshort$rmst+qnorm(1-alpha/2)*sqrt(Zshort$rmst.var))
    out.pointest = cbind(point.est, point.est.var, point.est.lci, point.est.uci)
    rownames(out.pointest) = c("RMST arm=1","RMST arm=0")
    colnames(out.pointest) = c("Est", "Var", "lower.95", "upper.95")

    out.compare=rbind(rmst.diff.result, rmst.ratio.result )
    rownames(out.compare)=c("RMST difference", "RMST ratio")
    colnames(out.compare)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")

    return(list(out.pointest,out.compare))

  }

  else{
    #====================================================================
    #---1 case 1:number of last observation is one
    #====================================================================
    if(time.max[2]==1){

      #==================================
      #  1.1 calculate RMST (arm=1:longer survival group)
      #==================================

      data.long = subset(data,arm==1)

      ft= survfit(Surv(time, event)~1,data.long)
      idx=ft$time<=tau

      wk.time=sort(c(ft$time[idx],tau))
      wk.surv=ft$surv[idx]
      wk.n.risk =ft$n.risk[idx]
      wk.n.event=ft$n.event[idx]

      time.diff <- diff(c(0, wk.time))
      areas <- time.diff * c(1, wk.surv)
      rmst = sum(areas)

      wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                       wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
      wk.var =c(wk.var,0)
      rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
      rmst.se  = sqrt(rmst.var)

      #--- output ---
      Z0=list()
      Z0$rmst = rmst
      Z0$rmst.var = rmst.var

      #==================================
      #  1.2 calculate RMST (arm=0:shorter survival group and last observation is event)
      #==================================

      nrow.maxtime = which.max(data.short$time)
      last.event=data.short$event[nrow.maxtime] #detect last=event or censor
      if(last.event==1){
        ft= survfit(Surv(time, event)~1,data.short)
        idx=ft$time<=tau

        wk.time=sort(c(ft$time[idx],tau))
        wk.surv=ft$surv[idx]
        wk.n.risk =ft$n.risk[idx]
        wk.n.event=ft$n.event[idx]

        time.diff <- diff(c(0, wk.time))
        areas <- time.diff * c(1, wk.surv)
        rmst = sum(areas)

        wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
        wk.var =c(wk.var,0)
        rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
        rmst.se  = sqrt(rmst.var)

        #--- output ---
        Z4=list()
        Z4$rmst = rmst
        Z4$rmst.var = rmst.var

        #--- calculate difference and ratio ---
        rmst.diff.04     = Z0$rmst-Z4$rmst
        rmst.diff.04.se  = sqrt(Z0$rmst.var + Z4$rmst.var)
        rmst.diff.04.low = rmst.diff.04 - qnorm(1-alpha/2)*rmst.diff.04.se
        rmst.diff.04.upp = rmst.diff.04 + qnorm(1-alpha/2)*rmst.diff.04.se
        rmst.diff.04.pval   = pnorm(-abs(rmst.diff.04)/rmst.diff.04.se)*2
        rmst.diff.04.result = c(rmst.diff.04, rmst.diff.04.low, rmst.diff.04.upp, rmst.diff.04.pval)

        rmst.log.ratio.04     = log(Z0$rmst) - log(Z4$rmst)
        rmst.log.ratio.04.se  = sqrt(Z0$rmst.var/Z0$rmst/Z0$rmst + Z4$rmst.var/Z4$rmst/Z4$rmst)
        rmst.log.ratio.04.low = rmst.log.ratio.04 - qnorm(1-alpha/2)*rmst.log.ratio.04.se
        rmst.log.ratio.04.upp = rmst.log.ratio.04 + qnorm(1-alpha/2)*rmst.log.ratio.04.se
        rmst.log.ratio.04.pval   = pnorm(-abs(rmst.log.ratio.04)/rmst.log.ratio.04.se)*2
        rmst.ratio.04.result     = c(exp(rmst.log.ratio.04), exp(rmst.log.ratio.04.low), exp(rmst.log.ratio.04.upp),rmst.log.ratio.04.pval)

        #--- organize for output ---
        point.est = rbind(Z0$rmst, Z4$rmst)
        point.est.var = rbind(Z0$rmst.var, Z4$rmst.var)
        out.pointest = cbind(point.est, point.est.var)
        rownames(out.pointest) = c("RMST arm=1","RMST arm=0")
        colnames(out.pointest) = c("Est", "Var")

        out.compare=rbind(rmst.diff.04.result, rmst.ratio.04.result )
        rownames(out.compare)=c("RMST difference", "RMST ratio")
        colnames(out.compare)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")

      }else{

        #==================================
        #  1.3 calculate RMST (arm=0:shorter survival group and last observation is censor)
        #==================================

        #--- 1.3.1 Convert censor to event ---

        data.long = subset(data,arm==1)
        data.short = subset(data,arm==0)
        event.short = data.short$event
        event.short.new = ifelse(data.short$time== max(data.short$time),1, data.short$event)
        data.short.new = data.frame(data.short,event.short.new)
        data.short.new = select(data.short.new,-event)
        data.short.new = rename(data.short.new,event=event.short.new)
        data.1.new = rbind(data.long,data.short.new)
        data.1.new.short = subset(data.1.new,arm==0)

        #--- calculate RMST ---
        ft= survfit(Surv(time, event)~1,data.1.new.short)
        idx=ft$time<=tau

        wk.time=sort(c(ft$time[idx],tau))
        wk.surv=ft$surv[idx]
        wk.n.risk =ft$n.risk[idx]
        wk.n.event=ft$n.event[idx]

        time.diff <- diff(c(0, wk.time))
        areas <- time.diff * c(1, wk.surv)
        rmst = sum(areas)

        wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
        wk.var =c(wk.var,0)
        rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
        rmst.se  = sqrt(rmst.var)

        #--- output ---
        Z1=list()
        Z1$rmst = rmst
        Z1$rmst.var = rmst.var


        #--- 1.3.2 Convert survival time to tau and event---

        data.long = subset(data,arm==1)
        data.short = subset(data,arm==0)
        event.short = data.short$event
        event.short.new = ifelse(data.short$time== max(data.short$time),1, data.short$event)
        time.short = data.short$time
        time.short.new = ifelse(time.short==max(time.short),tau,time.short)
        data.short.new = data.frame(data.short,event.short.new,time.short.new)
        data.short.new = select(data.short.new,-event)
        data.short.new = rename(data.short.new,event=event.short.new)
        data.short.new = select(data.short.new,-time)
        data.short.new = rename(data.short.new,time=time.short.new)
        data.2.new = rbind(data.long,data.short.new)
        data.2.new.short = subset(data.2.new,arm==0)

        #--- calculate RMST ---
        ft= survfit(Surv(time, event)~1,data.2.new.short)
        idx=ft$time<=tau

        wk.time=sort(c(ft$time[idx],tau))
        wk.surv=ft$surv[idx]
        wk.n.risk =ft$n.risk[idx]
        wk.n.event=ft$n.event[idx]

        time.diff <- diff(c(0, wk.time))
        areas <- time.diff * c(1, wk.surv)
        rmst = sum(areas)

        wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
        wk.var =c(wk.var,0)
        rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
        rmst.se  = sqrt(rmst.var)

        #--- output ---
        Z2=list()
        Z2$rmst = rmst
        Z2$rmst.var = rmst.var


        #--- 1.4 calculate difference and ratio ---
        #--- 1.4.1 contrast (RMST difference) ---
        rmst.diff.01     = Z0$rmst-Z1$rmst
        rmst.diff.01.se  = sqrt(Z0$rmst.var + Z1$rmst.var)
        rmst.diff.01.low = rmst.diff.01 - qnorm(1-alpha/2)*rmst.diff.01.se
        rmst.diff.01.upp = rmst.diff.01 + qnorm(1-alpha/2)*rmst.diff.01.se
        rmst.diff.01.pval   = pnorm(-abs(rmst.diff.01)/rmst.diff.01.se)*2
        rmst.diff.01.result = c(rmst.diff.01, rmst.diff.01.low, rmst.diff.01.upp, rmst.diff.01.pval)

        rmst.diff.02     = Z0$rmst-Z2$rmst
        rmst.diff.02.se  = sqrt(Z0$rmst.var + Z2$rmst.var)
        rmst.diff.02.low = rmst.diff.02 - qnorm(1-alpha/2)*rmst.diff.02.se
        rmst.diff.02.upp = rmst.diff.02 + qnorm(1-alpha/2)*rmst.diff.02.se
        rmst.diff.02.pval   = pnorm(-abs(rmst.diff.02)/rmst.diff.02.se)*2
        rmst.diff.02.result = c(rmst.diff.02, rmst.diff.02.low, rmst.diff.02.upp, rmst.diff.02.pval)


        #--- 1.4.2 contrast (RMST ratio) ---
        rmst.log.ratio.01     = log(Z0$rmst) - log(Z1$rmst)
        rmst.log.ratio.01.se  = sqrt(Z0$rmst.var/Z0$rmst/Z0$rmst + Z1$rmst.var/Z1$rmst/Z1$rmst)
        rmst.log.ratio.01.low = rmst.log.ratio.01 - qnorm(1-alpha/2)*rmst.log.ratio.01.se
        rmst.log.ratio.01.upp = rmst.log.ratio.01 + qnorm(1-alpha/2)*rmst.log.ratio.01.se
        rmst.log.ratio.01.pval   = pnorm(-abs(rmst.log.ratio.01)/rmst.log.ratio.01.se)*2
        rmst.ratio.01.result     = c(exp(rmst.log.ratio.01), exp(rmst.log.ratio.01.low), exp(rmst.log.ratio.01.upp),rmst.log.ratio.01.pval)

        rmst.log.ratio.02     = log(Z0$rmst) - log(Z2$rmst)
        rmst.log.ratio.02.se  = sqrt(Z0$rmst.var/Z0$rmst/Z0$rmst + Z2$rmst.var/Z2$rmst/Z2$rmst)
        rmst.log.ratio.02.low = rmst.log.ratio.02 - qnorm(1-alpha/2)*rmst.log.ratio.02.se
        rmst.log.ratio.02.upp = rmst.log.ratio.02 + qnorm(1-alpha/2)*rmst.log.ratio.02.se
        rmst.log.ratio.02.pval   = pnorm(-abs(rmst.log.ratio.02)/rmst.log.ratio.02.se)*2
        rmst.ratio.02.result     = c(exp(rmst.log.ratio.02), exp(rmst.log.ratio.02.low), exp(rmst.log.ratio.02.upp),rmst.log.ratio.02.pval)

        #--- 1.5 organize for output ---
        point.est = rbind(Z0$rmst, Z1$rmst, Z2$rmst)
        point.est.se = rbind(sqrt(Z0$rmst.var), sqrt(Z1$rmst.var), sqrt(Z2$rmst.var))
        rmst.Z0.CI = c(Z0$rmst-qnorm(1-alpha/2)*sqrt(Z0$rmst.var),Z0$rmst+qnorm(1-alpha/2)*sqrt(Z0$rmst.var))
        rmst.Z1.CI = c(Z1$rmst-qnorm(1-alpha/2)*sqrt(Z1$rmst.var),Z1$rmst+qnorm(1-alpha/2)*sqrt(Z1$rmst.var))
        rmst.Z2.CI = c(Z2$rmst-qnorm(1-alpha/2)*sqrt(Z2$rmst.var),Z2$rmst+qnorm(1-alpha/2)*sqrt(Z2$rmst.var))
        rmst.CI = rbind(rmst.Z0.CI,rmst.Z1.CI,rmst.Z2.CI)
        out.pointest = cbind(point.est, point.est.se,rmst.CI)
        rownames(out.pointest) = c("RMST arm=1","RMST arm=0 Method1","RMST arm=0 Method2")
        colnames(out.pointest) = c("Est", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))

        out.compare=rbind(rmst.diff.01.result, rmst.diff.02.result, rmst.ratio.01.result, rmst.ratio.02.result )
        rownames(out.compare)=c("RMST difference Method1 (arm=1)-(arm=0)","RMST difference Method2 (arm=1)-(arm=0)","RMST ratio Method1 (arm=1)/(arm=0)","RMST ratio Method2 (arm=1)/(arm=0)")
        colnames(out.compare)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")
      }

      return(list(out.pointest,out.compare))
    }else{

      #====================================================================
      #---2 case 2:number of last observation are more than one
      #====================================================================

      #--- make empty matrix for output ---
      list.pointest.pat3=matrix(,nrow=0,ncol=4)
      list.rmstdiff=matrix(,nrow=0,ncol=4)
      list.rmstratio=matrix(,nrow=0,ncol=4)

      #--- calculate every pattern of combination(x:convert event, y:extent survival time to tau, z: stay) ---
      x=0:time.max[2]
      y=0:time.max[2]
      z=0:time.max[2]
      ptrn=expand.grid(x,y,z) #make every pattern of x,y,z
      combination=subset(ptrn,(Var1+Var2+Var3)==time.max[2]) #only subset of change variable
      combination$flag=ifelse(combination$Var2==0&combination$Var3>0,1,0) #event comes earlier than censor, so exclude event and censor pattern on same time
      combi = subset(combination,flag==0)

      for (i in 1:nrow(combi)) {
        a = combi[i,1]
        b = combi[i,2]

        #--- 2.1.1 Convert censor to event ---
        if(a==0){
          data.1.new = data
        }else{
          data.long = subset(data,arm==1)
          data.short = subset(data,arm==0)
          event.short = data.short$event
          event.short.new = ifelse(data.short$time== max(data.short$time)& row(data.short)[,1]==nrow(data.short),1, data.short$event)
          #time=max and row=max-1 then convert to event
          data.short.new = data.frame(data.short,event.short.new)
          data.short.new = select(data.short.new,-event)
          data.short.new = rename(data.short.new,event=event.short.new)
          for(i in 0:a){
            event.short.new = ifelse(data.short.new$time== max(data.short.new$time)& row(data.short.new)[,1]==nrow(data.short.new)-(i-1),1, data.short.new$event)
            data.short.new = data.frame(data.short,event.short.new)
            data.short.new = select(data.short.new,-event)
            data.short.new = rename(data.short.new,event=event.short.new)}
          data.1.new = rbind(data.long,data.short.new)}

        #--- 2.1.2 Convert survival time to tau and event---
        if(b==0){
          data.2.new = data.1.new
        }else{
          data.long = subset(data.1.new,arm==1)
          data.short = subset(data.1.new,arm==0)
          time.short = data.short$time
          num.maxtime=table(data.short$time==max(data.short$time))
          time.short.new = ifelse(time.short==max(time.short)&row(data.short)[,1]==nrow(data.short)-(num.maxtime[2]-1),tau,time.short)
          data.short.new = data.frame(data.short,time.short.new)
          data.short.new = select(data.short.new,-time)
          data.short.new = rename(data.short.new,time=time.short.new)
          data.2.new = rbind(data.long,data.short.new)
          if(b>1){
            for(i in 0:b){
              time.short.new = ifelse(time.short.new==max(time.short)&row(data.short.new)[,1]==nrow(data.short.new)-(num.maxtime[2]-i),tau,time.short.new)
              data.short.new = data.frame(data.short,time.short.new)
              data.short.new = select(data.short.new,-time)
              data.short.new = rename(data.short.new,time=time.short.new)
              data.2.new = rbind(data.long,data.short.new)}}
          data.2.new = data.2.new[order(data.2.new$arm,data.2.new$time),]
        }

        if(b!=0){
          data.long = subset(data.2.new,arm==1)
          data.short = subset(data.2.new,arm==0)
          event.short = data.short$event
          event.short.new = ifelse(data.short$time== max(data.short$time)& row(data.short)[,1]==nrow(data.short),1, data.short$event)
          #time=max and row=max-1 then convert to event
          data.short.new = data.frame(data.short,event.short.new)
          data.short.new = select(data.short.new,-event)
          data.short.new = rename(data.short.new,event=event.short.new)
          for(i in 0:b){
            event.short.new = ifelse(data.short.new$time== max(data.short.new$time)& row(data.short.new)[,1]==nrow(data.short.new)-(i-1),1, data.short.new$event)
            data.short.new = data.frame(data.short,event.short.new)
            data.short.new = select(data.short.new,-event)
            data.short.new = rename(data.short.new,event=event.short.new)}
          data.2.new = rbind(data.long,data.short.new)}


        #==================================
        #  2.2 calculate RMST (arm=1:longer survival group)
        #==================================

        data.2.new.long = subset(data.2.new,arm==1)

        ft= survfit(Surv(time, event)~1,data.2.new.long)
        idx=ft$time<=tau

        wk.time=sort(c(ft$time[idx],tau))
        wk.surv=ft$surv[idx]
        wk.n.risk =ft$n.risk[idx]
        wk.n.event=ft$n.event[idx]

        time.diff <- diff(c(0, wk.time))
        areas <- time.diff * c(1, wk.surv)
        rmst = sum(areas)

        wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
        wk.var =c(wk.var,0)
        rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
        rmst.se  = sqrt(rmst.var)

        #--- output ---
        Z0=list()
        Z0$rmst = rmst
        Z0$rmst.var = rmst.var


        #==================================
        #  2.3 calculate RMST (arm=0:shorter survival group and last observations are censor)
        #==================================

        data.2.new.short = subset(data.2.new,arm==0)
        ft= survfit(Surv(time, event)~1,data.2.new.short)
        idx=ft$time<=tau

        wk.time=sort(c(ft$time[idx],tau))
        wk.surv=ft$surv[idx]
        wk.n.risk =ft$n.risk[idx]
        wk.n.event=ft$n.event[idx]

        time.diff <- diff(c(0, wk.time))
        areas <- time.diff * c(1, wk.surv)
        rmst = sum(areas)

        wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
        wk.var =c(wk.var,0)
        rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
        rmst.se  = sqrt(rmst.var)

        #--- output ---
        Z3=list()
        Z3$rmst = rmst
        Z3$rmst.var = rmst.var


        #--- 2.4 calculate difference and ratio ---
        #--- 2.4.1 contrast (RMST difference) ---
        rmst.diff.03     = Z0$rmst-Z3$rmst
        rmst.diff.03.se  = sqrt(Z0$rmst.var + Z3$rmst.var)
        rmst.diff.03.low = rmst.diff.03 - qnorm(1-alpha/2)*rmst.diff.03.se
        rmst.diff.03.up = rmst.diff.03 + qnorm(1-alpha/2)*rmst.diff.03.se
        rmst.diff.03.pval   = pnorm(-abs(rmst.diff.03)/rmst.diff.03.se)*2
        rmst.diff.03.result = c(rmst.diff.03, rmst.diff.03.low, rmst.diff.03.up, rmst.diff.03.pval)

        #--- 2.4.2 contrast (RMST ratio) ---
        rmst.log.ratio.03     = log(Z0$rmst) - log(Z3$rmst)
        rmst.log.ratio.03.se  = sqrt(Z0$rmst.var/Z0$rmst/Z0$rmst + Z3$rmst.var/Z3$rmst/Z3$rmst)
        rmst.log.ratio.03.low = rmst.log.ratio.03 - qnorm(1-alpha/2)*rmst.log.ratio.03.se
        rmst.log.ratio.03.up = rmst.log.ratio.03 + qnorm(1-alpha/2)*rmst.log.ratio.03.se
        rmst.log.ratio.03.pval   = pnorm(-abs(rmst.log.ratio.03)/rmst.log.ratio.03.se)*2
        rmst.ratio.03.result     = c(exp(rmst.log.ratio.03), exp(rmst.log.ratio.03.low), exp(rmst.log.ratio.03.up),rmst.log.ratio.03.pval)

        #--- 2.5 organize for output ---
        point.est = rbind(Z0$rmst, Z3$rmst)
        point.est.se = rbind(sqrt(Z0$rmst.var),  sqrt(Z3$rmst.var))
        rmst.Z0.CI = c(Z0$rmst-qnorm(1-alpha/2)*sqrt(Z0$rmst.var),Z0$rmst+qnorm(1-alpha/2)*sqrt(Z0$rmst.var))
        rmst.Z3.CI = c(Z3$rmst-qnorm(1-alpha/2)*sqrt(Z3$rmst.var),Z3$rmst+qnorm(1-alpha/2)*sqrt(Z3$rmst.var))
        rmst.CI = rbind(rmst.Z0.CI,rmst.Z3.CI)
        out.pointest = cbind(point.est, point.est.se, rmst.CI)
        rownames(out.pointest) = c("RMST arm=1","RMST arm=0 pattern3")
        colnames(out.pointest) = c("Est", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))

        out.compare=rbind(rmst.diff.03.result, rmst.ratio.03.result )
        rownames(out.compare)=c("RMST difference Method3 (arm=1)-(arm=0)", "RMST ratio Method3 (arm=1)/(arm=0)")
        colnames(out.compare)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")

        list.pointest.pat0 = out.pointest[1,]
        list.pointest.pat3 = rbind(list.pointest.pat3,out.pointest[2,])
        list.rmstdiff=rbind(list.rmstdiff,out.compare[1,])
        list.rmstratio=rbind(list.rmstratio,out.compare[2,])
      }

      #--- 2.6 pick up max and min p value result ---
      #--- list of all result ---
      result.pointest.treat = list.pointest.pat0
      result.list.pointest = data.frame(list.pointest.pat3)
      result.list.rmstdiff = data.frame(list.rmstdiff)
      result.list.rmstratio = data.frame(list.rmstratio)

      #--- 2.6.1 result of rmst difference ---
      result.diff.max = subset(result.list.rmstdiff,result.list.rmstdiff$p==max(result.list.rmstdiff$p))[1,]
      result.diff.min = subset(result.list.rmstdiff,result.list.rmstdiff$p==min(result.list.rmstdiff$p))[1,]
      result.diff=rbind(result.diff.max,result.diff.min)
      rownames(result.diff)=c("RMST difference (p value=max)","RMST difference (p value=min)")

      result.pointest.diff.max = result.list.pointest[which.max(result.list.rmstdiff$p),]
      result.pointest.diff.min = result.list.pointest[which.min(result.list.rmstdiff$p),]
      result.pointest.diff = rbind(result.pointest.treat,result.pointest.diff.max,result.pointest.diff.min)
      rownames(result.pointest.diff)=c("RMST arm=1","RMST (difference p value=max)","RMST (difference p value=min)")

      #--- 2.6.2 result of rmst ratio ---
      result.ratio.max = subset(result.list.rmstratio,result.list.rmstratio$p==max(result.list.rmstratio$p))[1,]
      result.ratio.min = subset(result.list.rmstratio,result.list.rmstratio$p==min(result.list.rmstratio$p))[1,]
      result.ratio=rbind(result.ratio.max,result.ratio.min)
      rownames(result.ratio)=c("RMST ratio (p value=max)","RMST ratio (p value=min)")

      result.pointest.ratio.max = result.list.pointest[which.max(result.list.rmstratio$p),]
      result.pointest.ratio.min = result.list.pointest[which.min(result.list.rmstratio$p),]
      result.pointest.ratio = rbind(result.pointest.ratio.max,result.pointest.ratio.min)
      rownames(result.pointest.ratio)=c("RMST (ratio p value=max)","RMST (ratio p value=min)")

      #--- output ---
      out.pointest = rbind(result.pointest.diff,result.pointest.ratio)
      out.compare = rbind(result.diff,result.ratio)

      return(list(out.pointest,out.compare))
    }
  }
}
