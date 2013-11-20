#******************************************************************************
# Model for 2013-14 AMS solar energy prediction
# Author: Owen Zhang Date: 2013-11-20
# Please see LICENSE file for licensing information
# Please see README (in raw format)
#******************************************************************************
setwd(".") #Please see readme regarding files required for input

require(reshape2)
require(sqldf)
require(gbm)
require(ncdf4)
require(Matrix)
require(data.table)
require(foreach)
require(doMC)
registerDoMC()

#simple utlity function
my_cap<-function(v, l, h){
  v[v<l]<-l
  v[v>h]<-h
  return(v)
}

#global paramter
reload_data<-T  #either to reload all data from scratch
train_tv<-T     #either traing on all data or save a portion (for my own testing)
n_trees<-3000   #how many trees


if(reload_data) {
  #function to read all the files in the base_folder into 5 dim arrays
  loadxx<-function(base_folder, t3) {
    t_fns<-dir(base_folder)
    print(t_fns)
    t_fns_short<-substr(t_fns,1,7)
    x0list<-list()
    for(j in 1:length(t_fns)) {
      nc<-nc_open(paste(base_folder, "/", t_fns[j], sep=""))
      x0list[[j]]<-ncvar_get(nc,  names(nc$var)[3])
      print(paste("file ", t_fns[j], ", var=", names(nc$var)[3], "dim=", paste(dim(x0list[[j]]), collapse=",")))
    }  
    
    xxlist<-foreach(pidx=1:4) %dopar% {
      xx<-data.table(matrix(0, nrow=dim(t3)[1], ncol=5*length(t_fns)))
      setnames(xx, paste("value", c(1:dim(xx)[2]), sep="_"))
      for(j in 1:length(t_fns)) {
        x0<-x0list[[j]]
        print(paste("file ", t_fns[j], ", var=", names(nc$var)[3], "dim=", paste(dim(x0), collapse=",")))
        col_range<-c(((j-1)*5+1):(j*5))
        setnames(xx, col_range, paste(t_fns_short[j], c(1:5), sep="_"))
        
        for(i in 1:98){
          if(pidx>2) lon_idx<-which(nc$dim$lon$vals==floor(stations$elon[i]+360)) else lon_idx<-which(nc$dim$lon$vals==ceiling(stations$elon[i]+360))
          if(pidx %% 2==0) lat_idx<-which(nc$dim$lat$vals==floor(stations$nlat[i])) else lat_idx<-which(nc$dim$lat$vals==ceiling(stations$nlat[i]))
          
          system.time(print(paste("station [", stations$stid[i], "], with lon_idx=", lon_idx,",lat_idx=", lat_idx)))
          system.time(x1<-x0[lon_idx,lat_idx,,,])
          #system.time(x2<-t(apply(x1, c(1,3), mean)))
          system.time(x2<-t(colMeans(aperm(x1,c(2,1,3))))) 
          system.time(for(k in 1:5) set(xx, which(as.character(t3$variable)==as.character(stations$stid[i])), as.integer((j-1)*5+k), x2[,k]))
        }
        
      }
      return(xx)
    }
    return(xxlist)
  }
  
  #function to average the 11 member and find the nearest 4 corner to convert 5 dim array to 15x5x4 features
  get_xxs<-function(xl, tt) {
    r<-list()
    xx1<-copy(xl[[1]])
    allones<-c(1:dim(xx1)[1])
    for(i in 1:dim(xx1)[2]) {
      print(i)
      x<-0
      for(j in 1:4) {
        x<-x+xl[[j]][,i,with=F]*tt[, paste("dw", j, sep="")]
      }
      set(xx1,allones,i,x)
    }
    r$xx<-xx1
    
    xx1_nearest<-copy(xl[[1]])
    allones<-c(1:dim(xx1_nearest)[1])
    for(i in 1:dim(xx1_nearest)[2]) {
      print(i)
      x<-rep(0, dim(xx1_nearest)[1])
      for(j in 1:4) {
        tmpx<-unlist(xl[[j]][tt$nearest==j,i,with=F])
        x[which(tt$nearest==j)]<-tmpx
      }
      set(xx1_nearest,allones,i,x)
    }
    r$xx_nearest<-xx1_nearest
    return(r)
  }
  
  
  #read station data and find the corners
  stations<-read.csv(file="station_info.csv")
  stations$intlon<-round(stations$elon+360)
  stations$intlat<-round(stations$nlat)
  stations$lon_diff<-with(stations, elon-intlon)
  stations$lat_diff<-with(stations, nlat-intlat)
  
  for(pidx in 1:4) {
    stations[, paste("d", pidx, sep="")]<-(-1)
    for(i in 1:dim(stations)[1]) {
      #print(i)
      if(pidx>2) lon2<-floor(stations$elon[i]+360) else lon2<-ceiling(stations$elon[i]+360)
      if(pidx %% 2==0) lat2<-floor(stations$nlat[i]) else lat2<-ceiling(stations$nlat[i])
      stations[i, paste("d", pidx, sep="")]<-sqrt(((stations$elon[i]+360-lon2)*.81)^2+(stations$nlat[i]-lat2)^2)
    }
  }
  
  #calculate the distance from the 4 corners in a very silly way
  d0<-0.01
  dsum<-with(stations, 1/(d1+d0)+1/(d2+d0)+1/(d3+d0)+1/(d4+d0))
  stations$dw1<-with(stations, 1/(d0+d1)/dsum)
  stations$dw2<-with(stations, 1/(d0+d2)/dsum)
  stations$dw3<-with(stations, 1/(d0+d3)/dsum)
  stations$dw4<-with(stations, 1/(d0+d4)/dsum)
  
  #calculate the nearest corner in a very silly way
  stations$dmin<-stations$d1
  stations$nearest<-1
  f_tmp<-stations$d2<stations$dmin
  stations$nearest[f_tmp]<-2
  stations$dmin[f_tmp]<-stations$d2[f_tmp]
  f_tmp<-stations$d3<stations$dmin
  stations$nearest[f_tmp]<-3
  stations$dmin[f_tmp]<-stations$d3[f_tmp]
  f_tmp<-stations$d4<stations$dmin
  stations$nearest[f_tmp]<-4
  stations$dmin[f_tmp]<-stations$d4[f_tmp]
  
  
  t1<-read.csv(file="train.csv")
  t2<-melt(t1, id.vars=c("Date"))
  t3<-sqldf("select a.*, b.intlon, b.intlat, b.elev, b.lon_diff, b.lat_diff, dw1, dw2, dw3, dw4, nearest from t2 a 
    left join stations b on a.variable=b.stid")
  t3<-t3[order(t3$variable, t3$Date),]
  #load ncdf4 files into 5 dim arrays
  xxlist<-loadxx("train", t3)
  #average 11 forecasts, then find the nearest 4 corners, so the data become 15x5x4 features for each station*day
  xxs<-get_xxs(xxlist, t3)
  
  
  sub0<-read.csv(file="sampleSubmission.csv")
  t2h<-melt(sub0, id.vars=c("Date"))
  t3h<-sqldf("select a.*, b.intlon, b.intlat, b.elev, b.lon_diff, b.lat_diff, dw1, dw2, dw3, dw4, nearest from t2h a 
    left join stations b on a.variable=b.stid")
  t3h<-t3h[order(t3h$variable, t3h$Date),]
  #load ncdf4 files into 5 dim arrays
  xxlisth<-loadxx("test", t3h)
  #average 11 forecasts, then find the nearest 4 corners, so the data become 15x5x4 features for each station*day
  xxsh<-get_xxs(xxlisth, t3h)
  
  #also add the 4 corners of dwsrf_4 raw values into the data frame
  xx_dswrfs<-cbind(xxlist[[1]][, 14, with=F], xxlist[[2]][, 14, with=F], xxlist[[3]][, 14, with=F] ,xxlist[[4]][, 14, with=F])
  setnames(xx_dswrfs, 1:4, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"))
  xxh_dswrfs<-cbind(xxlisth[[1]][, 14, with=F], xxlisth[[2]][, 14, with=F], xxlisth[[3]][, 14, with=F] ,xxlisth[[4]][, 14, with=F])
  setnames(xxh_dswrfs, 1:4, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"))
  
  
  #prepare 2 data frames, t4list[[1]] is the weighted average of 4 corners, t4list[[2]] is the nearest corner
  t4list<-list()
  for(i in 1:2) {
    if(i==1) {
      t4h<-cbind(t3h, xxsh$xx, xxh_dswrfs[, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"), with=F])
      t4<-cbind(t3, xxs$xx, xx_dswrfs[, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"), with=F])    
    } else {
      t4h<-cbind(t3h, xxsh$xx_nearest, xxh_dswrfs[, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"), with=F])
      t4<-cbind(t3, xxs$xx_nearest, xx_dswrfs[, c("dswrf_4_1", "dswrf_4_2", "dswrf_4_3", "dswrf_4_4"), with=F])    
    }
    
    #hard code the day idx
    t4h$day_idx<-c(5114:(5113+1796))
    y_idx<-t4h$day_idx/(5113/14)
    y_idx<-y_idx-floor(y_idx)
    t4h$doy_idx<-y_idx
    t4h$year<-as.integer(t4h$Date/10000L)
    t4h$split1<-2
    
  
    #hard code the day idx
    t4$day_idx<-c(1:5113)
    y_idx<-t4$day_idx/(5113/14)
    y_idx<-y_idx-floor(y_idx)
    t4$doy_idx<-y_idx
    t4$year<-as.integer(t4$Date/10000L)
    t4$split1<-0
    t4$split1[t4$year>2004]<-1
    
    #combind training and test data
    t4<-rbind(t4, t4h)
    t4<-t4[order(t4$split1),]
    t4$w<-1
    t4$w[t4$split1>1]<-0
    
    #linear combination of dswrf @ 5 time points
    t4$dswrf<-with(t4, dswrf_s_1*-.5+dswrf_s_2*-.1+dswrf_s_3+dswrf_s_4+dswrf_s_5*.8)
    
    t4list[[i]]<-t4
  }
  
  #save the data frames in case something goes wrong and we need to restart
  save(t4list, file="solar_final_t4list.RData")
} else {
  load(file="solar_final_t4list.RData")
}

#use all variables except for the ones listed below in the model
vns<-names(t4list[[1]])
vns<-vns[which(vns %in% c("Date", "w", "variable", "date2", "pred3", "pred4", "pred2", "pred5", "d2_station", "value", "day_idx", "split1", "pred", "dw1","dw2","dw3", "dw4")==F)]
all_vars<-paste(vns, collapse="+")

#fit 2 gbms glist[[1]] is based on the weighted average of 4 corners glist[[2]] is based on the nearest corner
glist<-list()
#plist has the predictions for testing purpose
plist<-list()
i_range<-c(1:2)
for(i in i_range){
  t4tmp<-t4list[[i]]
  if(train_tv) {
    #traing using all training data for submission
    t4t<-t4tmp[t4tmp$split1<=1, ]    
    t_frac<-1
  } else {
    #trainig use partial data for model tuning
    t4t<-t4tmp[t4tmp$split1<1, ]
    t_frac<-with(t4tmp, sum(split1==0)/length(split1))
  }
  #fit GBM, this is going to take a while!
  glist[[i]]<-gbm(as.formula(paste("value~", all_vars, "")), data=t4t
                  ,train.fraction=t_frac, shrinkage=0.05, n.trees=n_trees, n.minobsinnode=10
                  ,distribution="laplace", interaction.depth=10)
  plist[[i]]<-predict.gbm(glist[[i]], newdata=t4tmp[t4tmp$split1==1,], n.trees=glist[[i]]$n.trees)
  print(paste("*************************** ", glist[[i]]$n.trees, mean(abs(t4tmp$value[t4tmp$split1==1]-plist[[i]]))))
}
#testing the values to make sure nothing went terribly wrong
predv<-0
wk<-c(1,0.5)
for(i in i_range){
  predv<-predv+plist[[i]]*wk[i]/sum(wk)
}
print(paste("*************************** ", glist[[i_range[1]]]$n.trees, mean(abs(t4tmp$value[t4tmp$split1==1]-my_cap(predv, 1.3e6, 5e7)))))

#score test data
predh<-0
for(i in i_range){
  t4tmp<-t4list[[i]]
  t4h<-t4tmp[t4tmp$split1==2,]
  pred_tmp<-predict.gbm(glist[[i]], t4h, n.trees=glist[[i]]$n.trees)
  predh<-predh+pred_tmp*wk[i]/sum(wk)
}

t4h$pred<-my_cap(predh, 1.3e6, 5e7)

#output reslts
tmph<-t4h[, c("Date","variable","pred")]
sub1<-dcast(tmph, Date~variable, value.var="pred")
write.csv(sub1, file="sub_final.csv", row.names=F)

