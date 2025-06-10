load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

dq       = str_get_dq()
date_cmd = "date +" + dq + "%a %b %d %H:%M:%S %Z %Y" + dq
TotalTime     = systemfunc(date_cmd)


;;==========================
institute=(/"AS-RCEC","AWI","BCC","CAMS","CAS","CAS","CAS","CCCma","CCCma","CMCC","CMCC","CNRM-CERFACS","CNRM-CERFACS","CSIRO","CSIRO-ARCCSS",\
            "EC-Earth-Consortium","EC-Earth-Consortium","EC-Earth-Consortium","FIO-QLNM","INM","INM","IPSL","MIROC","MIROC","MOHC","MOHC",\
            "MPI-M","MRI","NASA-GISS","NCAR","NCC","NCC","NOAA-GFDL","NOAA-GFDL","THU","UA"/)


model = (/"TaiESM1","AWI-CM-1-1-MR","BCC-CSM2-MR","CAMS-CSM1-0", "CAS-ESM2-0","FGOALS-f3-L","FGOALS-g3","CanESM5","CanESM5-CanOE","CMCC-CM2-SR5", "CMCC-ESM2","CNRM-CM6-1","CNRM-ESM2-1","ACCESS-ESM1-5", "ACCESS-CM2",\
         "EC-Earth3","EC-Earth3-CC", "EC-Earth3-Veg","FIO-ESM-2-0","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","MIROC6","MIROC-ES2L","HadGEM3-GC31-LL","UKESM1-0-LL",\
          "MPI-ESM1-2-LR","MRI-ESM2-0","GISS-E2-1-G","CESM2","NorESM2-MM", "NorESM2-LM","GFDL-CM4","GFDL-ESM4","CIESM","MCM-UA-1-0"/)

nmodels = dimsizes(model)
;;--------------------------------------------------------------

;; for natural variability (N)--------
nyrs_pi       = 500
missing      = -999999999999.
NAJ_pi = new((/nmodels,nyrs_pi/), float,missing)
dir_pi   = "/home/users/shengc01/work11/scripts9/CMIP6_data_ts_pi_hist/v2/piControl/"
do im = 0, nmodels - 1 ,1 
    f_pi   = addfile(dir_pi+"NAJ.lat.JJA."+model(im)+".45W-0.based-maxU.piControl.nc", "r")
    NAJ_pi(im,:) = f_pi->NAJ_lat(0:nyrs_pi-1)
end do 
    ;NAJ_pi = dim_rmvmean_n(NAJ_pi, 0);;the departure from the ensemble mean for each member
    NAJ_pi = runave_n_Wrap(NAJ_pi,20,0,1)  
    NAJ_noise = dim_stddev_n_Wrap(NAJ_pi,1)
    NAJ_noise_MME = dim_avg(NAJ_noise)

print(NAJ_noise_MME)
printVarSummary(NAJ_noise)

print(NAJ_noise)
;;---------------------------------


;;For trend signal (S)-------------
YYYY  = ispan(1850,2014,1)

YEAR_strt = 1850
Yrstrt1 = ind(YYYY.eq.YEAR_strt)
Yrend1  = ind(YYYY.eq.2014)

YYYY  := ispan(2015,2100,1)
Yrstrt2 = ind(YYYY.eq.2015)
Yrend2  = ind(YYYY.eq.2099)

dir_hist   = "/home/users/shengc01/work11/scripts9/CMIP6_data_ts_pi_hist/v2/historical/"
dir_ssp245   = "/home/users/shengc01/work11/scripts9/CMIP6_data_ts_pi_hist/v2/ssp245/"
dir_ssp585   = "/home/users/shengc01/work11/scripts9/CMIP6_data_ts_pi_hist/v2/ssp585/"
nyrs          = 2099-YEAR_strt+1
;;--
dirout = "/home/users/shengc01/work11/scripts3/data/"
filout = "NAJ.lat.JJA.NCEP.v2.nc"
fout = addfile(dirout + filout, "r")
NAJ_obs = fout->NCEP_jet_lat1(YEAR_strt-1806:2015-1806)
NAJ_obs = runave_n_Wrap(NAJ_obs,20,0,0) - avg(NAJ_obs)
NN = dimsizes(NAJ_obs)

print("llll")
;;
NAJ_hist_ssp245 = new((/nmodels,nyrs/), float,missing)
NAJ_hist_ssp585 = new((/nmodels,nyrs/), float,missing)

do im = 0, nmodels - 1 ,1 
    print(model(im)+"..")
    f_hist   = addfile(dir_hist+"NAJ.lat.JJA."+model(im)+".45W-0.based-maxU.historical.nc", "r")
    f_ssp245 = addfile(dir_ssp245+"NAJ.lat.JJA."+model(im)+".45W-0.based-maxU.ssp245.nc", "r")
    f_ssp585 = addfile(dir_ssp585+"NAJ.lat.JJA."+model(im)+".45W-0.based-maxU.ssp585.nc", "r")
    NAJ_hist = f_hist->NAJ_lat(Yrstrt1:Yrend1)
    NAJ_ssp245 = f_ssp245->NAJ_lat(Yrstrt2:Yrend2)
    NAJ_ssp585 = f_ssp585->NAJ_lat(Yrstrt2:Yrend2)
    NAJ_hist_ssp245(im,:2014-YEAR_strt) = NAJ_hist - avg(NAJ_hist(1980-YEAR_strt:2014-YEAR_strt))
    NAJ_hist_ssp245(im,2014-YEAR_strt+1:) = NAJ_ssp245 - avg(NAJ_hist(1980-YEAR_strt:2014-YEAR_strt))
    NAJ_hist_ssp585(im,:2014-YEAR_strt) = NAJ_hist - avg(NAJ_hist(1980-YEAR_strt:2014-YEAR_strt))
    NAJ_hist_ssp585(im,2014-YEAR_strt+1:) = NAJ_ssp585  - avg(NAJ_hist(1980-YEAR_strt:2014-YEAR_strt))  
end do 

    NAJ_hist_ssp245   = runave_n_Wrap(NAJ_hist_ssp245,20,0,1)  
    NAJ_hist_ssp585   = runave_n_Wrap(NAJ_hist_ssp585,20,0,1)   

    printVarSummary(NAJ_hist_ssp245)
    NAJ_noiseC = conform(NAJ_hist_ssp585, NAJ_noise, 0)
    SNR_245   = NAJ_hist_ssp245/NAJ_noiseC  
    SNR_585   = NAJ_hist_ssp585/NAJ_noiseC     

number_245 = new(nyrs, float)
number_585 = new(nyrs, float)
do it = 0, nyrs - 1 ,1 
   number_245(it) =  num(SNR_245(:,it) .gt. 1.) * 1. / nmodels
   number_585(it) =  num(SNR_585(:,it) .gt. 1.) * 1. / nmodels
end do 

YEARS = ispan(YEAR_strt, 2099, 1)
ToE = new(nmodels,integer,-9999)

do im = 0, nmodels - 1 ,1 
    ToE_possible := ind(SNR_585(im,:) .gt. 1.)
    np = dimsizes(ToE_possible)
    do i = 0, np -1,1
        ntiming = ToE_possible(i)
        if (all (SNR_585(im,ntiming::) .gt. 1.) ) then
            ToE(im) = ntiming
            break
        else
            ToE(im) = -9999
        end if
    end do
  
end do 

print(ToE + YEAR_strt+ "..." + model)
ToE_yr  = ToE + YEAR_strt
ToE_avg = avg(ToE_yr)
ToE_median = dim_median(ToE_yr)
ToE_std = dim_stddev(ToE)

range1 = ToE_median - ToE_std
range2 = ToE_median + ToE_std

print("===========================================")
; print("ToE average: "+ToE_avg)
print("ToE median: "+ToE_median)
print("ToE range 1std: "+range1 + "---" + range2)


statb = stat_dispersion(ToE_yr, True) 
ALL_bot = statb(3)  ;; lower 1/6
ALL_bot1= statb(6)  ;; lower 1/4
ALL_mean = statb(8);;statb(8) ;; median         avg(ALL_bar)
ALL_top1 = statb(10) ;; upper 1/4
ALL_top = statb(13)  ;; upper 1/6

printVarSummary(statb)
print("ToE median: "+ALL_mean)
print("ToE range 10th-90th: "+ALL_bot + "---" + ALL_top)

NAJ_hist_ssp245_MME = dim_avg_n(NAJ_hist_ssp245, 0)  
NAJ_hist_ssp585_MME = dim_avg_n(NAJ_hist_ssp585, 0)

y_top = ALL_bot1
y_mid = ALL_mean
y_low = ALL_top1
;*****************************************************   
; CALCULATE MIN & MAX
;***************************************************** 
  mnmx      = new ((/6,nyrs/), float,missing)
  mnmx(0,:) = NAJ_hist_ssp245_MME - dim_stddev_n(NAJ_hist_ssp245, 0)
  mnmx(1,:) = NAJ_hist_ssp245_MME + dim_stddev_n(NAJ_hist_ssp245, 0)
  mnmx(2,:) = NAJ_hist_ssp585_MME - dim_stddev_n(NAJ_hist_ssp585, 0)
  mnmx(3,:) = NAJ_hist_ssp585_MME + dim_stddev_n(NAJ_hist_ssp585, 0)
  mnmx(4,:) = NAJ_hist_ssp245_MME
  mnmx(5,:) = NAJ_hist_ssp585_MME
  ; mnmx(5,0:NN-1) = NAJ_obs
  mnmx(0:3,:2014-YEAR_strt) = 0.
  hist      = new ((/2,nyrs/), float,missing)
  hist(0,:) = NAJ_hist_ssp245_MME
  hist(1,:) = NAJ_hist_ssp585_MME
  hist(:,2014-YEAR_strt+2:) = missing
  mnmx(4,:2014-YEAR_strt) = missing
  mnmx(5,:2014-YEAR_strt) = missing
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic9/Fig.4.ToE")
;; for line=========
time = ispan(YEAR_strt, 2099, 1)
font                              =0.018
  res                    = True              ; plot mods desired
  res@gsnDraw            = False             ; don't draw yet
  res@gsnFrame           = False             ; don't advance frame yet

  res@vpHeightF      = 0.3               ; change aspect ratio of plot
  res@vpWidthF           = 0.4

  res@trYMinF  =  - 2.                  
  res@trYMaxF  =   3.5   
  res@trXMinF            = 1980              ; set x-axis minimum
  res@trXMaxF            = 2100              ; set x-axis minimum
  res@tmXBLabelFontHeightF      = font
  res@xyMonoLineColor    = False             ; want colored lines
  res@xyLineColors       = (/"Blue","red"/) ; colors chosen
  res@xyLineThicknesses  = (/4.,4./)      ; line thicknesses
  res@xyDashPatterns     = (/0.,0./)      ; make all lines solid

  res@tiYAxisString  = ""      ; add an axis title  
  res@tiYAxisFontHeightF =  font 
  res@txFontHeightF  = font            ; change title font heights

  res@gsnYRefLine        = (/0,NAJ_noise_MME/)  
  res@gsnYRefLineThicknesses = (/3.,3./)
  res@gsnYRefLineColor  = (/"gray30","black"/)

  res@gsnYRefLineDashPatterns = (/0,2/)
  res@tmXTOn          = False
  res@tmYROn          = False
  ;;===小刻度=========================
  res@tmXBMinorOn     = True
  res@tmXBMinorValues  = time(0::10)
  ;;for XB label=======================
  res@tmXBMode         = "Explicit"
  res@tmXBValues        = time(0::20)
  res@tmXBLabels        = time(0::20)
  res@tmXBLabelFontHeightF = font

res@tmXTBorderOn = False
res@tmYRBorderOn = False
res@xyCurveDrawOrder = "PostDraw"
res@gsnLeftString="";;"(a)NAJ lat(CMIP6)"
  top_plot = gsn_csm_xy (wks,time,mnmx(4:5,:),res)       ; create line plot
  res@gsnLeftString=""
  res@xyLineColors       = (/"black","black"/) ; colors chosen
  res@xyLineThicknesses  = (/4.,4./)      ; line thicknesses
  hist_plot = gsn_csm_xy (wks,time,hist,res)       ; create line plot
  ; Create a plot with the area between both curves filled in blue.
  delete([/res@xyLineColors,res@gsnYRefLine,res@gsnYRefLineThicknesses,res@gsnYRefLineColor,res@gsnYRefLineDashPatterns/] )
  res@gsnLeftString=""
  res@gsnXYFillColors = "LightPink"
  res@gsnXYFillOpacities = 0.2
  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  bot_plot  = gsn_csm_xy (wks,time,mnmx(2:3,:),res)  ; Create filled XY plot.

  ; Create a plot with the area between both curves filled in pink.
  res@gsnXYFillColors = "LightBlue"
  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  mid_plot  = gsn_csm_xy (wks,time,mnmx(0:1,:),res)  ; Create another filled XY plot.

  ;
  ; Overlay the top and mid plots on the bottom plot.
  ;
  ; Don't draw anything yet, because we still need to
  ; attach a legend and some titles.
  ;
  overlay(bot_plot,mid_plot)
  overlay(bot_plot,top_plot)
  overlay(bot_plot,hist_plot)
  ;*****************************************************   
  ; Manually create and attach legend
  ;***************************************************** 
    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 5.                    ; line thicker
    res_lines@gsLineColor       = "blue"                 ; line color
    xx = (/1990,2000/)
    yy = (/1.9,1.9/)
    dum1 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)              ; add polyline
    dum2 = gsn_add_text(wks,bot_plot,"ssp245",2005,1.9,res_text); add text

    yy = (/2.3,2.3/)
    res_lines@gsLineColor       = "red"                                 ; change to blue
    dum3 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)                ; add polyline
    dum4 = gsn_add_text(wks,bot_plot,"ssp585",2005,2.3,res_text)       ; add text

    yy = (/2.7,2.7/)
    res_lines@gsLineColor       = "black"                                ; change to black
    dum5 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)                ; add polyline
    dum6 = gsn_add_text(wks,bot_plot,"hist",2005,2.7,res_text) ; add text

  ; pres = True
  ; maximize_output(wks,pres)
  ;   draw(bot_plot)
  ; frame(wks)
;; for pic 2 ========
xais = ispan(1, nmodels, 1)
xlabel = new(nmodels, "string")
xlabel = sprinti("%d", xais)+""
resM                   = True                     ; plot mods desired
resM@gsnLeftString     = ""
resM@gsnDraw            = False             ; don't draw yet
resM@gsnFrame           = False             ; don't advance frame yet

resM@vpHeightF      = res@vpHeightF                ; change aspect ratio of plot
resM@vpWidthF           = res@vpWidthF 
resM@tiMainString      = ""           ; add title
resM@xyMarkLineModes   = "Markers"                ; choose which have markers
resM@xyMarkers         =  16                      ; choose type of marker  
resM@xyMarkerColor     = "red"                    ; Marker color
resM@xyMarkerOpacityF  = 0.6
resM@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
; resM@xyMarkerThicknessF = 4.
resM@tmLabelAutoStride = True                     ; nice tick mark labels
resM@tmXTBorderOn = False
resM@tmYRBorderOn = False
resM@tmXTOn          = False
resM@tmYROn          = False


; resM@xyMarkerColor = "red"
; resM@tiMainString  = ""
; resM@xyMarkers     = NhlNewMarker(wks, "^", 19, 0.0, 0.0, 1.3125, 1.5, 0.0)

resM@tmXBLabelFontHeightF      = font
resM@tiYAxisFontHeightF =  font 
resM@txFontHeightF  = font            ; change title font heights
; resM@tmXBLabelFontHeightF = 0.01

resM@trXMinF            = min(xais) - 1             ; set x-axis minimum
resM@trXMaxF            = max(xais) + 7              ; set x-axis minimum

;;===小刻度=========================
resM@tmXBMinorOn     = True
resM@tmXBMinorValues  = xais
;;for XB label=======================
resM@tmXBMode         = "Explicit"
resM@tmXBValues        = (/1,5,10,15,20,25,30,36/)
resM@tmXBLabels        = (/1,5,10,15,20,25,30,36/)
resM@tmXBLabelFontHeightF = font
; resM@tmXBLabelAngleF      = 48
resM@tmXBTickSpacingF  = 1.

 plotM  = gsn_csm_xy (wks,ispan(1,nmodels,1),ToE_yr,resM) ; create plot

resM@gsnLeftString     = "";;"(b)TOE(CMIP6 ssp585)"
 resM@xyMarkers         =  4                      ; choose type of marker  
resM@xyMarkerColor     = "black"                    ; Marker color
resM@xyMarkerOpacityF  = 1.
resM@xyMarkerThicknessF = 2.
plotM1  = gsn_csm_xy (wks,ispan(1,nmodels,1),ToE_yr,resM) ; create plot
overlay(plotM,plotM1)


;;v2--填色的阴影---
Y1 = (/y_low/)
Y2 = (/y_top/)
X1 = (/36.7/)
X2 = (/38.7/)
dum_box = new(dimsizes(Y1), graphic)
COLOR = (/"orange","red","darkgreen"/)
do i = 0, dimsizes(Y1)-1
    min_lat = Y1(i)
    max_lat = Y2(i)
    min_lon = X1(i)   ;; 71.
    max_lon = X2(i)  ;; 83.
    boxlat  = (/max_lat,max_lat,min_lat,min_lat,max_lat/)
    boxlon  = (/min_lon,max_lon,max_lon,min_lon,min_lon/) 
    gonres    = True
    gonres@gsFillColor = COLOR(i)
    gonres@gsFillOpacityF = 0.4
    gonres@gsLineDashPattern = 14
    dum_box(i) = gsn_add_polygon(wks, plotM, boxlon, boxlat, gonres)
end do
    res_text@txFontHeightF      = 0.014
    ; res_text@txAngleF           = 30
    dumL1 = gsn_add_text(wks,plotM,toint(y_top),39,y_top,res_text)       ; add text
    dumL2 = gsn_add_text(wks,plotM,toint(y_mid),39,y_mid,res_text)       ; add text
    dumL3 = gsn_add_text(wks,plotM,toint(y_low),39,y_low,res_text)       ; add text
;;v1--线---
Y1 = (/y_mid/)
Y2 = (/y_mid/)
X1 = (/36.9/)
X2 = (/38.5/)
dum_line = new(dimsizes(Y1), graphic)
COLOR = (/"red","red","darkgreen"/)
do i = 0, dimsizes(Y1)-1
    min_lat = Y1(i)
    max_lat = Y2(i)
    min_lon = X1(i)   ;; 71.
    max_lon = X2(i)  ;; 83.
    boxlat  = (/max_lat,max_lat,min_lat,min_lat,max_lat/)
    boxlon  = (/min_lon,max_lon,max_lon,min_lon,min_lon/) 
    lnres    = True
    lnres@gsLineColor = "orange"
    lnres@gsLineThicknessF = 5.
    lnres@gsLineDashPattern = 0
    dum_line(i) = gsn_add_polyline(wks, plotM, boxlon, boxlat, lnres)
end do



print(ToE_yr)


gplt = new(2, graphic)
gplt(0) = bot_plot
gplt(1) = plotM
  pres                     =True
  pres@txString     =""  
  pres@gsnPanelXWhiteSpacePercent = 1.2
  ; pres@gsnPanelYWhiteSpacePercent = 1.5
  pres@gsnPanelLabelBar    = False
  pres@gsnPanelRowSpec     = False
  pres@pmLabelBarWidthF    = 0.6
  pres@pmLabelBarHeightF   = 0.07
  pres@lbLabelFontHeightF  = 0.023
  ; pres@lbTitleString      ="units:  m^2*PVU/s"
  gsn_panel(wks,gplt,(/1,2/),pres)

end