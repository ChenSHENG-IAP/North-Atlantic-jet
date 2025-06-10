load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
siglvl_a = 0.05 
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 20.
lat2 = 70.
lon1 = 0.
lon2 = 360.
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

; latent heat flux
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/surface_upward_heat_flux/post_cat/hfls/"
fE    = addfile(dir+"hfls_Amon_MPI-ESM_historical_rcp85_E15.ensmean_185001-209912.v2.nc", "r")
fP    = addfile(dir+"hfls_Amon_MPI-ESM_historical_rcp85_P15.ensmean_185001-209912.v2.nc", "r")

time   := fE->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
LH_Egrp   = fE->hfls(iYYYY,:,:)
time   := fP->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
LH_Pgrp   = fP->hfls(iYYYY,:,:)

LH_Egrp_JJA    = month_to_season(LH_Egrp, "JJA")
LH_Pgrp_JJA    = month_to_season(LH_Pgrp, "JJA")

;; sensible heat flux 
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/surface_upward_heat_flux/post_cat/hfss/"
fE    = addfile(dir+"hfss_Amon_MPI-ESM_historical_rcp85_E15.ensmean_185001-209912.v2.nc", "r")
fP    = addfile(dir+"hfss_Amon_MPI-ESM_historical_rcp85_P15.ensmean_185001-209912.v2.nc", "r")

time   := fE->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
SH_Egrp   = fE->hfss(iYYYY,:,:)
time   := fP->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
SH_Pgrp   = fP->hfss(iYYYY,:,:)

SH_Egrp_JJA    = month_to_season(SH_Egrp, "JJA")
SH_Pgrp_JJA    = month_to_season(SH_Pgrp, "JJA")

LH     = (LH_Egrp_JJA - LH_Pgrp_JJA)*-1.  ;;after checking with ERA5，I confirm that the upward flux is positive (see /home/users/shengc01/work11/scripts7/compare.heatflux.ERA5.MPIMME.ncl)
SH     = (SH_Egrp_JJA - SH_Pgrp_JJA)*-1.
SUM       = LH+SH
copy_VarMeta(SH_Egrp_JJA, LH)
copy_VarMeta(SH_Egrp_JJA, SH)
copy_VarMeta(SH_Egrp_JJA,SUM)


dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ts/"
fils  = dir+"ts_Amon_MPI-ESM_historical_rcp85_P15.ensmean_185001-209912.v2.nc"
fP    = addfile(fils,"r") 
time   := fP->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
tsP15     = fP->ts(iYYYY,:,:)
tsP15_JJA = month_to_season(tsP15, "JJA")

fils  = dir+"ts_Amon_MPI-ESM_historical_rcp85_E15.ensmean_185001-209912.v2.nc"
fE    = addfile(fils,"r") 
time   := fE->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
tsE15     = fE->ts(iYYYY,:,:)
tsE15_JJA = month_to_season(tsE15, "JJA")
meta   = tsE15(0,:,:)
SST_diff     = tsE15_JJA - tsP15_JJA 
rad = 3.1415926/180.
rad    = 4.0*atan(1.0)/180.0
re     = 6371220.0
rr     = re*rad
; SSTG_diff    = -1. * center_finite_diff_n(SST_diff, tsE15_JJA&lat, False, 0, 1)
; ;SSTG_diff    = -1. * center_finite_diff_n(SST_diff, tsE15_JJA&lat*rr, False, 0, 1)
GradLatLon = grad_latlon_cfd(SST_diff, tsE15_JJA&lat, tsE15_JJA&lon, False, False)
SSTG_diff = GradLatLon[0] * -1.
copy_VarMeta(tsE15_JJA,SSTG_diff)
copy_VarMeta(tsE15_JJA,SST_diff)

;;for mask --
; assume data is 3D (time,lat,lon)
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata,SSTG_diff&lat,SSTG_diff&lon)
; lsm is a 2D array, in order to use it in mask, we must conform it
; to the size of the 3D array "data". 
SSTG_diff = mask(SSTG_diff,conform(SSTG_diff,lsm,(/1,2/)).ge.1,False)
SST_diff  = mask(SST_diff,conform(SST_diff,lsm,(/1,2/)).ge.1,False)

;;Start making regression--------------------

x  = (TIME-1980) * 1.
RC := regCoef_n(x,SSTG_diff,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SSTG_diff = student_t(tval, df)
SSTG_diff_trd = RC * 10. *1000.*1000. ;; (K per 1000 km per 10a = 10-6 K per m per 10a)
copy_VarMeta(SST_diff(0,:,:), SSTG_diff_trd)
copy_VarMeta(SST_diff(0,:,:), prob_SSTG_diff)
SSTG_diff_trd = where(ismissing(SST_diff(0,:,:)), SST_diff@_FillValue, SSTG_diff_trd)
prob_SSTG_diff = where(ismissing(SST_diff(0,:,:)), SST_diff@_FillValue, prob_SSTG_diff)

x  = (TIME-1980) * 1.
RC := regCoef_n(x,SST_diff,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SST_diff = student_t(tval, df)
SST_diff_trd = RC * 10. 
copy_VarMeta(SST_diff(0,:,:), SST_diff_trd)
copy_VarMeta(SST_diff(0,:,:), prob_SST_diff)
SST_diff_trd = where(ismissing(SST_diff(0,:,:)), SST_diff@_FillValue, SST_diff_trd)
prob_SST_diff = where(ismissing(SST_diff(0,:,:)), SST_diff@_FillValue, prob_SST_diff)


x  = TIME-1980
RC := regCoef_n(x,SUM,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SUM = student_t(tval, df)
SUM_trd = RC * 10. 
copy_VarMeta(SUM(0,:,:), SUM_trd)
copy_VarMeta(SUM(0,:,:), prob_SUM)
SUM_trd = smth9_Wrap(SUM_trd, 0.5, 0.25, True)
;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic7/Fig.2.trd.SST.SSTG.heatflux.MPI.E15-P15.v2")
res                               =True
;res@gsnPolar                      ="NH"
font                              =0.03
res@gsnMaximize                   =False
res@gsnDraw                       =False
res@gsnFrame                      =False
res@tmXTOn                        =False
res@tmYROn                        =False 
ovalue = 0.03

res@gsnRightString                = ""
res@gsnLeftStringOrthogonalPosF   = ovalue
res@gsnRightStringOrthogonalPosF  = ovalue
res@gsnLeftStringFontHeightF      = font
res@gsnRightStringFontHeightF     = font
res@tmYLLabelFontHeightF          = font
res@tmXBLabelFontHeightF          = font

res@tmXBTickSpacingF              = 20
res@tmYLTickSpacingF              = 20

res@gsnAddCyclic                  = True
res@mpMinLatF                     =lat1
res@mpMaxLatF                     =lat2
res@mpMinLonF                     =-90
res@mpMaxLonF                     =20 ;;

;; for contour=============================================================
res@cnFillPalette                 ="cmp_b2r"   ;;GMT_polar  cmp_b2r   cmp_b2r
res@cnFillOn                      =True
res@cnLinesOn                     =False
res@cnLineThicknessF              = 1.5
res@cnMonoLineColor               = False
;res@cnLineColor                   = "white"
;res@gsnContourZeroLineThicknessF  = 3.
;;for line label==============================================================
res@cnLineLabelsOn                =False
res@cnLineLabelFontHeightF        =0.01  
res@cnLineLabelFontColor          ="black"
res@cnLineLabelBackgroundColor    ="white"
res@cnInfoLabelOn                 = False
res@cnLineLabelPlacementMode         = "Computed"
res@cnLineLabelInterval              =1
res@cnLineLabelDensityF              =0.5
;;for color bar===========================================================
res@lbLabelBarOn                  = True
res@pmLabelBarOrthogonalPosF      = 0.23
;res@pmLabelBarHeightF             = 0.07
res@lbLabelFontHeightF            = 0.025

;;;==zonal wind jet==========================specify the contour==============================================================
resh2 = True
resh2@gsnDraw                       = False
resh2@gsnFrame                      = False
resh2@gsnMaximize                   = False
resh2@cnLevelSelectionMode          = "ExplicitLevels"
resh2@gsnLeftString                 = ""
resh2@gsnRightString                = ""
resh2@cnFillOn                      = False
resh2@cnFillPalette                 = "MPL_PRGn"  
; resh2@cnLevels                      =  (/15.,20.,25.,30,40,50,60,70./)
resh2@cnLevels                      =   (/-2.,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0/) ;(/-0.3,-0.2,0.2,0.3,0.4/)
resh2@cnLinesOn                     = True
resh2@cnMonoLineColor               = False 
resh2@cnLineColors                   = (/"blue","blue","blue","blue","red","red","red","red"/)
; resh2@cnLineDashPatterns             = 16
resh2@cnLineThicknessF              = 1.
;;for line label==============================================================
resh2@cnLineLabelsOn                =True
resh2@cnLineLabelFontHeightF        =0.012  
resh2@cnLineLabelFormat             = "0@;*.2f"       
resh2@cnLineLabelFontColor          ="black"
resh2@cnLineLabelBackgroundColor    ="gray"
resh2@cnInfoLabelOn                 = False
resh2@cnLineLabelPlacementMode         = "Computed"
resh2@cnLineLabelInterval              =1.
resh2@cnLineLabelDensityF              =1.
resh2@cnLabelDrawOrder              = "PostDraw"
;;for confidence lev dot=====================================================================
opt = True
opt@gsnDraw                       =False
opt@gsnFrame                      =False
opt@gsnLeftString          =""
opt@gsnRightString         =""
opt@lbLabelBarOn           =False
opt@cnFillOn               = True
opt@cnLinesOn              = False               ; turn off contour lines
opt@cnLineLabelsOn         = False
opt@cnLevelSelectionMode   = "ExplicitLevels"      ; manually set cn levels
opt@cnLevels               = (/siglvl_a/)    ; 36y 95%
opt@cnFillPattern          = 17
opt@cnFillScales           = (/0.1,0.75/)    ; change densities  
opt@cnFillColors           = (/"white" ,"Transparent"/)
opt@cnFillDotSizeF         =  0.001
opt@cnInfoLabelOn          = False
;;============================================
plot = new(2, graphic)
dplot = new(2, graphic)
;;fig.a,b
max_vaule = 0.16
interval  = 0.02
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1

colors = read_colormap_file(res@cnFillPalette)
icol   = span_color_indexes(colors(2:,:), ncols) 
icol = icol+2  
icol(midcol1:midcol2) = -1
res@cnLevelSelectionMode = "ManualLevels"   
res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
res@cnMaxLevelValF       = max_vaule               
res@cnLevelSpacingF      = interval         ;;850->2
res@cnFillColors = icol

res@gsnRightString  = ""
res@gsnLeftString   = "" ;;"(a)trd of SST (K/10a)"
plot(0)   = gsn_csm_contour_map(wks, SST_diff_trd, res)  
dplot(0)  = gsn_csm_contour(wks, prob_SST_diff , opt)

max_vaule = 0.4
interval  = 0.05
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1
colors := read_colormap_file(res@cnFillPalette)
icol   := span_color_indexes(colors(2:,:), ncols) 
icol = icol+2  
icol(midcol1:midcol2) = -1
res@cnLevelSelectionMode = "ManualLevels"   
res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
res@cnMaxLevelValF       = max_vaule               
res@cnLevelSpacingF      = interval         ;;850->2
res@cnFillColors := icol
res@gsnRightString  = ""
res@gsnLeftString   = "";;"(b)trd of SSTG (10~S~-1~N~K/deg 10a~S~-1~N~)"
plot(1)   = gsn_csm_contour_map(wks, SSTG_diff_trd, res)  
dplot(1)  = gsn_csm_contour(wks, prob_SSTG_diff , opt)

do i = 0, 1
  overlay(plot(i), dplot(i))
end do

hplot   = gsn_csm_contour(wks, SUM_trd, resh2)
overlay(plot(0),hplot)  





lat_line = (/50.,50.,50./)   ;;65  (/40.,50,60/)
lon_line = (/-90.,-60.,20./)   ;;150
lnres = True
lnres@gsLineColor = "black"
lnres@gsLineThicknessF = 1.5
lnres@gsLineDashPattern = 14
dum1 = new(2, graphic)
do i = 0, 1
dum1(i) = gsn_add_polyline(wks, plot(i), lon_line, lat_line, lnres)  
end do   

; ;; Add box=======================================================
  min_lat = lat1
  max_lat = lat2
  min_lon = 315.   ;; 71.
  max_lon = 360.  ;; 83.
  boxlat  = (/max_lat,max_lat,min_lat,min_lat,max_lat/)
  boxlon  = (/min_lon,max_lon,max_lon,min_lon,min_lon/) 

  lnres    = True
  lnres@gsLineColor = "black"
  lnres@gsLineThicknessF = 1.0
  lnres@gsLineDashPattern = 0
  dum = gsn_add_polyline(wks, plot(1), boxlon, boxlat, lnres) 
    
gplt = new(2, graphic)
gplt(0) = plot(0)
gplt(1) = plot(1)


pres                     =True
pres@txString            =""
pres@gsnPanelXWhiteSpacePercent = 1.2
pres@gsnPanelYWhiteSpacePercent = 1.5
pres@gsnPanelLabelBar    = False
pres@gsnPanelRowSpec     = False
pres@pmLabelBarWidthF    = 0.6
pres@pmLabelBarHeightF   = 0.07
pres@lbLabelFontHeightF  = 0.013
gsn_panel(wks,gplt,(/1,2/),pres)
end
