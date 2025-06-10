load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
siglvl_a = 0.05
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 30.
lat2 = 70.
lon1 = 315.
lon2 = 360.
LEV  = 1000
lev1 = 1000
lev2 = 100  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)


dir   = "/data6/shengc2017/data/ERA5/PV_flux/season/JJA/"
fils  = systemfunc("ls "+dir + "*{"+YrStart+".."+YrEnd+"}.nc")
fs    = addfiles(fils, "r")
U     = fs[:]->U(:,{lev1:lev2},{lat1:lat2},{lon1:lon2})
U1000 = fs[:]->U(:,{LEV},:,:)
TEM     = fs[:]->TEM(:,{lev1:lev2},{lat1:lat2},{lon1:lon2})
U_NAmean = dim_avg_n_Wrap(dim_avg_n_Wrap(U,3),0)
U1000_NAmean = dim_avg_n_Wrap(U(:,{LEV},:,:),2)

printVarSummary(U_NAmean)
gradLatLon  = grad_latlon_cfd(TEM, TEM&latitude, TEM&longitude, True, False )   
MTG  = gradLatLon[0] * -1.  ;; -ptpy
MTG_NAmean = dim_avg_n_Wrap(MTG, 3)
MTG_NAmean_clim  = dim_avg_n_Wrap(MTG_NAmean, 0)   
copy_VarMeta(TEM(:,:,:,0),MTG_NAmean)
; U due to thermal wind 
plev = (TEM&level) * 100.
print(plev)
nz   = dimsizes(plev)    
lnp  = log(plev)
ptop = min(lnp)
dlnp   = dpres_plevel(lnp, log(100000), ptop, 0)  
dlnpc = conform(MTG_NAmean, dlnp, 1)   
print(dlnp)

; ;;===================
lat  = TEM&latitude
latc = conform(MTG_NAmean, lat, 2)
fprm  = coriolis_param(lat)
fprmc = conform(MTG_NAmean, fprm, 2)
fprmc@_FillValue = TEM@_FillValue
; fprmc = where(fprmc .eq. 0 .or. latc .le. 5., fprmc@_FillValue, fprmc) 
fprmc = where(fprmc .eq. 0 , fprmc@_FillValue, fprmc) 
Rprm    = 287.  ;; J/kg/K
kprm  = Rprm/fprmc
copy_VarMeta(fprmc, kprm)
Ut   = MTG_NAmean
Ut   = 0.
do iz = 0, nz - 1
Ut(:,iz,:) = tofloat(dim_sum_n_Wrap(MTG_NAmean(:,0:iz,:) * dlnpc(:,0:iz,:) * kprm(:,0:iz,:) , 1)  )
end do   
U1000_NAmean_C = conform(Ut,U1000_NAmean,(/0,2/))
Usum = Ut+U1000_NAmean_C
copy_VarMeta(TEM(:,:,:,0),Usum)
dir   = "/data6/shengc2017/data/ERA5/Eddy_Evector/season_dueto_daily/JJA/"
fils  = systemfunc("ls "+dir + "*{"+YrStart+".."+YrEnd+"}.nc")
fs    = addfiles(fils, "r")
Ex     = fs[:]->Ex(:,{lev1:lev2},{lat1:lat2},{320:350})
Ey     = fs[:]->Ey(:,{lev1:lev2},{lat1:lat2},{320:350})
Ep     = fs[:]->Ep(:,{lev1:lev2},{lat1:lat2},{320:350})

meta   = Ex
GradLatLon = grad_latlon_cfd(Ex, Ex&latitude, Ex&longitude, False, False)
Ex_diff = GradLatLon[1]
GradLatLon = grad_latlon_cfd(Ey, Ey&latitude, Ey&longitude, False, False)
Ey_diff = GradLatLon[0]
printVarSummary(Ex_diff)
;divH = uv2dv_cfd(Ex, Ey, Ex&latitude, Ex&longitude, 2)
Ep_diff = center_finite_diff_n(Ep, (Ep&level)*100., False, 0, 1)
divH = (   Ex_diff + Ey_diff  )  
copy_VarMeta(meta,divH)

Ey_NAmean = dim_avg_n_Wrap(Ey(:,{lev1:lev2},{lat1:lat2},{320:350}),3)
Ep_NAmean = dim_avg_n_Wrap(Ep(:,{lev1:lev2},{lat1:lat2},{320:350}),3) 
divH_NAmean=dim_avg_n_Wrap(divH(:,{lev1:lev2},{lat1:lat2},{320:350}),3)
copy_VarMeta(Ey_NAmean,Ep_NAmean)
Ey_NAmean=runave_n_Wrap(Ey_NAmean, 9, 1, 0)
Ep_NAmean=runave_n_Wrap(Ep_NAmean, 9, 1, 0)
divH_NAmean=runave_n_Wrap(divH_NAmean, 9, 1, 0)
;;Start making regression--------------------
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,MTG_NAmean,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_MTG_NAmean = student_t(tval, df)
MTG_NAmean_trd = RC * 10. * 10^6.
copy_VarMeta(MTG_NAmean(0,:,:), MTG_NAmean_trd)
copy_VarMeta(MTG_NAmean(0,:,:), prob_MTG_NAmean)
MTG_NAmean_trd = where(ismissing(MTG_NAmean(0,:,:)), MTG_NAmean@_FillValue, MTG_NAmean_trd)
prob_MTG_NAmean = where(ismissing(MTG_NAmean(0,:,:)), MTG_NAmean@_FillValue, prob_MTG_NAmean)
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,Ut,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_Ut = student_t(tval, df)
Ut_trd = RC * 10. 
copy_VarMeta(Ut(0,:,:), Ut_trd)
copy_VarMeta(Ut(0,:,:), prob_Ut)
Ut_trd = where(ismissing(Ut(0,:,:)), Ut@_FillValue, Ut_trd)
prob_Ut = where(ismissing(Ut(0,:,:)), Ut@_FillValue, prob_Ut)
; ;**********regression*************************

x  = (TIME-1980) * 1.
RC := regCoef_n(x,Usum,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_Usum = student_t(tval, df)
Usum_trd = RC * 10. 
copy_VarMeta(Usum(0,:,:), Usum_trd)
copy_VarMeta(Usum(0,:,:), prob_Usum)
Usum_trd = where(ismissing(Usum(0,:,:)), Usum@_FillValue, Usum_trd)
prob_Usum = where(ismissing(Usum(0,:,:)), Usum@_FillValue, prob_Usum)
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,U1000,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_U1000 = student_t(tval, df)
U1000_trd = RC * 10. 
copy_VarMeta(U1000(0,:,:), U1000_trd)
copy_VarMeta(U1000(0,:,:), prob_U1000)
U1000_trd = where(ismissing(U1000(0,:,:)), U1000@_FillValue, U1000_trd)
prob_U1000 = where(ismissing(U1000(0,:,:)), U1000@_FillValue, prob_U1000)

; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,Ey_NAmean,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_Ey_NAmean = student_t(tval, df)
Ey_NAmean_trd = RC * 10. * 10.
copy_VarMeta(Ey_NAmean(0,:,:), Ey_NAmean_trd)
copy_VarMeta(Ey_NAmean(0,:,:), prob_Ey_NAmean)
Ey_NAmean_trd = where(ismissing(Ey_NAmean(0,:,:)), Ey_NAmean@_FillValue, Ey_NAmean_trd)
prob_Ey_NAmean = where(ismissing(Ey_NAmean(0,:,:)), Ey_NAmean@_FillValue, prob_Ey_NAmean)
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,Ep_NAmean,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_Ep_NAmean = student_t(tval, df)
Ep_NAmean_trd = RC * 10. * 10. * 10.  
copy_VarMeta(Ep_NAmean(0,:,:), Ep_NAmean_trd)
copy_VarMeta(Ep_NAmean(0,:,:), prob_Ep_NAmean)
Ep_NAmean_trd = where(ismissing(Ep_NAmean(0,:,:)), Ep_NAmean@_FillValue, Ep_NAmean_trd)
prob_Ep_NAmean = where(ismissing(Ep_NAmean(0,:,:)), Ep_NAmean@_FillValue, prob_Ep_NAmean)

; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,divH_NAmean,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_divH_NAmean = student_t(tval, df)
divH_NAmean_trd = RC * 10. * 10000. * 10.
copy_VarMeta(divH_NAmean(0,:,:), divH_NAmean_trd)
copy_VarMeta(divH_NAmean(0,:,:), prob_divH_NAmean)
divH_NAmean_trd = where(ismissing(divH_NAmean(0,:,:)), divH_NAmean@_FillValue, divH_NAmean_trd)
prob_divH_NAmean = where(ismissing(divH_NAmean(0,:,:)), divH_NAmean@_FillValue, prob_divH_NAmean)
Ey_NAmean_trd= smth9_Wrap(Ey_NAmean_trd, 0.5, 0.25, False)
Ep_NAmean_trd = smth9_Wrap(Ep_NAmean_trd, 0.5, 0.25, False)
;;======================
wks   = gsn_open_wks ("eps", "/data6/shengc2017/work9_sub_Jet/NCC/pic1/Fig.E6.obs.trd.SST.SSTG.Ut.eddyDrivenJet.v3")
;;;==zonal wind jet==========================specify the contour==============================================================
resh = True
resh@gsnDraw                       = False
resh@gsnFrame                      = False
resh@gsnMaximize                   = False
resh@cnLevelSelectionMode          = "ExplicitLevels"
resh@gsnLeftString                 = ""
resh@gsnRightString                = ""
resh@cnFillOn                      = False
resh@cnFillPalette                 = "cmp_b2r"  

resh@cnLevels                      = (/10.,15.,20.,24.,30.,35.,40.,45.,50./)
resh@cnLinesOn                     = True
resh@cnMonoLineColor               = True
resh@cnLineColor                   = "black"
resh@cnLineThicknessF              = 2.
;;for line label==============================================================
resh@cnLineLabelsOn                =True
resh@cnLineLabelFontHeightF        =0.018  
resh@cnLineLabelFontColor          ="black"
resh@cnLineLabelBackgroundColor    ="white"
resh@cnInfoLabelOn                 = False
resh@cnLineLabelPlacementMode         = "Computed"
resh@cnLineLabelInterval              =1
resh@cnLineLabelDensityF              =0.5
resh@gsnPresHgtHeightLabelOn=False

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
resh2@cnLevels                      =  (/-0.4,-0.3,0.2,0.3,0.4/)
resh2@cnLinesOn                     = True
resh2@cnMonoLineColor               = False
resh2@cnLineColors                   = (/"blue","blue","red","red","red"/)
; resh2@cnLineDashPatterns             = 16
resh2@cnLineThicknessF              = 1.5
;;for line label==============================================================
resh2@cnLineLabelsOn                =True
resh2@cnLineLabelFontHeightF        =0.012  
resh2@cnLineLabelFormat             = "0@;*.1f"       
resh2@cnLineLabelFontColor          ="black"
resh2@cnLineLabelBackgroundColor    ="gray"
resh2@cnInfoLabelOn                 = False
resh2@cnLineLabelPlacementMode         = "Computed"
resh2@cnLineLabelInterval              =1
resh2@cnLineLabelDensityF              =0.5
resh2@cnLabelDrawOrder              = "PostDraw"
;;for confidence lev dot=====================================================================
opt = True
opt@gsnDraw                       =False
opt@gsnFrame                      =False
opt@gsnMaximize                   =False
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

;;   
plot   = new(4, graphic)
dplot  = new(4, graphic)
rest                               =True
; rest@gsnPolar                      ="NH"
font                              =0.03
rest@gsnMaximize                   =False
rest@gsnDraw                       =False
rest@gsnFrame                      =False
rest@tmXTOn                        =False
rest@tmYROn                        =False 
ovalue = 0.

rest@gsnRightString                = ""
rest@gsnLeftStringOrthogonalPosF   = ovalue
rest@gsnRightStringOrthogonalPosF  = ovalue
rest@gsnLeftStringFontHeightF      = font
rest@gsnRightStringFontHeightF     = font
rest@tmYLLabelFontHeightF          = font
rest@tmXBLabelFontHeightF          = font
rest@tmXBTickSpacingF              = 10

;; for contour=============================================================
rest@cnFillPalette                 ="cmp_b2r"    ;;BlueDarkRed18   cmp_b2r  GMT_polar
rest@cnFillOn                      =True
rest@cnLinesOn                     =False
rest@cnLineThicknessF              = 0.0000001 
rest@cnLineColor                   = "blue"
rest@gsnContourZeroLineThicknessF  = 3.
;;for line label==============================================================
rest@cnLineLabelsOn                =False
rest@cnLineLabelFontHeightF        =0.02  
rest@cnLineLabelFontColor          ="black"
rest@cnLineLabelBackgroundColor    ="white"
rest@cnInfoLabelOn                 = False
rest@cnLineLabelPlacementMode         = "Computed"
rest@cnLineLabelInterval              =2
rest@cnLineLabelDensityF              =0.5

;;for color bar===========================================================
rest@lbLabelBarOn                  = True
; rest@pmLabelBarOrthogonalPosF      = 0.15
rest@lbLabelFontHeightF            = font
rest@pmLabelBarHeightF             = 0.1
rest@gsnPresHgtHeightLabelOn=False
opt@gsnPresHgtHeightLabelOn=False
rest@tiYAxisString = ""
opt@tiYAxisString = ""

max_vaule = 0.6
interval  = 0.1
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1
colors := read_colormap_file(rest@cnFillPalette)
icol   := span_color_indexes(colors(2:,:), ncols)  
icol = icol+2  
icol(midcol1:midcol2) = -1
rest@cnLevelSelectionMode = "ManualLevels"   
rest@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
rest@cnMaxLevelValF       = max_vaule               
rest@cnLevelSpacingF      = interval         ;;850->2
rest@cnFillColors := icol
rest@gsnLeftString = "";;"(d)trd of MTG"
plot(0) = gsn_csm_pres_hgt(wks, MTG_NAmean_trd, rest)
dplot(0) = gsn_csm_pres_hgt(wks, prob_MTG_NAmean, opt)


max_vaule = 1.4
interval  = 0.2
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1
colors := read_colormap_file(rest@cnFillPalette)
icol   := span_color_indexes(colors(2:,:), ncols)  
icol = icol+2  
icol(midcol1:midcol2) = -1
rest@cnLevelSelectionMode = "ManualLevels"   
rest@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
rest@cnMaxLevelValF       = max_vaule               
rest@cnLevelSpacingF      = interval         ;;850->2
rest@cnFillColors := icol

rest@gsnLeftString = "";;"(g)trd of Ut + U1000"
plot(1) = gsn_csm_pres_hgt(wks, Usum_trd, rest)
dplot(1) = gsn_csm_pres_hgt(wks, prob_Usum, opt)

hplot1 = gsn_csm_pres_hgt(wks,  U_NAmean, resh)
hplot2 = gsn_csm_pres_hgt(wks,  U_NAmean, resh)


max_vaule = 1.4
interval  = 0.2
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1
colors := read_colormap_file(rest@cnFillPalette)
icol   := span_color_indexes(colors(2:,:), ncols)  
icol = icol+2  
icol(midcol1:midcol2) = -1
rest@cnLevelSelectionMode = "ManualLevels"   
rest@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
rest@cnMaxLevelValF       = max_vaule               
rest@cnLevelSpacingF      = interval         ;;850->2
rest@cnFillColors := icol
rest@gsnLeftString = "";;"(f)trd of Ut"
plot(2) = gsn_csm_pres_hgt(wks, Ut_trd, rest)
dplot(2) = gsn_csm_pres_hgt(wks, prob_Ut, opt)


max_vaule = 0.8
interval  = 0.1
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1
colors := read_colormap_file(rest@cnFillPalette)
icol   := span_color_indexes(colors(2:,:), ncols)  
icol = icol+2  
icol(midcol1:midcol2) = -1
rest@cnLevelSelectionMode = "ManualLevels"   
rest@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
rest@cnMaxLevelValF       = max_vaule               
rest@cnLevelSpacingF      = interval         ;;850->2
rest@cnFillColors := icol

rest@gsnLeftString = "";;"(e)trd of div(Ex,Ey)"
plot(3) = gsn_csm_pres_hgt(wks, divH_NAmean_trd, rest)
dplot(3) = gsn_csm_pres_hgt(wks, prob_divH_NAmean, opt)


;;for vector---
;;for vector ===================== 
vres                          = True
vres@gsnDraw                  = False
vres@gsnFrame                 = False
vres@gsnMaximize              = False
vres@gsnLeftString            = ""
vres@gsnRightString           = ""
vres@vcMinDistanceF           = 0.03           
vres@vcMinFracLengthF         = 0.3
vres@vcMinMagnitudeF          = 5.   

;;=========================================
;;设定LineArrow----------------------------
vres@vcGlyphStyle             = "CurlyVector"  ;;LineArrow
vres@vcMonoLineArrowColor     = True         
vres@vcLineArrowColor         = "black"
vres@vcLineArrowThicknessF    = 1.           

    vres@vcLineArrowHeadMaxSizeF = 0.012
    vres@vcLineArrowHeadMinSizeF = 0.012
;; set ref-------------------------------- 
value = 5.
vres@vcRefAnnoOn             = True
vres@vcRefMagnitudeF         = value            
vres@vcRefAnnoString1        = tostring_with_format(value, "%2.1f")
vres@vcRefAnnoString2On      = False
vres@vcRefLengthF            = 0.03          ; define length of vec ref
vres@vcRefAnnoOrthogonalPosF = -0.3
vres@vcRefAnnoParallelPosF   = 0.999
vres@vcRefAnnoSide           = "Bottom"
vres@vcRefAnnoBackgroundColor = "grey"
vres@vcRefAnnoFontHeightF    = 0.012
vres@vcRefAnnoPerimOn        = False
vplot = gsn_csm_vector(wks, Ey_NAmean_trd, Ep_NAmean_trd, vres)
;overlay(plot(4),vplot)
do i = 0, 3
overlay(plot(i), dplot(i))
end do

overlay(plot(1),hplot1)
overlay(plot(2),hplot2)


dum1 = new(4, graphic)
lat_line := (/1000.,100./)   ;;65  (/40.,50,60/)
lon_line := (/50.,50./)   ;;150
lnres = True
lnres@gsLineColor = "black"
lnres@gsLineThicknessF = 1.5
lnres@gsLineDashPattern = 14
do i = 0, 3
dum1(i) = gsn_add_polyline(wks, plot(i), lon_line, lat_line, lnres)  
end do 


dum2 = new(4, graphic)
lat_line := (/250.,250./)   ;;65  (/40.,50,60/)
lon_line := (/30.,70./)   ;;150
lnres = True
lnres@gsLineColor = "black"
lnres@gsLineThicknessF = 1.5
lnres@gsLineDashPattern = 14
do i = 0, 3
dum2(i) = gsn_add_polyline(wks, plot(i), lon_line, lat_line, lnres)  
end do 
gplt = new(4, graphic)
gplt(0) = plot(0) 
gplt(1) = plot(1)
gplt(2) = plot(2)
gplt(3) = plot(3)

pres                     =True
pres@gsnPanelDebug      = False                   
pres@txString            =""
pres@gsnPanelXWhiteSpacePercent = 1.2
pres@gsnPanelYWhiteSpacePercent = 5
pres@gsnPanelRowSpec     = False
gsn_panel(wks,gplt,(/2,2/),pres)
wallClockElapseTime(TotalTime, "totaltime", 0)  

end