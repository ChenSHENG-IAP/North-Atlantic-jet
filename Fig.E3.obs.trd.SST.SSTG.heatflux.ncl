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
lon1 = 0.
lon2 = 360.
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

dir   = "/home/users/shengc01/data/SST/"
f     = addfile(dir + "sst.mnmean.COBE.nc", "r")
;f     = addfile(dir + "HadISST_sst.nc", "r")   ;;sst.mnmean.v5.nc
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
SSTO   = f->sst(iYYYY,:,:)
SST    = month_to_season(SSTO, "JJA")


GradLatLon = grad_latlon_cfd(SST, SST&lat, SST&lon, False, False)
SSTG = GradLatLon[0] * -1.
copy_VarMeta(SST,SSTG)

dir   = "/home/users/shengc01/data/ERA5/"
f    = addfile(dir+"ERA5.surface.heatflux.1980.2023.nc", "r")

time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
mslhf   = f->mslhf(iYYYY,:,:)
msshf   = f->msshf(iYYYY,:,:)
mslhf_flt = short2flt(mslhf)
msshf_flt = short2flt(msshf)
EFLUX    = month_to_season(mslhf_flt, "JJA")
HFLUX    = month_to_season(msshf_flt, "JJA")
SUM       = EFLUX+HFLUX
copy_VarMeta(EFLUX,SUM)

;;for mask --
; assume data is 3D (time,lat,lon)
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata,SUM&latitude,SUM&longitude)
; lsm is a 2D array, in order to use it in mask, we must conform it
; to the size of the 3D array "data". 
SUM = mask(SUM,conform(SUM,lsm,(/1,2/)).ge.1,False)


; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,SSTG,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SSTG = student_t(tval, df)
SSTG_trd = RC * 10. * 10^6.
copy_VarMeta(SST(0,:,:), SSTG_trd)
copy_VarMeta(SST(0,:,:), prob_SSTG)
SSTG_trd = where(ismissing(SST(0,:,:)), SST@_FillValue, SSTG_trd)
prob_SSTG = where(ismissing(SST(0,:,:)), SST@_FillValue, prob_SSTG)
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,SST,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SST = student_t(tval, df)
SST_trd = RC * 10. 
copy_VarMeta(SST(0,:,:), SST_trd)
copy_VarMeta(SST(0,:,:), prob_SST)
SST_trd = where(ismissing(SST(0,:,:)), SST@_FillValue, SST_trd)
prob_SST = where(ismissing(SST(0,:,:)), SST@_FillValue, prob_SST)

; ;**********regression*************************
x  = TIME-1980
RC := regCoef_n(x,SUM,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SUM = student_t(tval, df)
Trd_SUM = RC * 10. * -1.
copy_VarMeta(SUM(0,:,:), Trd_SUM)
copy_VarMeta(SUM(0,:,:), prob_SUM)
;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic7/Fig.E3.obs.trd.SST.SSTG.heatflux")
res                               =True
;res@gsnPolar                      ="NH"
font                              =0.023
res@gsnMaximize                   =False
res@gsnDraw                       =False
res@gsnFrame                      =False
res@tmXTOn                        =False
res@tmYROn                        =False 
ovalue = 0.02

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
res@mpMaxLonF                     =20

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
res@pmLabelBarOrthogonalPosF      = 0.19
res@pmLabelBarHeightF             = 0.07
res@lbLabelFontHeightF            = 0.021

;;;==zonal wind jet==========================specify the contour===========
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
max_vaule = 0.5
interval  = 0.05
ncols = toint(2.*max_vaule/interval + 2.)
midcol2 = toint(ncols/2.) 
midcol1 = midcol2 - 1

colors = read_colormap_file(res@cnFillPalette)
icol   = span_color_indexes(colors(2:,:), ncols) 
icol = icol+2  
;icol(midcol1:midcol2) = -1
res@cnLevelSelectionMode = "ManualLevels"   
res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
res@cnMaxLevelValF       = max_vaule               
res@cnLevelSpacingF      = interval         ;;850->2
res@cnFillColors = icol
;;=============================================
res@gsnRightString  = ""

;;    
plot   = new(3, graphic)
hplot2 = new(3, graphic)
dplot  = new(3, graphic)
;;fig.a,b
res@gsnRightString  = ""
res@gsnLeftString   = "";;"(a)trd of SST (K/10a)"
plot(0)   = gsn_csm_contour_map(wks, SST_trd, res)  
dplot(0)  = gsn_csm_contour(wks, prob_SST , opt)

max_vaule = 1.
interval  = 0.1
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
res@gsnLeftString   = "";;"(b)trd of SSTG (10~S~-1~N~K/10deg 10a~S~-1~N~)"
plot(1)   = gsn_csm_contour_map(wks, SSTG_trd, res)  
dplot(1)  = gsn_csm_contour(wks, prob_SSTG, opt)
; hplot1    = gsn_csm_contour(wks,U1000_trd,resh2)
; overlay(plot(1),hplot1)

max_vaule = 5.
interval  = 1.
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
res@gsnLeftString   = "";;"(c)trd of turbulent heat flux (W m~S~-2~N~ 10a~S~-1~N~)"
plot(2)   = gsn_csm_contour_map(wks, Trd_SUM, res)  
dplot(2)  = gsn_csm_contour(wks, prob_SUM , opt)

do i = 0, 2
  overlay(plot(i), dplot(i))
end do


lat_line = (/50.,50.,50./)   ;;65  (/40.,50,60/)
lon_line = (/-90.,-60.,20./)   ;;150
lnres = True
lnres@gsLineColor = "black"
lnres@gsLineThicknessF = 1.5
lnres@gsLineDashPattern = 14
dum1 = new(3, graphic)
do i = 0, 2
dum1(i) = gsn_add_polyline(wks, plot(i), lon_line, lat_line, lnres)  
end do    


gplt = new(3, graphic)
gplt(0) = plot(0)
gplt(1) = plot(1)
gplt(2) = plot(2)
; gplt(3) = plot(3)
    pres                     =True
    pres@txString            =""
    pres@gsnPanelXWhiteSpacePercent = 1.2
    pres@gsnPanelYWhiteSpacePercent = 3.
    pres@gsnPanelLabelBar    = False
    pres@gsnPanelRowSpec     = True
    pres@pmLabelBarWidthF    = 0.6
    pres@pmLabelBarHeightF   = 0.07
    pres@lbLabelFontHeightF  = 0.013
    ; pres@lbTitleString      ="units:  m^2*PVU/s"
    gsn_panel(wks,gplt,(/2,1/),pres)
end
