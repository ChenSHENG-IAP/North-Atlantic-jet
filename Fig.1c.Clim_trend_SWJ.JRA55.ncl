load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
siglvl_a = 0.05
;;==========================
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 0.
lat2 = 90.
lon1 = 0.
lon2 = 360.
LEV  = 250
lev1 = 1000
lev2 = 70  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)
dir   = "/data6/shengc2017/data/JRA55/"
f     = addfile(dir+"uwnd.jra55.mon.195801-202212.nc","r")
time   = f->time
YYYY   = cd_calendar(time,-1)/100    ; entire file
iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)

UO     = f->uwnd(iYYYY,:,{lat1:lat2},{lon1:lon2})
U     = month_to_season(UO, "JJA")
U     =lonFlip(U)
U200  = U(:,{LEV},:,:)   
U200_clim = dim_avg_n_Wrap(U200, 0)    

U_zonalmean = dim_avg_n_Wrap(U(:,{lev1:lev2},:,{-45:0}),3)
U_zonalmean_clim = dim_avg_n_Wrap(U_zonalmean,0)

x  = TIME-1980
RC = regCoef_n(x,U200,0,0)
df = onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval = onedtond(RC@tval,dimsizes(RC))
prob = student_t(tval, df)
Trd = RC * 10.
copy_VarMeta(U200(0,:,:), Trd)
copy_VarMeta(U200(0,:,:), prob)


x  = TIME-1980
RC2 = regCoef_n(x,U_zonalmean,0,0)
df2 = onedtond(RC2@nptxy,dimsizes(RC2)) - 2 
tval2 = onedtond(RC2@tval,dimsizes(RC2))
prob2 = student_t(tval2, df2)
Trd2 = RC2 * 10.
copy_VarMeta(U_zonalmean(0,:,:), Trd2)
copy_VarMeta(U_zonalmean(0,:,:), prob2)

;;======================
wks   = gsn_open_wks ("eps", "/data6/shengc2017/work9_sub_Jet/NCC/pic1/Fig.1.Clim_trend_NAJ_JRA55")
  res                               =True
  ;res@gsnPolar                      ="NH"
  font                              =0.02
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

  res@tmXBTickSpacingF              = 10

  res@gsnAddCyclic                  = False
  res@mpMinLatF                     =35
  res@mpMaxLatF                     =75
  res@mpMinLonF                     =-70
  res@mpMaxLonF                     =5
  ; res@mpCenterLonF                  = 0 
; res@gsnPolarLabelFontHeightF      = font
  ;res@mpLandFillColor            = "white"                   
  res@mpShapeMode                   ="FreeAspect"
  res@vpWidthF                      = 0.8
  res@vpHeightF                     = 0.4
;; for contour=============================================================
  res@cnFillPalette                 ="cmp_b2r"   ;;GMT_polar  cmp_b2r   MPL_PRGn
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
  res@pmLabelBarOrthogonalPosF      = 0.17
  res@pmLabelBarHeightF             = 0.07
  res@lbLabelFontHeightF            = 0.019
;;;============================specify the contour==============================================================
  resh                               = True
  resh@gsnDraw                       = False
  resh@gsnFrame                      = False
  resh@gsnMaximize                   = False
  resh@cnLevelSelectionMode          = "ExplicitLevels"
  resh@cnLevels                      = (/3000/)
  resh@cnLineLabelsOn                = False
  resh@cnMonoLineThickness           = True
  resh@cnLineThicknessF              = 2 
  resh@cnMonoLineColor               = True
  resh@cnLineColor                   = "navyblue"
  resh@cnLineDashPattern             = 0 
  resh@gsnLeftString                 = ""
  resh@gsnRightString                = ""
  resh@cnInfoLabelOn                 = False   
;;====================================================================================================
;;;==zonal wind jet==========================specify the contour==============================================================
    resh2 = True
    resh2@gsnDraw                       = False
    resh2@gsnFrame                      = False
    resh2@gsnMaximize                   = False
    resh2@cnLevelSelectionMode          = "ExplicitLevels"
    resh2@gsnLeftString                 = ""
    resh2@gsnRightString                = ""
    resh2@cnFillOn                      = False
    resh2@cnFillPalette                 = "cmp_b2r"  

    resh2@cnLevels                      = (/8.,16.,24.,30,40,50,60,70./)
    resh2@cnLinesOn                     = True
    resh2@cnMonoLineColor               = True
    resh2@cnLineColor                   = "black"
    resh2@cnLineThicknessF              = 2.
;;for line label==============================================================
    resh2@cnLineLabelsOn                =True
    resh2@cnLineLabelFontHeightF        =0.012  
    resh2@cnLineLabelFontColor          ="black"
    resh2@cnLineLabelBackgroundColor    ="white"
    resh2@cnInfoLabelOn                 = False
    resh2@cnLineLabelPlacementMode         = "Computed"
    resh2@cnLineLabelInterval              =1
    resh2@cnLineLabelDensityF              =0.5

;;=================================================
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
    max_vaule = 2.
    interval  = 0.25
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
;;=============================================
  res@gsnRightString  = "";;"JJA(0.05)"
  res@gsnLeftString   = "";;"(a)U@"+LEV+"hPa"
    plot   = gsn_csm_contour_map(wks, Trd, res)  
    hhplot = gsn_csm_contour(wks, U200_clim, resh2)
    dplot  = gsn_csm_contour(wks, prob, opt)
    overlay(plot, dplot)
    overlay(plot, hhplot)

;;for pres-lat ----------------------------------
delete(res)
  res                               =True
  ; res@gsnPolar                      ="NH"
  font                              =0.02
  res@gsnMaximize                   =False
  res@gsnDraw                       =False
  res@gsnFrame                      =False
  res@tmXTOn                        =False
  res@tmYROn                        =False 
  ovalue = 0.

  res@gsnRightString                = ""
  res@gsnLeftStringOrthogonalPosF   = ovalue
  res@gsnRightStringOrthogonalPosF  = ovalue
  res@gsnLeftStringFontHeightF      = font
  res@gsnRightStringFontHeightF     = font
  res@tmYLLabelFontHeightF          = font
  res@tmXBLabelFontHeightF          = font
  res@tmXBTickSpacingF              = 10
  res@tiYAxisString                 = ""
;; for contour=============================================================
  res@cnFillPalette                 ="cmp_b2r"    ;;BlueDarkRed18   cmp_b2r  GMT_polar   MPL_PRGn
  res@cnFillOn                      =True
  res@cnLinesOn                     =False
  res@cnLineThicknessF              = 0.0000001 
  res@cnLineColor                   = "blue"
  res@gsnContourZeroLineThicknessF  = 3.
;;for line label==============================================================
  res@cnLineLabelsOn                =False
  res@cnLineLabelFontHeightF        =0.02  
  res@cnLineLabelFontColor          ="black"
  res@cnLineLabelBackgroundColor    ="white"
  res@cnInfoLabelOn                 = False
  res@cnLineLabelPlacementMode         = "Computed"
  res@cnLineLabelInterval              =2
  res@cnLineLabelDensityF              =0.5


;;for color bar===========================================================
  res@lbLabelBarOn                  = True
  res@pmLabelBarOrthogonalPosF      = 0.04
  res@pmLabelBarHeightF             = 0.07
  res@lbLabelFontHeightF            = 0.019
;;============================================
    max_vaule = 1.4
    interval  = 0.2
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors = read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols) 
    icol = icol+2  
    icol(midcol1:midcol2) = -1

    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors = icol
  ; res@mpShapeMode                   ="FreeAspect"
  res@vpWidthF                      = 0.5
  res@vpHeightF                     = 0.4  
  res@gsnRightString  = "";;"JJA(0.05)"
  res@gsnLeftString   = "";;"(a)U@"+LEV+"hPa"
  res@gsnPresHgtHeightLabelOn=False
  plot2 = gsn_csm_pres_hgt(wks, Trd2(:,{35:75}), res)

  opt@gsnPresHgtHeightLabelOn=False
  dplot2 = gsn_csm_pres_hgt(wks, prob2, opt)

  resh2@gsnPresHgtHeightLabelOn=False
  resh2@cnLevels  := (/10.,15.,20.,24.,30.,35.,40.,45.,50./)
  resh2@cnLineLabelFontHeightF        =0.018 
  hhplot2 = gsn_csm_pres_hgt(wks, U_zonalmean_clim(:,{35:75}), resh2)
  overlay(plot2,dplot2)
  overlay(plot2,hhplot2)

  gplt = new(2, graphic)
  gplt(0) = plot 
  gplt(1) = plot2
  pres                     =True
  pres@gsnPanelDebug      = True                   
  pres@gsnPanelXF         = (/ 0.0556078, 0.575608/)
  pres@txString            =""
  gsn_panel(wks,gplt,(/1,2/),pres)
  wallClockElapseTime(TotalTime, "totaltime", 0)  
end