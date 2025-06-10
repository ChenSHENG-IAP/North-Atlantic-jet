load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
siglvl_a = 0.05
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = -30.
lat2 = 90.
lon1 = 0.
lon2 = 360.
lev  = (/250/)
LEV  =  lev * 100

TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/ua_Amon_MPI-ESM_historical_rcp85_ensmean_185001-209912.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ua_mme     = f->ua(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    ua_MME      = month_to_season(ua_mme, "JJA")
    U_clim  = dim_avg_n_Wrap(ua_MME, 0)

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/zg/zg_Amon_MPI-ESM_historical_rcp85_SSTGmin15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    zg_ori     = f->zg(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    zg_min15      = month_to_season(zg_ori, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/ua_Amon_MPI-ESM_historical_rcp85_SSTGmin15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ua_ori     = f->ua(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    ua_min15      = month_to_season(ua_ori, "JJA")
    ; U_clim  = dim_avg_n_Wrap(ua_min15, 0)
    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/va/va_Amon_MPI-ESM_historical_rcp85_SSTGmin15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    va_ori     = f->va(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    va_min15      = month_to_season(va_ori, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ts/ts_Amon_MPI-ESM_historical_rcp85_SSTGmin15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ts_ori     = f->ts(iYYYY,{lat1:lat2},{lon1:lon2})
    ts_min15      = month_to_season(ts_ori, "JJA")
    ; delete([/time,YYYY,iYYYY/])

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/zg/zg_Amon_MPI-ESM_historical_rcp85_SSTGmax15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    zg_ori     = f->zg(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    zg_max15      = month_to_season(zg_ori, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/ua_Amon_MPI-ESM_historical_rcp85_SSTGmax15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ua_ori     = f->ua(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    ua_max15      = month_to_season(ua_ori, "JJA")
    ; U_clim  = dim_avg_n_Wrap(ua_max15, 0)
    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/va/va_Amon_MPI-ESM_historical_rcp85_SSTGmax15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    va_ori     = f->va(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    va_max15      = month_to_season(va_ori, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ts/ts_Amon_MPI-ESM_historical_rcp85_SSTGmax15.ensmean_185001-209912.v3.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ts_ori     = f->ts(iYYYY,{lat1:lat2},{lon1:lon2})
    ts_max15      = month_to_season(ts_ori, "JJA")
    zg_plot =  zg_max15 - zg_min15 
    ua_plot =  ua_max15 - ua_min15
    va_plot =  va_max15 - va_min15
    ts_plot =  ts_max15 - ts_min15 
    
    copy_VarMeta(zg_max15,zg_plot)
    copy_VarMeta(ua_max15,ua_plot)
    copy_VarMeta(va_max15,va_plot)
    copy_VarMeta(va_max15,ts_plot)
    GradLatLon = grad_latlon_cfd(ts_plot, ts_plot&lat, ts_plot&lon, False, False)
    SSTG_plot = GradLatLon[0] * -1.
    copy_VarMeta(zg_max15,SSTG_plot)
    ;;for mask --
    ; assume data is 3D (time,lat,lon)
    a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
    lsdata = a->LSMASK
    lsm = landsea_mask(lsdata,ts_plot&lat,ts_plot&lon)
    ; lsm is a 2D array, in order to use it in mask, we must conform it
    ; to the size of the 3D array "data". 
    ts_plot = mask(ts_plot,conform(ts_plot,lsm,(/1,2/)).ge.1,False)
    SSTG_plot = mask(SSTG_plot,conform(SSTG_plot,lsm,(/1,2/)).ge.1,False)

    meta = zg_max15(0,:,:)
    ; ;**********regression*************************
    x  = TIME*1.
    RC = regCoef_n(x,zg_plot,0,0)
    df = onedtond(RC@nptxy,dimsizes(RC)) - 2 
    tval = onedtond(RC@tval,dimsizes(RC))
    prob_plot = student_t(tval, df)
    zgtrd_plot = RC * 10.
    copy_VarMeta(meta, zgtrd_plot)
    copy_VarMeta(meta, prob_plot)
    ; Trd = dim_rmvmean_n_Wrap(Trd, 2)
    printVarSummary(zgtrd_plot)


    ; ;**********regression*************************
    x  = TIME*1.
    RU = regCoef_n(x,ua_plot,0,0)
    df = onedtond(RU@nptxy,dimsizes(RU)) - 2 
    tval = onedtond(RU@tval,dimsizes(RU))
    prob_U_plot = student_t(tval, df)
    Utrd_plot = RU * 10.
    copy_VarMeta(meta, Utrd_plot)
    copy_VarMeta(meta, prob_U_plot)
    ; Utrd = dim_rmvmean_n_Wrap(Utrd, 2)

    ; ;**********regression*************************
    x  = TIME*1.
    RV = regCoef_n(x,va_plot,0,0)
    df = onedtond(RV@nptxy,dimsizes(RV)) - 2 
    tval = onedtond(RV@tval,dimsizes(RV))
    prob_V_plot = student_t(tval, df)
    Vtrd_plot = RV * 10.
    copy_VarMeta(meta, Vtrd_plot)
    copy_VarMeta(meta, prob_V_plot)
    ; Vtrd = dim_rmvmean_n_Wrap(Vtrd, 2)

    uv_prob_plot = where(prob_U_plot.ge.prob_V_plot,prob_V_plot,prob_U_plot)
    copy_VarMeta(meta, uv_prob_plot)


    ; ;**********regression*************************
    x  = TIME*1.
    RT = regCoef_n(x,ts_plot,0,0)
    df = onedtond(RT@nptxy,dimsizes(RT)) - 2 
    tval = onedtond(RT@tval,dimsizes(RT))
    prob_T_plot = student_t(tval, df)
    Tstrd_plot = RT * 10.
    copy_VarMeta(meta, Tstrd_plot)
    copy_VarMeta(meta, prob_T_plot)
    ; ;**********regression*************************
    x  = TIME*1.
    RSSTG = regCoef_n(x,SSTG_plot,0,0)
    df = onedtond(RSSTG@nptxy,dimsizes(RSSTG)) - 2 
    tval := onedtond(RSSTG@tval,dimsizes(RSSTG))
    prob_SSTG_plot = student_t(tval, df)
    SSTGtrd_plot = RSSTG * 10. * 10^6
    copy_VarMeta(meta, SSTGtrd_plot)
    copy_VarMeta(meta, prob_SSTG_plot)
;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic8/Fig.E7.diff.trend.uv.zg.ts.MPI.NAWH_based.v3."+YrStart+"-"+YrEnd)
  res                               =True
  ;res@gsnPolar                      ="NH"
  font                              =0.03
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

  res@gsnAddCyclic                  = True
  res@mpMinLatF                     =20
  res@mpMaxLatF                     =70
  res@mpMinLonF                     =-90
  res@mpMaxLonF                     =20

;; for contour=============================================================
  res@cnFillPalette                 ="cmp_b2r"   ;;GMT_polar  cmp_b2r   MPL_PRGn
  res@cnFillOn                      =True
  res@cnLinesOn                     =False
  res@cnLineThicknessF              = 1.5
  res@cnMonoLineColor               = False

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
  res@pmLabelBarOrthogonalPosF      = 0.20
  res@pmLabelBarHeightF             = 0.12
  res@lbLabelFontHeightF            = 0.03
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

  resh2@cnLevels                      =  (/13.,16./)
  resh2@cnLinesOn                     = True
  resh2@cnMonoLineColor               = True
  resh2@cnLineColor                   = "black"
  resh2@cnLineThicknessF              = 2.
  ;;for line label==============================================================
  resh2@cnLineLabelsOn                =True
  resh2@cnLineLabelFontHeightF        =0.015  
  resh2@cnLineLabelFontColor          ="black"
  resh2@cnLineLabelBackgroundColor    ="white"
  resh2@cnInfoLabelOn                 = False
  resh2@cnLineLabelPlacementMode         = "Computed"
  resh2@cnLineLabelInterval              =1
  resh2@cnLineLabelDensityF              =0.5
  resh2@cnLabelDrawOrder              ="PostDraw"
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

;;for vector ===================== 

vres                          = True
vres@gsnDraw                  = False
vres@gsnFrame                 = False
vres@gsnMaximize              = False
vres@gsnLeftString            = ""
vres@gsnRightString           = ""
vres@vcMinDistanceF           = 0.03           
;vres@vcMinFracLengthF         = 0.3
vres@vcMinMagnitudeF          = 0.3    

;;LineArrow----------------------------
 vres@vcGlyphStyle             = "CurlyVector"
;  vres@vcMonoLineArrowColor     = "True"         
 vres@vcLineArrowColor         = "black"
 vres@vcLineArrowThicknessF    = 1.             

;;========================================
vres@vcLevelSelectionMode     = "ExplicitLevels"                 
vres@vcLevels                 = (/siglvl_a/)        ;95%
vres@vcLevelColors            = (/"blue","black"/)
vres@lbLabelBarOn             =  False
;; set ref-------------------------------- 

value = 0.6
vres@vcRefAnnoOn             = True
vres@vcRefMagnitudeF         = value            
vres@vcRefAnnoString1        = 0.8
vres@vcRefAnnoString2On      = False
vres@vcRefLengthF            = 0.04          ; define length of vec ref
vres@vcRefAnnoOrthogonalPosF = -0.19
vres@vcRefAnnoParallelPosF   = 0.999
vres@vcRefAnnoSide           = "Bottom"
vres@vcRefAnnoBackgroundColor = "grey"
vres@vcRefAnnoFontHeightF    = 0.015
vres@vcRefAnnoPerimOn        = False

  max_vaule = 1.
  interval  = 0.2
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
  res@gsnLeftString   = ""
  plot  = new(3, graphic)
  dplot = new(3, graphic)
  vplot = new(3, graphic)
  dum = new(3, graphic)

  res@gsnLeftString = "";;"(a)U250 trend"
  res@gsnRightString = "";;"max15-min15"
  plot(0) = gsn_csm_contour_map(wks, Utrd_plot, res)
  dplot(0)= gsn_csm_contour(wks, prob_U_plot, opt)
  vplot(0)  = gsn_csm_vector_scalar(wks, Utrd_plot, Vtrd_plot, uv_prob_plot, vres)



ncols = 14
  midcol2 = toint(ncols/2.) 
  midcol1 = midcol2 - 1

  colors := read_colormap_file(res@cnFillPalette)
  icol   := span_color_indexes(colors(2:,:), ncols)  
  icol = icol+2  
  icol(midcol1:midcol2) = -1
  res@cnLevelSelectionMode = "ExplicitLevels"   
  res@cnLevels             = (/-0.5,-0.4,-0.3,-0.2,-0.1,-0.02,0,0.02,0.1,0.2,0.3,0.4,0.5/)

  res@cnFillColors := icol
  res@gsnLeftString = "";;"(b)SST trend"
  plot(1) = gsn_csm_contour_map(wks, Tstrd_plot, res)
  dplot(1)= gsn_csm_contour(wks, prob_T_plot, opt)


  max_vaule = 0.6
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
  plot(2) = gsn_csm_contour_map(wks, SSTGtrd_plot, res)
  dplot(2)= gsn_csm_contour(wks, prob_SSTG_plot, opt)

do i = 0, 0
  overlay(plot(i), vplot(i))
end do

do i = 1, 1
  overlay(plot(i), dplot(i))
end do

do i = 2, 2
  overlay(plot(i), dplot(i))
end do
  h2plot = gsn_csm_contour(wks,U_clim,resh2)
  overlay(plot(0),h2plot)

    

  pres                     =True
  pres@gsnPanelDebug      = True                   
  ; pres@gsnPanelXF         = (/ 0.0556078, 0.575608/)
  pres@txString            =""
  pres@gsnPanelLabelBar    = False
  pres@gsnPanelRowSpec     = False
  pres@pmLabelBarWidthF    = 0.6
  pres@pmLabelBarHeightF   = 0.07
  pres@lbLabelFontHeightF  = 0.013
  ; ; pres@lbTitleString      ="units:  m^2*PVU/s"
  gsn_panel(wks,plot,(/1,3/),pres)
  wallClockElapseTime(TotalTime, "totaltime", 0)  
end