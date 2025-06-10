load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
siglvl_a = 0.1
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2014
lat1 = -30.
lat2 = 90.
lon1 = 0.
lon2 = 360.
lev  = (/250/)
LEV  =  lev * 100

TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/historical-ensmean-10mems/zg_Amon_IPSL-CM6A-LR_historical_ensmean_gr_185001-201412.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    zg_mme     = f->zg(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    zg_MME      = month_to_season(zg_mme, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/historical-ensmean-10mems/ua_Amon_IPSL-CM6A-LR_historical_ensmean_gr_185001-201412.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ua_mme     = f->ua(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    ua_MME      = month_to_season(ua_mme, "JJA")
    U_clim  = dim_avg_n_Wrap(ua_MME, 0)
    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/historical-ensmean-10mems/va_Amon_IPSL-CM6A-LR_historical_ensmean_gr_185001-201412.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    va_mme     = f->va(iYYYY,{LEV},{lat1:lat2},{lon1:lon2})
    va_MME      = month_to_season(va_mme, "JJA")

    dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/historical-ensmean-10mems/ts_Amon_IPSL-CM6A-LR_historical_ensmean_gr_185001-201412.nc"
    f     = addfile(dir,"r")
    time   := f->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    ts_mme     = f->ts(iYYYY,{lat1:lat2},{lon1:lon2})
    ts_MME      = month_to_season(ts_mme, "JJA")

    delete([/time,YYYY,iYYYY/])

    expt  = "atl"  ;; atl  pac
    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/"+expt+"-pacemaker/"
    f   = addfile(dir+"zg_Amon_IPSL-CM6A-LR_dcppC-"+expt+"-pacemaker_s1920-ensmean_gr_193001-201412.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    hgt  = f->zg(iYYYY,{LEV},{lat1:lat2},{lon1:lon2}) 
    zg_atl = month_to_season(hgt, "JJA")


    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/"+expt+"-pacemaker/"
    f   = addfile(dir+"ua_Amon_IPSL-CM6A-LR_dcppC-"+expt+"-pacemaker_s1920-ensmean_gr_193001-201412.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    U  = f->ua(iYYYY,{LEV},{lat1:lat2},{lon1:lon2}) 
    ua_atl = month_to_season(U, "JJA")

    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/"+expt+"-pacemaker/"
    f   = addfile(dir+"va_Amon_IPSL-CM6A-LR_dcppC-"+expt+"-pacemaker_s1920-ensmean_gr_193001-201412.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    V  = f->va(iYYYY,{LEV},{lat1:lat2},{lon1:lon2}) 
    va_atl = month_to_season(V, "JJA")

    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/pacemaker/DCPP/IPSL-CM6A-LR/"+expt+"-pacemaker/"
    f   = addfile(dir+"ts_Amon_IPSL-CM6A-LR_dcppC-"+expt+"-pacemaker_s1920-ensmean_gr_193001-201412.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    TS  = f->ts(iYYYY,{lat1:lat2},{lon1:lon2}) 
    ts_atl = month_to_season(TS, "JJA")

    zg_atl = zg_atl - zg_MME 
    ua_atl = ua_atl - ua_MME
    va_atl = va_atl - va_MME
    ts_atl = ts_atl - ts_MME 

    copy_VarMeta(zg_MME,zg_atl)
    copy_VarMeta(ua_MME,ua_atl)
    copy_VarMeta(va_MME,va_atl)
    copy_VarMeta(va_MME,ts_atl)
    ; SSTG_atl    = -1. * center_finite_diff_n(ts_atl, ts_atl&lat, False, 0, 1)
    GradLatLon = grad_latlon_cfd(ts_atl, ts_atl&lat, ts_atl&lon, False, False)
    SSTG_atl = GradLatLon[0] * -1.
    copy_VarMeta(zg_atl,SSTG_atl)
;;for mask --
; assume data is 3D (time,lat,lon)
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata,ts_atl&lat,ts_atl&lon)
; lsm is a 2D array, in order to use it in mask, we must conform it
; to the size of the 3D array "data". 
ts_atl = mask(ts_atl,conform(ts_atl,lsm,(/1,2/)).ge.1,False)
SSTG_atl = mask(SSTG_atl,conform(SSTG_atl,lsm,(/1,2/)).ge.1,False)    
meta = zg_MME(0,:,:)
    ; ;**********regression*************************
    x  = TIME*1.
    RC = regCoef_n(x,zg_atl,0,0)
    df = onedtond(RC@nptxy,dimsizes(RC)) - 2 
    tval = onedtond(RC@tval,dimsizes(RC))
    prob_atl = student_t(tval, df)
    Trd_atl = RC * 10.
    copy_VarMeta(meta, Trd_atl)
    copy_VarMeta(meta, prob_atl)
    ; Trd = dim_rmvmean_n_Wrap(Trd, 2)
    printVarSummary(Trd_atl)


    ; ;**********regression*************************
    x  = TIME*1.
    RU = regCoef_n(x,ua_atl,0,0)
    df = onedtond(RU@nptxy,dimsizes(RU)) - 2 
    tval = onedtond(RU@tval,dimsizes(RU))
    prob_U_atl = student_t(tval, df)
    Utrd_atl = RU * 10.
    copy_VarMeta(meta, Utrd_atl)
    copy_VarMeta(meta, prob_U_atl)
    ; Utrd = dim_rmvmean_n_Wrap(Utrd, 2)

    ; ;**********regression*************************
    x  = TIME*1.
    RV = regCoef_n(x,va_atl,0,0)
    df = onedtond(RV@nptxy,dimsizes(RV)) - 2 
    tval = onedtond(RV@tval,dimsizes(RV))
    prob_V_atl = student_t(tval, df)
    Vtrd_atl = RV * 10.
    copy_VarMeta(meta, Vtrd_atl)
    copy_VarMeta(meta, prob_V_atl)
    ; Vtrd = dim_rmvmean_n_Wrap(Vtrd, 2)


             
    uv_prob_atl = where(prob_U_atl.ge.prob_V_atl,prob_V_atl,prob_U_atl)
    copy_VarMeta(meta, uv_prob_atl)


    ; ;**********regression*************************
    x  = TIME*1.
    RT = regCoef_n(x,ts_atl,0,0)
    df = onedtond(RT@nptxy,dimsizes(RT)) - 2 
    tval = onedtond(RT@tval,dimsizes(RT))
    prob_T_atl = student_t(tval, df)
    Tstrd_atl = RT * 10.
    copy_VarMeta(meta, Tstrd_atl)
    copy_VarMeta(meta, prob_T_atl)

    ; ;**********regression*************************
    x  = TIME*1.
    RSSTG = regCoef_n(x,SSTG_atl,0,0)
    df = onedtond(RSSTG@nptxy,dimsizes(RSSTG)) - 2 
    tval := onedtond(RSSTG@tval,dimsizes(RSSTG))
    prob_SSTG_atl = student_t(tval, df)
    SSTGtrd_atl = RSSTG * 10. * 10^6
    copy_VarMeta(meta, SSTGtrd_atl)
    copy_VarMeta(meta, prob_SSTG_atl)
;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic8/Fig.E6.fgh.IPSL.atl.pacemaker.uv.v3."+YrStart+"-"+YrEnd)
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
;;;============================specify the contour========================
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
;;;==zonal wind jet==========================specify the contour==========
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
;;=========================================
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
  interval  = 0.1
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
  hplot = new(3, graphic)
  dplot = new(3, graphic)
  vplot = new(3, graphic)
  dum = new(3, graphic)

  res@gsnLeftString = "";;"(a)U250 trend"
  res@gsnRightString = "";;"atl-pacemaker"
  plot(0) = gsn_csm_contour_map(wks, Utrd_atl, res)
  dplot(0)= gsn_csm_contour(wks, prob_U_atl, opt)
  vplot(0)  = gsn_csm_vector_scalar(wks, Utrd_atl, Vtrd_atl, uv_prob_atl, vres)


  max_vaule = 0.5
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
  res@gsnLeftString = "";;"(b)SST trend"

  plot(1) = gsn_csm_contour_map(wks, Tstrd_atl, res)
  dplot(1)= gsn_csm_contour(wks, prob_T_atl, opt)

  max_vaule = 0.8
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
  res@gsnLeftString = "";;"(c)SSTG trend"
  plot(2) = gsn_csm_contour_map(wks, SSTGtrd_atl, res)
  dplot(2)= gsn_csm_contour(wks, prob_SSTG_atl, opt)
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
  pres@txString            =""
  pres@gsnPanelLabelBar    = False
  pres@gsnPanelRowSpec     = False
  pres@pmLabelBarWidthF    = 0.6
  pres@pmLabelBarHeightF   = 0.07
  pres@lbLabelFontHeightF  = 0.013
  gsn_panel(wks,plot,(/1,3/),pres)
  wallClockElapseTime(TotalTime, "totaltime", 0)  
end