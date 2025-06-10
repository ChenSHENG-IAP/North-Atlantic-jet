load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
siglvl_a = 5.  ;;making dot when more than 5 models have the same sign.
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2014
lat1 = 0.
lat2 = 90.
lon1 = 0.
lon2 = 360.
LEV  = 25000
lev1 = 100000
lev2 = 7000  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)
model = (/"BCC-CSM2-MR","CanESM5","CNRM-CM6-1","HadGEM3-GC31-LL","IPSL-CM6A-LR","MRI-ESM2-0"/)
expt="historical"  ;;"hist-GHG"  "hist-nat" "hist-aer" "historical"
expt_out = expt
if (expt .eq. "historical") then
  expt_out = "All"
end if
dir   = "/gws/nopw/j04/ncas_climate_vol1/users/shengc/DAMIP/multi-model-2x2/ua/"+expt+"/"
  do imodel = 0,5 
    fils  = systemfunc("ls "+dir+"ua_Amon_"+model(imodel)+"_"+expt+"_ensmean_g?_????01-????12.2x2.nc")
    fs     = addfiles(fils,"r")
    time   := fs[:]->time
    YYYY   := cd_calendar(time,-1)/100    ; entire file
    iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    UO     = fs[:]->ua(iYYYY,:,{lat1:lat2},{lon1:lon2})
    U      = month_to_season(UO, "JJA")
    U      =lonFlip(U)
    U200  = U(:,{LEV},:,:)    ;;U(:,{LEV},:,:)     dim_avg_n_Wrap(U(:,{lev1:lev2},:,:), 1)
    U200_clim = dim_avg_n_Wrap(U200, 0)    
    U_zonalmean = dim_avg_n_Wrap(U(:,{lev1:lev2},:,{-45:0}),3)
    U_zonalmean_clim = dim_avg_n_Wrap(U_zonalmean,0)
    ; ;**********regression map*************************
    x  = TIME-1980
    RC = regCoef_n(x,U200,0,0)
    df = onedtond(RC@nptxy,dimsizes(RC)) - 2 
    tval = onedtond(RC@tval,dimsizes(RC))
    prob = student_t(tval, df)
    Trd = RC * 10.
    copy_VarMeta(U200(0,:,:), Trd)
    copy_VarMeta(U200(0,:,:), prob)

    ; ;**********regression cross-section*************************
    x  = TIME-1980
    RC2 = regCoef_n(x,U_zonalmean,0,0)
    df2 = onedtond(RC2@nptxy,dimsizes(RC2)) - 2 
    tval2 = onedtond(RC2@tval,dimsizes(RC2))
    prob2 = student_t(tval2, df2)
    Trd2 = RC2 * 10.
    copy_VarMeta(U_zonalmean(0,:,:), Trd2)
    copy_VarMeta(U_zonalmean(0,:,:), prob2)
    if ( imodel .eq. 0) then
      U200_clim_array           = U200(0:5,:,:)
      Trd_array                 = U200(0:5,:,:)
      U_zonalmean_clim_array    = U_zonalmean(0:5,:,:)
      Trd2_array                = U_zonalmean(0:5,:,:)
    end if
    U200_clim_array(imodel,:,:) = U200_clim
    Trd_array(imodel,:,:)       = Trd

    U_zonalmean_clim_array(imodel,:,:)=U_zonalmean_clim
    Trd2_array(imodel,:,:)      = Trd2
  end do 

    U200_clim_MME = dim_avg_n_Wrap(U200_clim_array,0)
    Trd_MME       = dim_avg_n_Wrap(Trd_array,0)
    U_zonalmean_clim_MME = dim_avg_n_Wrap(U_zonalmean_clim_array,0)
    Trd2_MME      = dim_avg_n_Wrap(Trd2_array,0)

    sig_flg = Trd_array
    sig2_flg = Trd2_array
    do imodel = 0, 5
      sig_flg(imodel,:,:) = where((Trd_array(imodel,:,:)*Trd_MME) .ge. 0., 1., 0.)
      sig2_flg(imodel,:,:) = where((Trd2_array(imodel,:,:)*Trd2_MME) .ge. 0., 1., 0.)
    end do
    SameSign_num  = dim_sum_n_Wrap(sig_flg, 0)
    SameSign2_num = dim_sum_n_Wrap(sig2_flg,0)

;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic1/Fig.E2.Clim_trend_NAJ_MME_"+expt_out)
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

    res@tmXBTickSpacingF              = 60

    res@gsnAddCyclic                  = True
    res@mpMinLatF                     =lat1
    res@mpMaxLatF                     =lat2
    ; res@mpMinLonF                     =-70
    ; res@mpMaxLonF                     =5
    ; res@mpCenterLonF                  = -90. 
    ; res@gsnPolarLabelFontHeightF      = font
    res@mpLandFillColor            = "white"                   
    res@mpShapeMode                   ="FreeAspect"
    res@vpWidthF                      = 0.8
    res@vpHeightF                     = 0.4
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
    res@lbLabelBarOn                  = False
    res@pmLabelBarOrthogonalPosF      = 0.17
    res@pmLabelBarHeightF             = 0.07
    res@lbLabelFontHeightF            = 0.019
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
    ;;====================================================================================================
    ;;;==zonal wind jet==========================specify the contour======================================
    resh2 = True
    resh2@gsnDraw                       = False
    resh2@gsnFrame                      = False
    resh2@gsnMaximize                   = False
    resh2@cnLevelSelectionMode          = "ExplicitLevels"
    resh2@gsnLeftString                 = ""
    resh2@gsnRightString                = ""
    resh2@cnFillOn                      = False
    resh2@cnFillPalette                 = "cmp_b2r"  

    resh2@cnLevels                      =  (/8.,16.,24.,30,40,50,60,70./)
    resh2@cnLinesOn                     = True
    resh2@cnMonoLineColor               = True
    resh2@cnLineColor                   = "black"
    resh2@cnLineThicknessF              = 1.5
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
    opt@cnFillColors           = (/"Transparent" ,"white"/)
    opt@cnFillDotSizeF         =  0.003
    opt@cnInfoLabelOn          = False
    ;;============================================
    max_vaule = 0.5
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
    res@gsnRightString  = "";;"JJA"
    res@gsnLeftString   = "";;"(a)trd of U250 ALL(MME)"
    plot   = gsn_csm_contour_map(wks, Trd_MME, res)  
    hhplot = gsn_csm_contour(wks, U200_clim_MME, resh2)
    dplot  = gsn_csm_contour(wks, SameSign_num, opt)
    overlay(plot, dplot)
    overlay(plot, hhplot)

    ;;for pres-lat ----------------------------------
    delete(res)
    res                               =True
    ; res@gsnPolar                      ="NH"
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
    res@cnFillPalette                 ="cmp_b2r"    ;;BlueDarkRed18   cmp_b2r  GMT_polar
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
    res@lbLabelBarOn                  = False
    res@pmLabelBarOrthogonalPosF      = 0.04
    res@pmLabelBarHeightF             = 0.07
    res@lbLabelFontHeightF            = 0.019
    ;;============================================
    res@tmYLMinorOn     = True
    res@tmYLMinorValues  = Trd2_MME&plev/100.
    res@tmYLMode         = "Explicit"
    res@tmYLValues        = (/1000,700,500,250,150,70/)
    res@tmYLLabels        = (/1000,700,500,250,150,70/)*1.
    res@tmXBTickSpacingF              = 20
    max_vaule = 0.5
    interval  = 0.1
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
    res@gsnRightString  = "";;"JJA"
    res@gsnLeftString   = "";;"(b)NAJ trend"
    res@gsnPresHgtHeightLabelOn=False
    plot2 = gsn_csm_pres_hgt(wks, Trd2_MME(:,{35:75}), res)

    opt@gsnPresHgtHeightLabelOn=False
    dplot2 = gsn_csm_pres_hgt(wks, SameSign2_num, opt)

    resh2@gsnPresHgtHeightLabelOn=False
    resh2@cnLevels  := (/10.,15.,20.,24.,30.,35.,40.,45.,50./)
    resh2@cnLineLabelFontHeightF        =0.018 
    hhplot2 = gsn_csm_pres_hgt(wks, U_zonalmean_clim_MME(:,{35:75}), resh2)
    overlay(plot2,dplot2)
    overlay(plot2,hhplot2)


  min_lat = 30.
  max_lat = 70.
  min_lon = 315.   ;; 71.
  max_lon = 360.  ;; 83.
  boxlat  = (/max_lat,max_lat,min_lat,min_lat,max_lat/)
  boxlon  = (/min_lon,max_lon,max_lon,min_lon,min_lon/) 

  lnres    = True
  lnres@gsLineColor = "black"
  lnres@gsLineThicknessF = 1.0
  lnres@gsLineDashPattern = 0
  dum = gsn_add_polyline(wks, plot, boxlon, boxlat, lnres) 

    gplt = new(2, graphic)
    gplt(0) = plot 
    gplt(1) = plot2
    pres                     =True
    pres@gsnPanelDebug      = True                  
    pres@gsnPanelXF         = (/ 0.0556078, 0.575608/)+0.1
    pres@txString            =""
    pres@gsnPanelLabelBar    = True
    pres@lbLabelFontHeightF  = 0.02

    gsn_panel(wks,gplt,(/1,2/),pres)
    wallClockElapseTime(TotalTime, "totaltime", 0)  
end