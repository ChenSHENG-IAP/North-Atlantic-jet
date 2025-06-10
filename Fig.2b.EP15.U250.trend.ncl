load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
siglvl_a = 0.05 
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 0.
lat2 = 90.
lon1 = 0.
lon2 = 360.
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)


dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/"
fils  = systemfunc("ls -v "+dir+"ua_Amon_MPI-ESM_historical_rcp85_ensmean_185001-209912.nc")
fP    = addfiles(fils,"r") 
time   := fP[:]->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
U_ensmean     = fP[:]->ua(iYYYY,{25000},:,:)
U_ensmean_JJA = month_to_season(U_ensmean, "JJA")
U_clim  = dim_avg_n_Wrap(U_ensmean_JJA, 0)
print("111")
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/"
fils  = systemfunc("ls -v "+dir+"ua_Amon_MPI-ESM_historical_rcp85_P15.ensmean_185001-209912.v2.nc")
fP    = addfiles(fils,"r") 
time   := fP[:]->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
UP     = fP[:]->ua(iYYYY,{25000},:,:)
UP_JJA = month_to_season(UP, "JJA")
print("222")

fils  = systemfunc("ls -v "+dir+"ua_Amon_MPI-ESM_historical_rcp85_E15.ensmean_185001-209912.v2.nc")
fE    = addfiles(fils,"r") 
time   := fE[:]->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
UE     = fE[:]->ua(iYYYY,{25000},:,:)
UE_JJA = month_to_season(UE, "JJA")
meta   = UE(0,:,:)
print("333")

Utrd_E15 = regCoef_n(TIME,UE_JJA-UP_JJA,0,0) 
df_E15 := onedtond(Utrd_E15@nptxy,dimsizes(Utrd_E15)) - 2 
tval_E15 := onedtond(Utrd_E15@tval,dimsizes(Utrd_E15))
Uprob_E15 = student_t(tval_E15, df_E15)



Utrd_E15 = Utrd_E15 * 10.
Utrd_E15 = smth9(Utrd_E15, 0.50, 0.25, False)
copy_VarMeta(meta, Utrd_E15)
copy_VarMeta(meta, Uprob_E15)

 
;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic2/Fig.2.EP15.U250.trend")
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
    res@mpMinLatF                     =30
    res@mpMaxLatF                     =70
    res@mpMinLonF                     =-90
    res@mpMaxLonF                     =20

    ;; for contour===============================================================
    res@cnFillPalette                 ="cmp_b2r"   ;;GMT_polar  cmp_b2r   cmp_b2r
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
    ;;for color bar==============================================================
    res@lbLabelBarOn                  = True
    res@pmLabelBarOrthogonalPosF      = 0.23
    ;res@pmLabelBarHeightF             = 0.07
    res@lbLabelFontHeightF            = 0.025
    ;;==========================specify the contour==============================
    resh1 = True
    resh1@gsnDraw                       = False
    resh1@gsnFrame                      = False
    resh1@gsnMaximize                   = False
    resh1@cnLevelSelectionMode          = "ExplicitLevels"
    resh1@gsnLeftString                 = ""
    resh1@gsnRightString                = ""
    resh1@cnFillOn                      = False
    resh1@cnFillPalette                 = "MPL_PRGn"  

    resh1@cnLevels                      =  (/15.,20.,25.,30,40,50,60,70./)
    resh1@cnLinesOn                     = True
    resh1@cnMonoLineColor               = True
    resh1@cnLineColor                   = "black"
    resh1@cnLineThicknessF              = 1.5
    ;;for line label==============================================================
    resh1@cnLineLabelsOn                =True
    resh1@cnLineLabelFontHeightF        =0.012  
    resh1@cnLineLabelFontColor          ="black"
    resh1@cnLineLabelBackgroundColor    ="gray"
    resh1@cnInfoLabelOn                 = False
    resh1@cnLineLabelPlacementMode         = "Computed"
    resh1@cnLineLabelInterval              =1
    resh1@cnLineLabelDensityF              =0.5   
    resh1@cnLabelDrawOrder              = "PostDraw"

    ;;for confidence lev dot=======================================================
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
    max_vaule = 1.5
    interval  = 0.3
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
    res@gsnRightString  = ""
    
    plot   = new(2, graphic)
    hplot2 = new(2, graphic)
    dplot  = new(2, graphic)



    res@gsnRightString  = "";;"U250(m/s 10a~S~-1~N~)"
    res@gsnLeftString   = "";;"(a)E15-P15"
    plot(0)   = gsn_csm_contour_map(wks, Utrd_E15, res)  
    dplot(0)  = gsn_csm_contour(wks, Uprob_E15, opt)
    hplot2(0) = gsn_csm_contour(wks, U_clim, resh1)


    plot(1)   = gsn_csm_contour_map(wks, Utrd_E15, res)  
    dplot(1)  = gsn_csm_contour(wks, Uprob_E15, opt)
    hplot2(1) = gsn_csm_contour(wks, U_clim, resh1)
    do i = 0, 1
      overlay(plot(i), dplot(i))
      overlay(plot(i), hplot2(i))
    end do
    pres                     =True
    pres@txString            =""
    pres@gsnPanelXWhiteSpacePercent = 1.2
    pres@gsnPanelYWhiteSpacePercent = 1.5
    pres@gsnPanelLabelBar    = False
    pres@gsnPanelRowSpec     = False
    pres@pmLabelBarWidthF    = 0.6
    pres@pmLabelBarHeightF   = 0.07
    pres@lbLabelFontHeightF  = 0.018
    gsn_panel(wks,plot,(/1,2/),pres)
end
