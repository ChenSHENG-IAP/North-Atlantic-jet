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
llev = (/250,500/)
LEV  = llev*100.
TIME = ispan(YrStart, YrEnd, 1) * 1.

dir_in = "/home/users/shengc01/work11/scripts10/data/"
fil_in = "ERA5.trend.Z.1980.2019.FromQOMOtoJasmin.WORK9.NCC.scripts1.nc"
fin = addfile(dir_in+fil_in, "r")
Trd_Z_obs = fin->Trd
prob_Z_obs= fin->prob 

dir   = "/home/users/shengc01/data/SST/"
f     = addfile(dir + "sst.mnmean.COBE.nc", "r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
SSTO   := f->sst(iYYYY,{lat1:lat2},:)
SST_COBE    = month_to_season(SSTO, "JJA")
meta = SST_COBE

GMSST   := f->sst(iYYYY,{-60:60},:)
GMSST_JJA   = month_to_season(GMSST, "JJA")

rad    = 4.0*atan(1.0)/180.0
clat  := cos(GMSST_JJA&lat * rad)
SST_glb_ts = wgt_areaave_Wrap(GMSST_JJA, clat, 1.0, 0)
SST_glb_ts_C := conform(SST_COBE, SST_glb_ts, 0)
SST_COBE   = SST_COBE - SST_glb_ts_C
copy_VarMeta(meta,SST_COBE)

; ;**********regression*************************
x  = (TIME-1980)*1.
RC := regCoef_n(x,SST_COBE,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_SST_COBE := student_t(tval, df)
SST_trd_COBE = RC * 10. 
copy_VarMeta(SST_COBE(0,:,:), SST_trd_COBE)
copy_VarMeta(SST_COBE(0,:,:), prob_SST_COBE)
SST_trd_COBE = where(ismissing(SST_COBE(0,:,:)), SST_COBE@_FillValue, SST_trd_COBE)
print("OBS...")
delete([/RC,df,tval/])

;;for MPI MME-----------------------------
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/zg/zg_Amon_MPI-ESM_historical_rcp85_ensmean_185001-209912.nc"
f     = addfile(dir,"r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
zg_MME     = f->zg(iYYYY,{LEV},:,:)
zg_MME_JJA = month_to_season(zg_MME, "JJA")
x  = (TIME-1980)*1.
RC := regCoef_n(x,zg_MME_JJA,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_zg_MME_JJA = student_t(tval, df)
Trd_zg_MME_JJA = RC * 10.
copy_VarMeta(zg_MME_JJA(0,:,:,:), Trd_zg_MME_JJA)
copy_VarMeta(zg_MME_JJA(0,:,:,:), prob_zg_MME_JJA)


dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ts/ts_Amon_MPI-ESM_historical_rcp85_ensmean_185001-209912.nc"
f     = addfile(dir,"r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
ts_MME     = f->ts(iYYYY,{lat1:lat2},{lon1:lon2})
ts_MME_JJA     = month_to_season(ts_MME, "JJA")
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata,ts_MME_JJA&lat,ts_MME_JJA&lon)
ts_MME_JJA = mask(ts_MME_JJA,conform(ts_MME_JJA,lsm,(/1,2/)).ge.1,False)
meta := ts_MME_JJA(0,:,:)
RC := regCoef_n(x,ts_MME_JJA,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_ts_MME_JJA = student_t(tval, df)
Trd_ts_MME_JJA = RC * 10.
copy_VarMeta(meta, Trd_ts_MME_JJA)
copy_VarMeta(meta, prob_ts_MME_JJA)
print("MPI MME...")

;; for MPI #037--------------------------
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/zg/zg_Amon_MPI-ESM_historical_rcp85_r037i1850p3_185001-209912.nc"
f     = addfile(dir,"r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
zg037     := f->zg(iYYYY,{LEV},:,:)
zg037_JJA      = month_to_season(zg037, "JJA")
zg_diff = zg037_JJA - zg_MME_JJA
copy_VarMeta(zg037_JJA,zg_diff)
RC := regCoef_n(x,zg_diff,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_zg_diff := student_t(tval, df)
Trd_zg_diff = RC * 10.
copy_VarMeta(zg_diff(0,:,:,:), Trd_zg_diff)
copy_VarMeta(zg_diff(0,:,:,:), prob_zg_diff)


dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ts/ts_Amon_MPI-ESM_historical_rcp85_r037i1850p3_185001-209912.nc"
f    = addfile(dir,"r") 
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
ts037     = f->ts(iYYYY,{lat1:lat2},{lon1:lon2})
ts037_JJA = month_to_season(ts037, "JJA")

ts_diff     = ts037_JJA - ts_MME_JJA 
copy_VarMeta(ts037_JJA,ts_diff) 
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata,ts_diff&lat,ts_diff&lon)
ts_diff = mask(ts_diff,conform(ts_diff,lsm,(/1,2/)).ge.1,False)
meta = ts_diff(0,:,:)
RC := regCoef_n(x,ts_diff,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_ts_diff = student_t(tval, df)
Trd_ts_diff = RC * 10.
copy_VarMeta(meta, Trd_ts_diff)
copy_VarMeta(meta, prob_ts_diff)
delete([/zg_MME,ts_MME/])
print("MPI #037...")

;;for CESM MME--
dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/CESM2/post/Z3/b.e21.BHIST.SSP370.cmip6.f09_g17.LE2-ensmean.cam.h0.Z3.198001-202412.shift.plev.nc"
f     = addfile(dir,"r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
Z3_MME     = f->Z3(iYYYY,{llev},:,:)
Z3_MME_JJA = month_to_season(Z3_MME, "JJA")
x  = (TIME-1980)*1.
RC := regCoef_n(x,Z3_MME_JJA,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_Z3_MME_JJA = student_t(tval, df)
Trd_Z3_MME_JJA = RC * 10.
copy_VarMeta(Z3_MME_JJA(0,:,:,:), Trd_Z3_MME_JJA)
copy_VarMeta(Z3_MME_JJA(0,:,:,:), prob_Z3_MME_JJA)


dir   = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/CESM2/post/TS/b.e21.BHIST.SSP370.cmip6.f09_g17.LE2-ensmean.cam.h0.TS.198001-202412.shift.nc"
f     = addfile(dir,"r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
ts_cesmMME     = f->TS(iYYYY,{lat1:lat2},{lon1:lon2})
ts_cesmMME_JJA = month_to_season(ts_cesmMME, "JJA")

;;for mask --
; assume data is 3D (time,lat,lon)
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata := a->LSMASK
lsm := landsea_mask(lsdata,ts_cesmMME_JJA&lat,ts_cesmMME_JJA&lon)
ts_cesmMME_JJA = mask(ts_cesmMME_JJA,conform(ts_cesmMME_JJA,lsm,(/1,2/)).ge.1,False)
meta := ts_cesmMME_JJA(0,:,:)
RC := regCoef_n(x,ts_cesmMME_JJA,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_ts_cesmMME_JJA = student_t(tval, df)
Trd_ts_cesmMME_JJA = RC * 10.
copy_VarMeta(meta, Trd_ts_cesmMME_JJA)
copy_VarMeta(meta, prob_ts_cesmMME_JJA)

print("CESM MME...")

;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic10/Fig.E4.V2.obs.MPI.CESM2")
  res                               =True
  ;res@gsnPolar                      ="NH"
  font                              =0.025
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
  res@mpMinLatF                     =30
  res@mpMaxLatF                     =70.
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
  res@lbLabelBarOn                  = False
  res@pmLabelBarOrthogonalPosF      = 0.24
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
    opt@cnMonoFillScale        = False
    opt@cnFillScales           = (/0.7,0.75/)    ; change densities  
    opt@cnFillColors           = (/"white" ,"Transparent"/)
    opt@cnFillDotSizeF         =  0.001
    opt@cnInfoLabelOn          = False
;;============================================
    max_vaule = 20.
    interval  = 2.
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors = read_colormap_file(res@cnFillPalette)
    icol   = span_color_indexes(colors(2:,:), ncols) ;; 获取12个色号 
    icol = icol+2  
    icol(midcol1:midcol2) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors = icol
;;=============================================
  res@gsnRightString  = ""
  res@gsnLeftString   = ""
  plot_obs  = new(3, graphic)
  dplot_obs = new(3, graphic)

do i = 0, 1
    res@gsnLeftString = ""
    plot_obs(i)   = gsn_csm_contour_map(wks, Trd_Z_obs(i,:,:), res)  
    dplot_obs(i)  = gsn_csm_contour(wks, prob_Z_obs(i,:,:), opt)
    overlay(plot_obs(i), dplot_obs(i))
end do

    max_vaule = 0.4
    interval  = 0.05
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols) 
    icol = icol+2  
    icol(midcol1) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
    res@gsnLeftString = ""
    plot_obs(2)   = gsn_csm_contour_map(wks, SST_trd_COBE, res)  
    dplot_obs(2)  = gsn_csm_contour(wks, prob_SST_COBE, opt)
    overlay(plot_obs(2), dplot_obs(2))
;;for MPI===================================


    max_vaule = 20.
    interval  = 2.
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols) 
    icol := icol+2  
    icol(midcol1:midcol2) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
;;=============================================
  res@gsnRightString  = ""
  res@gsnLeftString   = ""
  plot_MPIMME  = new(3, graphic)
  dplot_MPIMME = new(3, graphic)

do i = 0, 1
    res@gsnLeftString = ""
    plot_MPIMME(i)   = gsn_csm_contour_map(wks, Trd_zg_MME_JJA(i,:,:), res)  
    dplot_MPIMME(i)  = gsn_csm_contour(wks, prob_zg_MME_JJA(i,:,:), opt)
    overlay(plot_MPIMME(i), dplot_MPIMME(i))
end do

    max_vaule = 0.4
    interval  = 0.05
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols) 
    icol = icol+2  
    icol(midcol1) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
    res@gsnLeftString = ""
    plot_MPIMME(2)   = gsn_csm_contour_map(wks, Trd_ts_MME_JJA, res)  
    dplot_MPIMME(2)  = gsn_csm_contour(wks, prob_ts_MME_JJA, opt)
    overlay(plot_MPIMME(2), dplot_MPIMME(2))

;;for CESM===================================
    max_vaule = 20.
    interval  = 2.
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols)  
    icol := icol+2  
    icol(midcol1:midcol2) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
;;=============================================
  res@gsnRightString  = ""
  res@gsnLeftString   = ""
  plot_CESMMME  = new(3, graphic)
  dplot_CESMMME = new(3, graphic)

do i = 0, 1
    res@gsnLeftString = ""
    plot_CESMMME(i)   = gsn_csm_contour_map(wks, Trd_Z3_MME_JJA(i,:,:), res)  
    dplot_CESMMME(i)  = gsn_csm_contour(wks, prob_Z3_MME_JJA(i,:,:), opt)
    overlay(plot_CESMMME(i), dplot_CESMMME(i))
end do

    max_vaule = 0.4
    interval  = 0.05
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols)  
    icol = icol+2  
    icol(midcol1) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
    res@gsnLeftString = ""
    plot_CESMMME(2)   = gsn_csm_contour_map(wks, Trd_ts_cesmMME_JJA, res)  
    dplot_CESMMME(2)  = gsn_csm_contour(wks, prob_ts_cesmMME_JJA, opt)
    overlay(plot_CESMMME(2), dplot_CESMMME(2))

;;for #037===================================
    max_vaule = 20.
    interval  = 2.
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols)  
    icol := icol+2  
    icol(midcol1:midcol2) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
;;=============================================
  res@gsnRightString  = ""
  res@gsnLeftString   = ""
  plot_diff  = new(3, graphic)
  dplot_diff = new(3, graphic)

do i = 0, 1
    res@gsnLeftString = ""
    plot_diff(i)   = gsn_csm_contour_map(wks,   Trd_zg_diff(i,:,:), res)  
    dplot_diff(i)  = gsn_csm_contour(wks, prob_zg_diff(i,:,:), opt)
    overlay(plot_diff(i), dplot_diff(i))
end do

    max_vaule = 0.4
    interval  = 0.05
    ncols = toint(2.*max_vaule/interval + 2.)
    midcol2 = toint(ncols/2.) 
    midcol1 = midcol2 - 1

    colors := read_colormap_file(res@cnFillPalette)
    icol   := span_color_indexes(colors(2:,:), ncols)  
    icol = icol+2  
    icol(midcol1) = -1
    res@cnLevelSelectionMode = "ManualLevels"   
    res@cnMinLevelValF       = -max_vaule       ;;850->10  ;500      
    res@cnMaxLevelValF       = max_vaule               
    res@cnLevelSpacingF      = interval         ;;850->2
    res@cnFillColors := icol
    res@gsnLeftString = ""
    plot_diff(2)   = gsn_csm_contour_map(wks, Trd_ts_diff, res)  
    dplot_diff(2)  = gsn_csm_contour(wks, prob_ts_diff, opt)
    overlay(plot_diff(2), dplot_diff(2))

gplt = new(12, graphic)

gplt((/0,4,8/)) = plot_obs
gplt((/1,5,9/)) = plot_MPIMME
gplt((/2,6,10/)) = plot_CESMMME
gplt((/3,7,11/)) = plot_diff
  pres                     =True
  pres@gsnPanelDebug      = False                   
  ; pres@gsnPanelXF         = (/ 0.0556078, 0.575608/)
  pres@txString            =""
  pres@gsnPanelXWhiteSpacePercent = 3.
  pres@gsnPanelYWhiteSpacePercent = 3.
  pres@gsnPanelLabelBar    = False
  pres@gsnPanelRowSpec     = False
  pres@pmLabelBarWidthF    = 0.6
  pres@pmLabelBarHeightF   = 0.07
  pres@lbLabelFontHeightF  = 0.013
  ; ; pres@lbTitleString      ="units:  m^2*PVU/s"
  gsn_panel(wks,gplt,(/3,4/),pres)
  ; gsn_panel(wks,plot_obs,(/3,1/),pres)
  wallClockElapseTime(TotalTime, "totaltime", 0)  
end