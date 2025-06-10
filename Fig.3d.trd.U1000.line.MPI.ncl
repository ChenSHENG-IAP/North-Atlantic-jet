load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 30.
lat2 = 70.
lon1 = 315.
lon2 = 360.
lev1  = 1000. * 100
lev2  = 100. * 100
LEV   = 1000. * 100 
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

;;for U


    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/"
    f   = addfile(dir+"ua_Amon_MPI-ESM_historical_rcp85_E15.ensmean_185001-209912.v2.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    UEO     = f->ua(iYYYY,{lev1:lev2},{lat1:lat2},{lon1:lon2})  
    UE     = month_to_season(UEO,"JJA")

    dir = "/gws/nopw/j04/ncas_climate_vol2/users/shengc02/LE-post/MPI/ua/"
    f   = addfile(dir+"ua_Amon_MPI-ESM_historical_rcp85_P15.ensmean_185001-209912.v2.nc", "r")
    time   = f->time
    YYYY   = cd_calendar(time,-1)/100    ; entire file
    iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
    UPO     = f->ua(iYYYY,{lev1:lev2},{lat1:lat2},{lon1:lon2})  
    UP    = month_to_season(UPO,"JJA")
    Udiff     = UE - UP 
    copy_VarMeta(UE,Udiff)
    U1000_NAmean = dim_avg_n_Wrap(Udiff(:,{LEV},:,:),2)

print("U..")
; ;**********regression*************************
x  = (TIME-1980) * 1.
RC := regCoef_n(x,U1000_NAmean,0,0)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob_U1000_NAmean = student_t(tval, df)
U1000_NAmean_trd = RC * 10. 
copy_VarMeta(U1000_NAmean(0,:), U1000_NAmean_trd)
copy_VarMeta(U1000_NAmean(0,:), prob_U1000_NAmean)
U1000_NAmean_trd = where(ismissing(U1000_NAmean(0,:)), U1000_NAmean@_FillValue, U1000_NAmean_trd)
prob_U1000_NAmean = where(ismissing(U1000_NAmean(0,:)), U1000_NAmean@_FillValue, prob_U1000_NAmean)
xais  = U1000_NAmean&lat
;;======================
wks   = gsn_open_wks ( "eps", "/home/users/shengc01/work11/pic7/Fig.3.trd.U1000.line.MPI")
font                              =0.015
;; for first plot==========
resL = True
resL@gsnMaximize                   =False
resL@gsnDraw                       =False
resL@gsnFrame                      =False
resL@vpHeightF= 0.5                    ; change aspect ratio of plot
resL@vpWidthF = 0.6
resL@gsnLeftString = ""
resL@gsnRightString =  "";;sprintf("%3.2f",escorc(N3,SI))+"***"
resL@gsnStringFontHeightF = font

resL@tiYAxisOn   = False
resL@tiYAxisString = ""
resL@tmXTOn     = False
resL@tmYROn     = True
resL@trXMinF           = min(xais)
resL@trXMaxF           = max(xais)
resL@trYMinF  =  -0.4
resL@trYMaxF  =   0.4
resL@tmYRFormat        ="0@;*.1f"       


resL@tmYLLabelsOn = False 
resL@tmYLOn     = False
resL@tmYRLabelFontHeightF          = font
resL@tmYRLabelsOn = True
resL@tmYUseLeft    = False
resL@tmYRAutoPrecision    = False          ; no auto precision
; resL@tmYRPrecision        = 1              ; set the precision
;----------------------------------
resL@xyLineThicknessF  = 0.001
resL@xyLineColor = "black"
resL@xyDashPattern=0
resL@gsnYRefLine        = 0.
resL@gsnYRefLineThicknessF = 2.
resL@gsnYRefLineColor  = "black"
resL@gsnYRefLineDashPattern = 0

resL@gsnXRefLine        = 50
resL@gsnXRefLineThicknessF = 1.5
resL@gsnXRefLineColor  = "black"
resL@gsnXRefLineDashPattern = 14
;;===小刻度=========================
resL@tmXBMinorOn     = False
resL@tmXBMinorValues  = xais(::2)
;;for XB label=======================
resL@tmXBLabelFontHeightF = font
resL@tmXBTickSpacingF    = 10
; ;;for YL label=======================
 resl = resL

;;for bar =====
 resL@gsnXYBarChart         = False            ; create bar chart

 resL@gsnAboveYRefLineColor = "LightPink"           ; above ref line fill red
 resL@gsnBelowYRefLineColor = "LightBlue"           ; above ref line fill red
 resL@gsnXYBarChartOutlineThicknessF = 0.


 resL@gsnLeftString = ""
 resL@gsnRightString =  ""
 plot = gsn_csm_xy(wks,xais,U1000_NAmean_trd,resL)

 resl@gsnLeftString = ""
 resl@gsnRightString = ""
 resl@xyLineThicknesses = 5.
 resl@xyLineColor = "gray"
 resl@xyDashPattern = 0

 plot_1 = gsn_csm_xy (wks,xais,U1000_NAmean_trd,resl)
 ; overlay(plot,plot_1)
 


 draw(plot)
 frame(wks)

end