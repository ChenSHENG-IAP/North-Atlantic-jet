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
LEV  = 1000
lev1 = 1000
lev2 = 100  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)


dir   = "/data6/shengc2017/data/ERA5/PV_flux/season/JJA/"
fils  = systemfunc("ls "+dir + "*{"+YrStart+".."+YrEnd+"}.nc")
fs    = addfiles(fils, "r")
U     = fs[:]->U(:,{lev1:lev2},{lat1:lat2},{lon1:lon2})
U1000_NAmean = dim_avg_n_Wrap(U(:,{LEV},:,:),2)

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
xais  = U1000_NAmean&latitude
;;======================
wks   = gsn_open_wks ("eps", "/data6/shengc2017/work9_sub_Jet/NCC/pic1/Fig.E6.obs.trd.U1000line.v3")
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
resL@trYMinF  =  -0.5
resL@trYMaxF  =   0.2
;resL@tmYLFormat        ="0@;*.1f"       ; 调整数字1  有效数字


resL@tmYLLabelsOn = False 
resL@tmYLOn     = False
resL@tmYRLabelFontHeightF          = font
resL@tmYRLabelsOn = True
resL@tmYUseLeft    = False
resL@tmYRAutoPrecision    = False          ; no auto precision
resL@tmYRPrecision        = 1              ; set the precision
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
resL@tmXBMinorOn     = True
resL@tmXBMinorValues  = xais(::2)
;;for XB label=======================
resL@tmXBLabelFontHeightF = font
resL@tmXBTickSpacingF    = 10
 resl = resL

;;for bar =====
 resL@gsnXYBarChart         = False            ; create bar chart
 ; resL@gsnAboveYRefLineColor = "white"           ; above ref line fill red
 ; resL@gsnBelowYRefLineColor = "LightBlue"          ; below ref line fill blue

 resL@gsnAboveYRefLineColor = "LightPink"           ; above ref line fill red
 resL@gsnBelowYRefLineColor = "LightBlue"           ; above ref line fill red
 resL@gsnXYBarChartOutlineThicknessF = 0.


 ; resL@gsnLeftString = "(e)Us(shading)&eddy conv vint(line)"
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