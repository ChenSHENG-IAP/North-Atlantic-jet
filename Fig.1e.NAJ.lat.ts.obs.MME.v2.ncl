load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;;==========================
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
LEV   = 25000
lat1 = 30.
lat2 = 70.
lon1 = 0.
lon2 = 360.
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)


filout = "NAJ.lat.JJA.models.1-6.historical.45W-0.v2.nc"
dirout = "/home/users/shengc01/work11/scripts1/data/"
fout = addfile(dirout + filout, "r")
NAJ_lat = fout->NAJ_lat 

dirout = "/home/users/shengc01/work11/scripts1/data/"
filout = "NAJ.lat.JJA.v2.nc"
fout = addfile(dirout + filout, "r")

MERRA2_jet_lat1 = fout->MERRA2_jet_lat1(0:nt-1)
ERA5_jet_lat1   = fout->ERA5_jet_lat1(0:nt-1)
JRA55_jet_lat1  = fout->JRA55_jet_lat1(0:nt-1)
missing = -9999.
Y_jet = new((/4,nt/), float,missing)
Y_jet(0,:) = tofloat(MERRA2_jet_lat1)
Y_jet(1,:) = tofloat(JRA55_jet_lat1)
Y_jet(2,:) = tofloat(ERA5_jet_lat1)
Y_jet(3,0:34) = tofloat(dim_avg_n(NAJ_lat, 1))
Y_jet(3,35:)  = missing 
; Y_jet = dim_rmvmean_n(Y_jet, 1)

; ;**********regression*************************
x  = (TIME-1980)*1.
RC = regCoef_n(x,Y_jet,0,1)
df = onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval = onedtond(RC@tval,dimsizes(RC))
prob = student_t(tval, df)
a  = RC
b  = RC@yintercept

ytrd = Y_jet

printVarSummary(x)
printVarSummary(a)

ytrd(0,:) = a(0) * x + b(0)
ytrd(1,:) = a(1) * x + b(1)
ytrd(2,:) = a(2) * x + b(2)
ytrd(3,:) = a(3) * x + b(3)
ytrd(3,35:) = missing
print(a * 10.)
ymean = dim_avg_n_Wrap(dim_avg_n_Wrap(Y_jet(0:2,:), 0),0)
print(prob)

;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic1/Fig1.NAJ.lat.ts.obs.MME.v2")
;; for line=========
xais = ispan(1, nt, 1)
yr   = ispan(YrStart, YrEnd, 1)
font                              =0.015
;; for first plot==========
resL = True 
resL@gsnMaximize                   =False
resL@gsnDraw                       =False
resL@gsnFrame                      =False
resL@vpHeightF= 0.4                    ; change aspect ratio of plot
resL@vpWidthF = 0.8  
resL@gsnLeftString = ""
resL@vpXF     = 0.1
; resL@gsnRightString =  "trend(P="+tostring_with_format(prob, "%3.2f")+")";;sprintf("%3.2f",escorc(N3,SI))+"***"
resL@gsnStringFontHeightF = font
resL@tmYLLabelFontHeightF          = font
resL@tiYAxisOn   = True
resL@tiYAxisString = ""
resL@tmXTOn     = False  
resL@tmYROn     = False                                  
resL@trXMinF           = min(xais)             
resL@trXMaxF           = max(xais) 
resL@trYMinF  =  - 7. + ymean                  
resL@trYMaxF  =   12.  + ymean                 
;----------------------------------
resL@xyLineThicknesses  = (/4.,4.,4.,4./)
; resL@xyLineOpacities = (/0.7,0.7,0.7,0.7/)
; resL@xyLineColors = (/"black","pink1","cornflowerblue","red"/)
; resL@xyLineColors = (/"cornflowerblue","black","red","orange"/)
resL@xyLineColors = (/"cornflowerblue","orange","firebrick1","black"/)
resL@xyDashPattern=0
resL@gsnYRefLine        = ymean  
resL@gsnYRefLineThicknessF = 2.
resL@gsnYRefLineColor  = "black"
resL@gsnYRefLineDashPattern = 0
;;===小刻度=========================
resL@tmXBMinorOn     = True
resL@tmXBMinorValues  = xais
;;for XB label=======================
resL@tmXBMode         = "Explicit"
resL@tmXBValues        = xais(0::5)
resL@tmXBLabels        = yr(0::5)
resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
; ;;for YL label=======================
; resL@tmXBMode         = "Explicit"
; resL@tmXBValues        = ispan(29, 35, 1)
; resL@tmXBLabels        = ""+resL@tmXBValues+"N"
; resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
resL@xyCurveDrawOrder = "PostDraw"
resl = resL
; ;;for bar =====
;  resL@gsnXYBarChart         = True            ; create bar chart 
;  resL@gsnAboveYRefLineColor = "red"           ; above ref line fill red
;  resL@gsnBelowYRefLineColor = "Lightblue"          ; below ref line fill blue


resL@pmLegendDisplayMode="Always"
resL@lgPerimOn          = False
resL@xyExplicitLegendLabels= (/" MERRA2    "," JRA55        "," ERA5         "," ALL(MME) "/) + "("+sprintf("%3.2f", a*10. )+"~S~o~N~ 10a~S~-1~N~)" 
resL@pmLegendWidthF=0.1
resL@pmLegendHeightF=0.12
resL@pmLegendOrthogonalPosF=-1.2 
resL@pmLegendParallelPosF=0.65
resL@lgLabelFontHeightF = 0.015  
resL@lgLabelFontColor   = "black"
resL@gsnLeftString = "";;"NAJ lat"
resL@gsnRightString = "";;"P<0.05"
plot0 = gsn_csm_xy (wks,xais,Y_jet,resL)
 resl@gsnLeftString = ""
 resl@gsnRightString = ""
 resl@xyLineThicknesses = (/2.,2.,2.,2./)
 resl@xyDashPattern = 1

;--------legend
 plot2_1 = gsn_csm_xy (wks,xais,ytrd,resl)
overlay(plot0, plot2_1)
draw(plot0)
frame(wks)

end