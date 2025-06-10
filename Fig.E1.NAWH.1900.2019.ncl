load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
  ihp    = 0     ;A scalar indicating the low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
  sigma  = 1.0     ;A scalar indicating the power of the sigma factor (nsigma >= 0). Note: nsigma=1. is common.
  nWgt   = 11     
  fca    = 1./9.  
  fcb    =-999.     ;A scalar used only when a band-pass filter is desired. It is the second cut-off frequency (fca < fcb < 0.5).
  wgt    = filwgts_lanczos(nWgt,ihp,fca,fcb,sigma)
;;==========================

;;==========================
TotalTime = systemfunc("date")
YrStart = 1900
YrEnd   = 2019
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)





fil = "NAWH.4index.JJA.COBE.1891-2022.nc"
dir = "/home/users/shengc01/work11/scripts10/data/"
fout = addfile(dir + fil, "r")
ist = YrStart - 1891
ied = ist + nt -1
NAWH_Box = fout->NAWH_Box(ist:ied)
NAWH_BoxMinNH = fout->NAWH_BoxMinNH(ist:ied)
NAWH_BoxMinNA = fout->NAWH_BoxMinNA(ist:ied)
NAWH_SSTGbased = fout->NAWH_SSTGbased(ist:ied) 



NAWHI = new((/4,nt/), float)
NAWHI(0,:) = tofloat(NAWH_BoxMinNA) * -1.
NAWHI(1,:) = tofloat(NAWH_BoxMinNH) * -1.
NAWHI(2,:) = tofloat(NAWH_SSTGbased) * -1.
NAWHI(0:2,:)  = dim_standardize(NAWHI(0:2,:), 1)

mean = dim_avg_n_Wrap(NAWHI(0:2,:), 0)
NAWHI(3,:) = tofloat(mean) 

NAWHI     = runave_n_Wrap(NAWHI, 10, 1, 1)
NAWHI_dtrd     = dtrend_msg_n(TIME, NAWHI, True, False, 1)

x_recent = ispan(1980, 2019, 1)*1.
RC_raw = regCoef_n(x_recent,NAWHI(:,80::),0,1)
RC_dtr = regCoef_n(x_recent,NAWHI_dtrd(:,80::),0,1)

print(RC_raw)
print(RC_dtr)
print((RC_dtr/RC_raw)*100+"%")

NAWHI     = dim_standardize(NAWHI, 1)
NAWHI_dtrd     = dim_standardize(NAWHI_dtrd, 1)


; ;**********regression*************************
x  = (TIME-YrStart)*1.
RC := regCoef_n(x,NAWHI,0,1)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob := student_t(tval, df)
a  = RC
b  = RC@yintercept
ytrd = NAWHI
ytrd(0,:) = a(0) * x + b(0)
ytrd(1,:) = a(1) * x + b(1)
ytrd(2,:) = a(2) * x + b(2)
ytrd(3,:) = a(3) * x + b(3)


;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic10/Fig.E1.NAWH."+YrStart+"."+YrEnd)

;; for line=========
xais = ispan(1, nt, 1)
yr   = ispan(YrStart, YrEnd, 1)
font = 0.02
;; for first plot==========
resL = True 
resL@gsnMaximize                   =False
resL@gsnDraw                       =False
resL@gsnFrame                      =False
resL@vpHeightF= 0.4                    ; change aspect ratio of plot
resL@vpWidthF = 0.7  
resL@gsnLeftString = ""
; resL@vpXF     = 0.1
; resL@gsnRightString =  "trend(P="+tostring_with_format(prob, "%3.2f")+")";;sprintf("%3.2f",escorc(N3,SI))+"***"
resL@gsnStringFontHeightF = font
resL@tmYLLabelFontHeightF          = font
resL@tiYAxisOn   = True
resL@tiYAxisString = ""
resL@tmXTOn     = False  
resL@tmYROn     = False                                  
resL@trXMinF           = min(xais)             
resL@trXMaxF           = max(xais) 
resL@trYMinF  =  - 4.                  
resL@trYMaxF  =   4.                  
;----------------------------------
resL@xyLineThicknesses  = (/4.,4.,4.,4./)
; resL@xyLineOpacities = (/0.7,0.7,0.7,1./)
; resL@xyLineColors = (/"black","pink1","cornflowerblue","red"/)
resL@xyLineColors = (/"gray","gray80","gray60","red"/)
;resL@xyLineColors = (/"Transparent","Transparent","red","black"/)
; resL@gsnXYFillOpacities = 0.1
resL@xyDashPattern=0
resL@gsnYRefLine        = 0.  
resL@gsnYRefLineThicknessF = 1.
resL@gsnYRefLineColor  = "black"
resL@gsnYRefLineDashPattern = 0
;;===小刻度=========================
resL@tmXBMinorOn     = True
resL@tmXBMinorValues  = xais(0::10)
;;for XB label=======================
resL@tmXBMode         = "Explicit"
resL@tmXBValues        =  (/1900,1920,1940,1960,1980,2000,2019/)-1900 + 1  ;;xais(0::20)
resL@tmXBLabels        = (/1900,1920,1940,1960,1980,2000,2019/) 
resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
; ;;for YL label=======================
; resL@tmXBMode         = "Explicit"
; resL@tmXBValues        = ispan(29, 35, 1)
; resL@tmXBLabels        = ""+resL@tmXBValues+"N"
; resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
resl = resL

resL@pmLegendDisplayMode="Always"
resL@lgPerimOn          = False
resL@xyExplicitLegendLabels= (/" NAWH_NA"," NAWH_NH"," NAWH"," MEAN"/) ;;+ "("+sprintf("%3.2f", a*10. )+"~S~o~N~ 10a~S~-1~N~)" 
; resL@xyExplicitLegendLabels= (/" "," ","   NAWH","  NAJ lat"/) ;;+ "("+sprintf("%3.2f", a*10. )+"~S~o~N~ 10a~S~-1~N~)" 
resL@pmLegendWidthF=0.1
resL@pmLegendHeightF=0.12
resL@pmLegendOrthogonalPosF=-1.2 
resL@pmLegendParallelPosF=0.2
resL@lgLabelFontHeightF = 0.016  
resL@lgLabelFontColor   = "black"
resL@gsnLeftString = "(a)summer NAWH";;"NAWH intensity"
resL@gsnRightString = "";;"JJA"
plot0 = gsn_csm_xy (wks,xais,NAWHI,resL)

resl@gsnLeftString = ""
resl@gsnRightString = ""
resl@xyLineThicknesses = (/2.,2.,2.,2./)
resl@xyDashPattern = 1

;--------legend  
delete(resl@xyLineColors)
resl@xyLineColor = "red"
plot0_1 = gsn_csm_xy (wks,xais,ytrd(3,:),resl)
overlay(plot0, plot0_1)

resL@gsnLeftString = "(b)summer NAWH detrend";;"NAWH intensity"
plot1  = gsn_csm_xy (wks,xais,NAWHI_dtrd,resL)
resl@xyLineColor = "Transparent"
resl@gsnBelowYRefLineColor = "LightBlue"
resl@gsnAboveYRefLineColor = "LightPink"
plot1_1  = gsn_csm_xy (wks,xais,NAWHI_dtrd(3,:),resl)
overlay(plot1,plot1_1)
gplt   = new(2, graphic)
gplt(0) = plot0 
gplt(1) = plot1

  pres                     =True
  gsn_panel(wks,gplt,(/1,2/),pres)
draw(plot0)
frame(wks)
end
