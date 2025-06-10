load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;;==========================
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019

TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)

model = "MPI"   ;;CanESM2  MK36  "CM3"  MPI   CESM1
dir = "/home/users/shengc01/work11/scripts2/data/"
fil = "NAJ.lat.JJA."+model+".LE.45W-0.v2.nc"
f = addfile(dir + fil, "r")

time   = f->time
YYYY   = cd_calendar(time,-1)/100    ; entire file
iYYYY  = ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
NAJ_lat  = f->NAJ_lat(iYYYY,:)
nmem  = dimsizes(NAJ_lat(0,:))
nt  = dimsizes(NAJ_lat(:,0))
member = ispan(0, nmem-1, 1)
print(iYYYY)


NAJ_internal = dim_rmvmean_n_Wrap(NAJ_lat, 1)
NAJ_ensmean  = dim_avg_n_Wrap(NAJ_lat,1)
; ;**********sort*************************
;x  = (TIME-1980)*1.
x = TIME(1980- YrStart::)
RC = regCoef_n(x,NAJ_internal(1980- YrStart::,: ),0,0)
RC!0 = "member"
RC&member = member
; print("Raw:"+RC)

qsort(RC)
member_sort = RC&member
print("Sort:"+RC+".."+member_sort)

Nselet = 15
E15 = member_sort(0:Nselet-1)
P15 = member_sort(nmem-Nselet:nmem-1)

print("E15:" + sprinti("%03d", (E15+1)))
print("P15:" + sprinti("%03d", (P15+1)))
maxRC = max(RC) * 10.
minRC = min(RC) * 10.
NAJ_MME = new((/3,nt/), float)
NAJ_MME(0,:) = NAJ_ensmean
NAJ_MME(1,:)  = dim_avg_n_Wrap(NAJ_internal(:,(/P15/)), 1)
NAJ_MME(2,:)  = dim_avg_n_Wrap(NAJ_internal(:,(/E15/)), 1)

opt   = 1

NAJ_internal = runave_n(NAJ_internal, 9, opt, 0) ;;NAJinternal(time,member)
NAJ_MME = runave_n(NAJ_MME, 9, opt, 1) ;;NAJMME(3,time)

NAJ_MME = dim_rmvmean_n(NAJ_MME, 1)  ;; NAJ
NAJ_internal = dim_rmvmean_n_Wrap(NAJ_internal,0)

NAJ_MME!0="group"
NAJ_MME&group = ispan(0, 2, 1)


;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic2/Fig.2.NAJlat.ts.LE.v2."+model)
;; for line=========
xais = ispan(1, nt, 1)
font                              =0.028
;; for first plot==========
resL = True 
resL@gsnMaximize                   =False
resL@gsnDraw                       =False
resL@gsnFrame                      =False
resL@vpHeightF= 0.4                    ; change aspect ratio of plot
resL@vpWidthF = 0.8  
resL@gsnLeftString = ""
resL@vpXF = 0.15
; resL@gsnRightString =  "trend(P="+tostring_with_format(prob, "%3.2f")+")";;sprintf("%3.2f",escorc(N3,SI))+"***"
resL@gsnStringFontHeightF = font
resL@tmYLLabelFontHeightF          = font
resL@tiYAxisOn   = True
resL@tiYAxisString = ""
resL@tmXTOn     = False  
resL@tmYROn     = False                                  
resL@trXMinF           = min(xais)             
resL@trXMaxF           = max(xais) 
; resL@trYMinF  =  - 5. + ymean                  
; resL@trYMaxF  =   7.  + ymean                 
;----------------------------------
; resL@xyLineThicknesses  = (/3.,3.,3./)

opa = new(nmem, float)
color = new(nmem,string)
opa = 0.2
color="grey70"
; color(P15) = "red"
; color(E15) = "blue"
resL@xyLineColors = color
; resL@xyLineOpacities = opa

resL@xyDashPattern=0
resL@gsnYRefLine        = 0.  
resL@gsnYRefLineThicknessF = 3.
resL@gsnYRefLineColor  = "black"
resL@gsnYRefLineDashPattern = 0
;;===小刻度=========================
resL@tmXBMinorOn     = True
resL@tmXBMinorValues  = xais
;;for XB label=======================
resL@tmXBMode         = "Explicit"
resL@tmXBValues        = xais(0::5)
resL@tmXBLabels        = TIME(0::5)
resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
; ;;for YL label=======================
; resL@tmXBMode         = "Explicit"
; resL@tmXBValues        = ispan(29, 35, 1)
; resL@tmXBLabels        = ""+resL@tmXBValues+"N"
; resL@tmXBLabelFontHeightF = font
; resL@tmXBLabelAngleF      = 45
resl = resL
; ;;for bar =====
;  resL@gsnXYBarChart         = True            ; create bar chart 
;  resL@gsnAboveYRefLineColor = "red"           ; above ref line fill red
;  resL@gsnBelowYRefLineColor = "Lightblue"          ; below ref line fill blue


plot0 = gsn_csm_xy (wks,xais,NAJ_internal({member|:},{time|:}),resl)
 resl@gsnLeftString = ""
 resl@gsnRightString = ""
 ; resl@xyLineThicknesses = (/3.,3.,3./)
 resl@xyDashPattern = 1

;--------legend
resl@pmLegendDisplayMode="Always"
resl@lgPerimOn          = False 
resl@xyExplicitLegendLabels= (/" mean"," P"+Nselet+"  ("+sprintf("%2.1f",maxRC)+"~S~o~N~ 10a~S~-1~N~)"," E"+Nselet+"  ("+sprintf("%2.1f",minRC)+"~S~o~N~ 10a~S~-1~N~)"/) 
resl@xyLineThicknesses  = (/6.,6.,6./)
resl@xyDashPatterns=(/0.,0.,0./)
resl@xyLineColors := (/"black","red","cornflowerblue"/)
resl@pmLegendWidthF=0.1
resl@pmLegendHeightF=0.1
resl@pmLegendOrthogonalPosF=-1.28 
resl@pmLegendParallelPosF=0.7
resl@lgLabelFontHeightF = 0.02  
resl@lgLabelFontColor   = "black"
resl@xyLineColor=  "red"
resl@gsnLeftString = "";;"NAJ lat "+ model +"("+nmem+")"
resl@gsnRightString = ""
resl@xyCurveDrawOrder = "PostDraw"
;;resl@xyLineOpacities :=(/1.,0.7,0.7/)
 plot2_1 = gsn_csm_xy (wks,xais,NAJ_MME({group|:},{time|:}),resl)
overlay(plot0, plot2_1)

draw(plot0)
frame(wks)




end