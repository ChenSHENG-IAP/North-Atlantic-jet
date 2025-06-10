load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;;==========================
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2019
lat1 = 0.
lat2 = 90.
lon1 = 0.
lon2 = 360.
LEV  = 25000
lev1 = 100000
lev2 = 7000  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)
;; select member
model = "MPI"   ;;CanESM2  MK36  "CM3"  MPI  CESM1
dir = "/home/users/shengc01/work11/scripts2/data/"
fil = "NAJ.lat.JJA."+model+".LE.45W-0.v2.nc"
f = addfile(dir + fil, "r")
time   := f->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
NAJ_lat  = f->NAJ_lat(iYYYY,:)  ;; year, member
NAJ_lat = dim_rmvmean_n(NAJ_lat, 1)
member = ispan(0,dimsizes(NAJ_lat(0,:))-1,1)
nmem   = dimsizes(member)
; ;**********sort*************************
x  = (TIME-1980)*1.
RC = regCoef_n(x,NAJ_lat,0,0)
RCO = RC
RC!0 = "member"
RC&member = member
qsort(RC)
member_sort = RC&member
Nselet = 15
E15 = member_sort(0:Nselet-1)
P15 = member_sort(nmem-Nselet:nmem-1)


dirout = "/home/users/shengc01/work11/scripts4/data/"
filout = "SSTG.index.JJA.MPI.LE.1860-2019.v3.nc"
fout = addfile(dirout + filout, "r")
time   := fout->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
SSTG = fout->SSTG2(iYYYY,:)

RC_SSTG = regCoef_n(x,SSTG,0,0)
RC_SSTG!0 = "member"
RC_SSTG&member = member

dirout = "/home/users/shengc01/work11/scripts4/data/"
filout = "IPO.index.JJA.MPI.LE.1860-2019.v3.nc"
fout = addfile(dirout + filout, "r")
time   := fout->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
IPOI = fout->IPO2(iYYYY,:)
IPOI = dim_rmvmean_n(IPOI,1)
RC_IPOI = regCoef_n(x,IPOI,0,0)
RC_IPOI!0 = "member"
RC_IPOI&member = member


NAJ = RCO * 10.
SSTG_trd   = RC_SSTG * 100.



NAJ = dim_standardize(NAJ, 1)
SSTG_trd = dim_standardize(SSTG_trd, 1)
yy = new((/2,nmem/), float)
yy(0,:) = (/NAJ/)
yy(1,:) = tofloat(SSTG_trd)
corr   = escorc_n(SSTG_trd,NAJ,0,0)
pvalue_cor = rtest(corr, dimsizes(NAJ), 0)
print("CC between SSTG&NAJ="+corr)
print("P-value of fig(a)="+pvalue_cor)
x1 = fspan(-2., 2, 40)
y1 = corr*x1 



T2M_num_p1 = num(SSTG_trd .ge. 0 .and. NAJ .ge. 0)
T2M_num_p2 = num(SSTG_trd .ge. 0 .and. NAJ .le. 0)
T2M_num_p3 = num(SSTG_trd .le. 0 .and. NAJ .le. 0)
T2M_num_p4 = num(SSTG_trd .le. 0 .and. NAJ .ge. 0)

pct_T2M    = toint(100. * (T2M_num_p2 + T2M_num_p4)/nmem )



;;=============================================

wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic8/Fig.E7.scatter.plt.NAJ.NAWH.Line.NAJ.NAWH.piControl" )      
res                    = True


  res@gsnDraw            = False       ; Don't draw plot
  res@gsnFrame           = False       ; Don't advance frame
  font                              = 0.015
  res@gsnLeftStringFontHeightF      = font
  res@gsnRightStringFontHeightF     = font
  res@tmYLLabelFontHeightF          = font
  res@tmXBLabelFontHeightF          = font
;
; Don't use gsnMaximize. Instead, we will control 
; exactly where plot goes.
;
  res@vpXF              = 0.057
  res@vpYF              = 0.85
  res@vpWidthF          = 0.32
  res@vpHeightF         = 0.2

;---Resources for X/Y axes.

  res@trYMinF            =  -3.               ; min value on y-axis
  res@trYMaxF            =  3.               ; max value on y-axis
  res@trXMinF            =  -3.                  ; min value on x-axis
  res@trXMaxF            =  3.                  ; max value on x-axis

res@tmYLFormat        ="0@;*.1f"       
res@tmXBFormat        ="0@;*.1f"       
res@tmYLLabelFontHeightF = font
res@tmXBLabelFontHeightF = font
ndot  = nmem
color = new(ndot, string)
marker = new(ndot, integer)


irec_follow = ind(NAJ*SSTG_trd .lt. 0 )
irec_unfollow = ind(NAJ*SSTG_trd .ge. 0)

color(irec_follow) =  "red" 
color(irec_unfollow) = "black"

marker = 16

;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.005
  res@xyMarkerThicknessF = 0.5
  res@xyMonoMarkerColor  = False          ; Allow different colors for markers
  ; res@xyMarkerColors    = color           ; It's okay to list more than
  res@xyMarkerColors    = "black"           ; It's okay to list more than
                                          ; you need here
  res@xyMarkers          = marker

;---Tickmark resources
  res@tmXTOn           = False
  res@tmYROn           = False
  res@tmXBMode         = "Manual"
  res@tmYLMode         = "Manual"
  res@tmXBTickSpacingF = 1.
  res@tmYLTickSpacingF = 1.

;---Title resources
  res@tiXAxisString      = "";;"NAWH trend"
  res@tiYAxisString      = "";;"NAJ trend"
  res@tiMainFontHeightF  = font
  res@gsnXRefLine        = 0.0
  res@gsnYRefLine        = 0.0
  res@gsnXRefLineThicknessF = 1.5
  res@gsnYRefLineThicknessF = 1.5
;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------
  res@gsnLeftString = "";;"(a)trend of NAWH&NAJ"
  res@gsnRightString = "";;"corr="+tostring_with_format(corr(0), "%3.2f")+"***"
  plot1 = gsn_csm_xy (wks,transpose((/SSTG_trd,SSTG_trd/)),transpose((/NAJ,NAJ/)),res)

; for y=kx
delete([/res@xyMarkLineMode,res@gsnLeftString,res@gsnRightString/])
res@xyLineThicknessF  = 5.
res@xyLineColor = "gray"
res@xyDashPattern=0
res@xyCurveDrawOrder = "PreDraw"
plot1_1 = gsn_csm_xy (wks,x1,y1,res)
overlay(plot1, plot1_1)
draw(plot1)

;; for pic 2 ---
;;==========================
TotalTime = systemfunc("date")
YrStart = 0
YrEnd   = 1000
TIME  := ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME) 

dirout = "/home/users/shengc01/work11/scripts3/data/"
filout = "NAJ.lat.JJA.MPI.piControl.nc"
fout = addfile(dirout + filout, "r")
NAJ  := fout->NAJ_lat_MPI(YrStart:YrEnd)
nt   = dimsizes(NAJ)

dirout = "/home/users/shengc01/work11/scripts6/data/"
filout = "SSTG.index.JJA.MPI.piControl.nc"
fout = addfile(dirout + filout, "r")
NAWHI  = fout->SSTG(YrStart:YrEnd)
nt   = dimsizes(NAJ)

missing = -9999.
yyy = new((/2,nt/), float,missing)
yyy(0,:) = NAJ
yyy(1,:) = NAWHI

yyy = dim_rmvmean_n_Wrap(yyy, 1)
yyy = runave_n_Wrap(yyy, 40, 1, 1)
; yyy    = wgt_runave_n_Wrap(yyy,wgt,1,1)
yyy = dim_standardize_n_Wrap(yyy, 1, 1)
R_NAWHI_NAJ = escorc(yyy(1,:), yyy(0,:))
p_R_NAWHI_NAJ = rtest(R_NAWHI_NAJ, nt, 0)

print("--->CC bwt NAWHI and NAJ="+R_NAWHI_NAJ)
print("--->P-value="+p_R_NAWHI_NAJ)
;; for line=========
xais = ispan(YrStart, YrEnd, 1)
yr   = xais

;; for fourth plot===========================
    resL = True  
    resL@gsnDraw            = False       ; Don't draw plot
    resL@gsnFrame           = False       ; Don't advance frame
  resL@vpXF              = res@vpXF + res@vpWidthF + 0.08
  resL@vpYF              = res@vpYF
  resL@vpWidthF         = res@vpWidthF + 0.2
  resL@vpHeightF        = res@vpHeightF   
;;===小刻度=========================
    resL@tmXBMinorOn     = False
    resL@tmXBMinorValues  = xais
;;for XB label=======================
    resL@tmXBMode         = "Explicit"
    resL@tmXBValues        = xais(0::200)
    resL@tmXBLabels        = xais(0::200)
    resL@tmXBLabelAngleF   = 45
;;==============================
    resL@gsnStringFontHeightF = font 
    resL@gsnDraw     = False
    resL@gsnFrame    = False
    resL@tiYAxisOn   = False
    resL@tmXTOn      = False
    resL@tmYROn      = False
    resL@tmYLLabelFontHeightF = font 
    resL@tmXBLabelFontHeightF = font
    ; resL@gsnYRefLine = 0.
  
    resL@gsnLeftString = "";;"(b)corr="+sprintf("%3.2f", R_NAWHI_NAJ)+"***"
    resL@gsnRightString = "";;;;"piControl" 
                
    resL@trXMinF           = min(xais)             
    resL@trXMaxF           = max(xais) 
    resL@trYMinF  = -3.5                   ; min value on y-axis
    resL@trYMaxF  =  3.5                   ; max value on y-axis
    ; resL@tmYLTickSpacingF = 1. 
;;----------------------------------
    resL@xyLineThicknesses=(/2,2/)
    resL@xyLineColors=(/"firebrick1","grey30"/)
    ; resL@xyLineColors=(/"firebrick1","cornflowerblue"/)
    resL@xyDashPatterns=(/0,0/)
    ;resL@xyLineOpacities = (/1.,1./)
    resL@gsnYRefLine        = (/0./)  
    resL@gsnYRefLineThicknesses = (/1./)
    resL@gsnYRefLineColors  = (/"black"/)
    resL@gsnYRefLineDashPatterns = (/0/) 

    ;--------legend
    resL@pmLegendDisplayMode="Always"
    resL@xyExplicitLegendLabels=(/" NAJ"," NAWH"/)
    resL@lgLabelFontHeightF = font
    resL@lgPerimOn = False
    resL@pmLegendWidthF=0.05
    resL@pmLegendHeightF=0.055
    resL@pmLegendOrthogonalPosF=-0.65  ;;-0.58
    resL@pmLegendParallelPosF=0.85
    resL@tmXBLabelFontHeightF = font
    resL@tmYLLabelFontHeightF = font
    plotb = gsn_csm_xy(wks,xais,yyy,resL) 
draw(plotb)
frame(wks)
exit
  gplt    = new(2, graphic)
  gplt(0) = plot1
  gplt(1) = plotb
  pres                     =True
  pres@gsnPanelDebug      = True                   
  pres@gsnPanelXF         = (/ 0.07, 0.56/)
  pres@txString            =""
  ; pres@gsnPanelXWhiteSpacePercent = 1.2
  ; pres@gsnPanelYWhiteSpacePercent = 1.5
  pres@gsnPanelLabelBar    = False
  pres@gsnPanelRowSpec     = False
  pres@pmLabelBarWidthF    = 0.6
  pres@pmLabelBarHeightF   = 0.07
  pres@lbLabelFontHeightF  = 0.013
  ; ; pres@lbTitleString      ="units:  m^2*PVU/s"
  gsn_panel(wks,gplt,(/1,2/),pres)
end
