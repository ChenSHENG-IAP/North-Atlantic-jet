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
print(nmem)
; ;**********sort*************************
x  = (TIME-1980)*1.
RC_NAJ = regCoef_n(x,NAJ_lat,0,0)
RC_NAJ!0 = "member"
RC_NAJ&member = member



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

dirout = "/home/users/shengc01/work11/scripts4/data/"
filout = "AMO.index.JJA.MPI.LE.1860-2019.v3.nc"
fout = addfile(dirout + filout, "r")
time   := fout->time
YYYY   := cd_calendar(time,-1)/100    ; entire file
iYYYY  := ind(YYYY.ge.YrStart .and. YYYY.le.YrEnd)
AMOI = fout->AMO3(iYYYY,:)
AMOI = dim_rmvmean_n(AMOI,1)
RC_AMOI = regCoef_n(x,AMOI,0,0)
RC_AMOI!0 = "member"
RC_AMOI&member = member





NAJ_trd = dim_standardize(RC_NAJ, 1)
IPO_trd = dim_standardize(RC_IPOI, 1)
AMO_trd = dim_standardize(RC_AMOI, 1)

R_NAJ_IPO = escorc(NAJ_trd, IPO_trd)
P_NAJ_IPO = rtest(R_NAJ_IPO, nmem, 0)

R_NAJ_AMO = escorc(NAJ_trd, AMO_trd)
P_NAJ_AMO = rtest(R_NAJ_AMO, nmem, 0)



print("P-value IN CC(NAJ,IPO):"+P_NAJ_IPO)
print("P-value IN CC(NAJ,AMO):"+P_NAJ_AMO)

x1 = fspan(-2., 2, 40)
y1 = R_NAJ_IPO*x1 
y2 = R_NAJ_AMO*x1 


wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic8/Fig.E8.scatter.plt.NAJ.AMO.IPO.MPI_GE" )      
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
  res@vpXF              = 0.065
  res@vpYF              = 0.85
  res@vpWidthF          = 0.32
  res@vpHeightF         = 0.32

;---Resources for X/Y axes.

  res@trYMinF            =  -3.               ; min value on y-axis
  res@trYMaxF            =  3.               ; max value on y-axis
  res@trXMinF            =  -3.                  ; min value on x-axis
  res@trXMaxF            =  3.                  ; max value on x-axis

res@tmYLFormat        ="0@;*.1f"       
res@tmXBFormat        ="0@;*.1f"       
res@tmYLLabelFontHeightF = font
res@tmXBLabelFontHeightF = font


marker = 16
;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.005
  res@xyMarkerThicknessF = 0.5
  res@xyMonoMarkerColor  = True          ; Allow different colors for markers
  res@xyMarkerColors    = "black"           ; It's okay to list more than
                                          ; you need here
  res@xyMonoMarker       = True                                        
  res@xyMarker          = marker

;---Tickmark resources
  res@tmXTOn           = False
  res@tmYROn           = False
  res@tmXBMode         = "Manual"
  res@tmYLMode         = "Manual"
  res@tmXBTickSpacingF = 1.
  res@tmYLTickSpacingF = 1.

;---Title resources
  res@tiXAxisString      = "";;"AMO trend"
  res@tiYAxisString      = "";;"NAJ trend"
  res@tiMainFontHeightF  = 0.03
  res@gsnXRefLine        = 0.0
  res@gsnYRefLine        = 0.0
  res@gsnXRefLineThicknessF = 1.5
  res@gsnYRefLineThicknessF = 1.5
;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------

  res@gsnLeftString = "";;"(a)trend of IPO&NAJ"
  res@gsnRightString = "";;"corr="+tostring_with_format(R_NAJ_IPO, "%3.2f")
  plot1 = gsn_csm_xy (wks,transpose((/IPO_trd,IPO_trd/)),transpose((/NAJ_trd,NAJ_trd/)),res)

  res@gsnLeftString = "";;"(b)trend of AMO&NAJ"
  res@gsnRightString = "";;"corr="+tostring_with_format(R_NAJ_AMO, "%3.2f")
  plot2 = gsn_csm_xy (wks,transpose((/AMO_trd,AMO_trd/)),transpose((/NAJ_trd,NAJ_trd/)),res)

; for y=kx
delete([/res@xyMarkLineMode,res@gsnLeftString,res@gsnRightString/])
res@xyLineThicknessF  = 5.
res@xyLineColor = "gray"
res@xyDashPattern=0
res@xyCurveDrawOrder = "PreDraw"
plot1_1 = gsn_csm_xy (wks,x1,y1,res)
overlay(plot1, plot1_1)

plot2_2 = gsn_csm_xy (wks,x1,y2,res)
overlay(plot2, plot2_2)

  gplt    = new(2, graphic)
  gplt(0) = plot1
  gplt(1) = plot2
  pres                     =True
  pres@gsnPanelDebug      = True                   
  ; pres@gsnPanelXF         = (/ 0.07, 0.56/)
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