load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
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
ied = YrEnd   - 1891

NAWH_SSTGbased = fout->NAWH_SSTGbased(ist:ied) 

fil2 = "NAJ.lat.JJA.v2.20V3.1806-2015.nc"
dir = "/home/users/shengc01/work11/scripts10/data/"
fout = addfile(dir + fil2, "r")
jet_lat1 = new(nt, float, -999999.)
jet_lat1 = -999999.

ist = YrStart - 1806
ied = 2015 - 1806
jet_lat1(0:115) = fout->jet_lat1(ist:ied)


ist0  = 1900 - 1900
iend0 = 1940 - 1900

ist1  = 1940 - 1900
iend1 = 1980 - 1900

ist2  = 1980 - 1900
iend2 = 2019 - 1900

NAWHI = new((/3,41/), float,-999999.)
NAWHI = -999999.
NAWHI(0,:) = NAWH_SSTGbased(ist0:iend0) * -1.
NAWHI(1,:) = NAWH_SSTGbased(ist1:iend1) * -1.
NAWHI(2,0:39) = NAWH_SSTGbased(ist2:iend2) * -1.
NAWHI  = dim_standardize(NAWHI, 1)


ist2  = 1980 - 1900
iend2 = 2015 - 1900

NAJI = new((/3,41/), float,-999999.)
NAJI = -999999.
NAJI(0,:) = jet_lat1(ist0:iend0) 
NAJI(1,:) = jet_lat1(ist1:iend1) 
NAJI(2,0:35) = jet_lat1(ist2:iend2) 
NAJI  = dim_standardize(NAJI, 1)

y_P1 = new((/2,41/), float,-999999.)
y_P2 = new((/2,41/), float,-999999.)
y_P3 = new((/2,41/), float,-999999.)
y_P1(0,:) = NAWHI(0,:)
y_P1(1,:) = NAJI(0,:)

y_P2(0,:) = NAWHI(1,:)
y_P2(1,:) = NAJI(1,:)

y_P3(0,:) = NAWHI(2,:)
y_P3(1,:) = NAJI(2,:)

; ;**********regression*************************
x  = ispan(0,40,1)*1.
RC := regCoef_n(x,NAWHI,0,1)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob := student_t(tval, df)
a  = RC
b  = RC@yintercept
NAWH_trd = NAWHI
NAWH_trd(0,:) = a(0) * x + b(0)
NAWH_trd(1,:) = a(1) * x + b(1)
NAWH_trd(2,:) = a(2) * x + b(2)


; ;**********regression*************************
x  = ispan(0,40,1)*1.
RC := regCoef_n(x,NAJI,0,1)
df := onedtond(RC@nptxy,dimsizes(RC)) - 2 
tval := onedtond(RC@tval,dimsizes(RC))
prob := student_t(tval, df)
a  = RC
b  = RC@yintercept
NAJ_trd = NAJI
NAJ_trd(0,:) = a(0) * x + b(0)
NAJ_trd(1,:) = a(1) * x + b(1)
NAJ_trd(2,:) = a(2) * x + b(2)

ytrd_P1 = new((/2,41/), float,-999999.)
ytrd_P2 = new((/2,41/), float,-999999.)
ytrd_P3 = new((/2,41/), float,-999999.)
ytrd_P1(0,:) = NAWH_trd(0,:)
ytrd_P1(1,:) = NAJ_trd(0,:)

ytrd_P2(0,:) = NAWH_trd(1,:)
ytrd_P2(1,:) = NAJ_trd(1,:)

ytrd_P3(0,:) = NAWH_trd(2,:)
ytrd_P3(1,:) = NAJ_trd(2,:)

;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic10/Fig.E5.NAWH."+YrStart+"."+YrEnd)

;; for line=========
xais = ispan(0,40,1)*1.
yr   = ispan(1900, 1940, 1)
font = 0.03
;; for first plot==========
resL = True 
resL@gsnMaximize                   =False
resL@gsnDraw                       =False
resL@gsnFrame                      =False
resL@vpHeightF= 0.4                    
resL@vpWidthF = 0.7  
resL@gsnLeftString = ""
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
resL@xyLineThicknesses  = (/3.,3./)
resL@xyLineColors = (/"gray60","red"/)
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
resL@tmXBValues        =  (/0.,20.,40./)  
resL@tmXBLabels        = (/1900,1920,1940/) 
resL@tmXBLabelFontHeightF = font
resl = resL

resL@pmLegendDisplayMode="Always"
resL@lgPerimOn          = False
resL@xyExplicitLegendLabels= (/" NAWH"," NAJ"/) 
resL@pmLegendWidthF=0.1
resL@pmLegendHeightF=0.12
resL@pmLegendOrthogonalPosF=-1.3 
resL@pmLegendParallelPosF=0.6
resL@lgLabelFontHeightF = 0.019  
resL@lgLabelFontColor   = "black"
resL@gsnLeftString = "(a)1900-1940"
resL@gsnRightString = "";;"JJA"
plot1 = gsn_csm_xy (wks,xais,y_P1,resL)

resl@gsnLeftString = ""
resl@gsnRightString = ""
resl@xyLineThicknesses = (/2.,2./)
resl@xyDashPattern = 1


plot1_1 = gsn_csm_xy (wks,xais,ytrd_P1,resl)
overlay(plot1, plot1_1)
resL@gsnLeftString = "(b)1940-1980"
resL@tmXBLabels        = (/1940,1960,1980/) 
plot2 = gsn_csm_xy (wks,xais,y_P2,resL)
plot2_1 = gsn_csm_xy (wks,xais,ytrd_P2,resl)
overlay(plot2, plot2_1)
resL@gsnLeftString = "(c)1980-2019"
resL@trXMaxF           = max(xais-1) 
resL@tmXBValues        =  (/0.,20.,39./)  
resL@tmXBLabels        = (/1980,2000,2019/) 
plot3 = gsn_csm_xy (wks,xais,y_P3,resL)
plot3_1 = gsn_csm_xy (wks,xais,ytrd_P3,resl)
overlay(plot3, plot3_1)

gplt   = new(3, graphic)
gplt(0) = plot1 
gplt(1) = plot2
gplt(2) = plot3
pres                     =True
gsn_panel(wks,gplt,(/1,3/),pres)

end
