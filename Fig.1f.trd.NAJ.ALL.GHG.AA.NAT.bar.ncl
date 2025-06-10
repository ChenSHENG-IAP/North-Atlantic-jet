load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;;==========================
TotalTime = systemfunc("date")
YrStart = 1980
YrEnd   = 2014
lat1 = 0.
lat2 = 90.
lon1 = 0.
lon2 = 360.
LEV  = 25000
lev1 = 100000
lev2 = 7000  
TIME  = ispan(YrStart, YrEnd, 1)
nt    = dimsizes(TIME)
model = (/"BCC-CSM2-MR","CanESM5","CNRM-CM6-1","HadGEM3-GC31-LL","IPSL-CM6A-LR","MRI-ESM2-0"/)
dir = "/home/users/shengc01/work11/scripts3/data/"
nmodel = dimsizes(model)

print(nmodel)
ALL_baro = new(200, double)
GHG_baro = new(200, double)
AA_baro  = new(200, double)
NAT_baro = new(200, double)
ist  = 0
do i = 0, nmodel-1
     model_str = model(i)
     f = addfile(dir+"NAJ.lat.JJA."+model_str+".historical.45W-0.v2.nc", "r")
     BCC_historical := f->NAJ_lat_trd 
     nm = dimsizes(BCC_historical)

     ALL_baro(ist:ist+nm-1) = BCC_historical
     ist = ist+nm
     ALL_num = ist 
end do 

ist  = 0
do i = 0, nmodel-1
     model_str = model(i)
     f = addfile(dir+"NAJ.lat.JJA."+model_str+".hist-GHG.45W-0.v2.nc", "r")
     BCC_GHG := f->NAJ_lat_trd 
     nm = dimsizes(BCC_GHG)

     GHG_baro(ist:ist+nm-1) = BCC_GHG
     ist = ist+nm
     GHG_num = ist 
end do 

ist  = 0
do i = 0, nmodel-1
     model_str = model(i)
     f = addfile(dir+"NAJ.lat.JJA."+model_str+".hist-aer.45W-0.v2.nc", "r")
     BCC_aer := f->NAJ_lat_trd 
     nm = dimsizes(BCC_aer)

     AA_baro(ist:ist+nm-1) = BCC_aer
     ist = ist+nm
     AA_num = ist 
end do 

ist  = 0
do i = 0, nmodel-1
     model_str = model(i)
     f = addfile(dir+"NAJ.lat.JJA."+model_str+".hist-nat.45W-0.v2.nc", "r")
     BCC_nat := f->NAJ_lat_trd 
     nm = dimsizes(BCC_nat)

     NAT_baro(ist:ist+nm-1) = BCC_nat
     ist = ist+nm
     NAT_num = ist 
end do 

ALL_bar = ALL_baro(0:ALL_num-1) * 10.
GHG_bar = GHG_baro(0:GHG_num-1) * 10.
AA_bar  = AA_baro(0:AA_num-1)   * 10.
NAT_bar = NAT_baro(0:NAT_num-1) * 10.



statb = stat_dispersion(ALL_bar, True)
ALL_bot = statb(3)  ;; lower 1/6
ALL_bot1= statb(6)  ;; lower 1/4
ALL_mean = avg(ALL_bar);;statb(8) ;; median         avg(ALL_bar)
ALL_top1 = statb(10) ;; upper 1/4
ALL_top = statb(13)  ;; upper 1/6

statb = stat_dispersion(GHG_bar, True)
GHG_bot = statb(3)  ;; lower 1/6
GHG_bot1= statb(6)  ;; lower 1/4
GHG_mean = avg(GHG_bar) ;;statb(8) ;; median         avg(GHG_bar)
GHG_top1 = statb(10) ;; upper 1/4
GHG_top = statb(13)  ;; upper 1/6

statb = stat_dispersion(AA_bar, True)
AA_bot = statb(3)  ;; lower 1/6
AA_bot1= statb(6)  ;; lower 1/4
AA_mean = avg(AA_bar) ;;statb(8) ;; median         avg(AA_bar)
AA_top1 = statb(10) ;; upper 1/4
AA_top = statb(13)  ;; upper 1/6


statb = stat_dispersion(NAT_bar, True)
NAT_bot = statb(3)  ;; lower 1/6
NAT_bot1= statb(6)  ;; lower 1/4
NAT_mean = avg(NAT_bar) ;;statb(8) ;; median         avg(NAT_bar)
NAT_top1 = statb(10) ;; upper 1/4
NAT_top = statb(13)  ;; upper 1/6



OBS = (1.06 + 1.0 + 1.08) / 3. * -1.
yval = new((/5,5/),"double",-999.)

yval(0,0) = 0.
yval(0,1) = 0.
yval(0,2) = OBS
yval(0,3) = 0.
yval(0,4) = 0.

yval(1,0) = ALL_bot
yval(1,1) = ALL_bot1
yval(1,2) = ALL_mean
yval(1,3) = ALL_top1
yval(1,4) = ALL_top

yval(2,0) = GHG_bot
yval(2,1) = GHG_bot1
yval(2,2) = GHG_mean
yval(2,3) = GHG_top1
yval(2,4) = GHG_top

yval(3,0) = AA_bot
yval(3,1) = AA_bot1
yval(3,2) = AA_mean
yval(3,3) = AA_top1
yval(3,4) = AA_top

yval(4,0) = NAT_bot
yval(4,1) = NAT_bot1
yval(4,2) = NAT_mean
yval(4,3) = NAT_top1
yval(4,4) = NAT_top



;;======================
wks   = gsn_open_wks ("eps", "/home/users/shengc01/work11/pic3/Fig1.trd.NAJ.ALL.GHG.AA.NAT.bar")

;====================================================================================== 

     sres = True
     sres@vpWidthF = 0.8
     sres@vpHeightF = 0.4
     sres@vpXF = .15
     sres@trXMinF = 0.4
     sres@trXMaxF = 5.6
     sres@trYMinF = -1.5

     sres@trYMaxF = 1.5
     sres@gsnYRefLine = 0.
     sres@gsnYRefLineThicknesses = 1.
     sres@gsnDraw = False
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     barwidth                   = 0.2
     sres@gsnXYBarChartBarWidth = barwidth           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/1,2,3,4,5/)
     sres@tmXBLabels = (/"OBS","ALL","GHG","AA","NAT"/)
     fht                       = 0.023
     sres@tmXBLabelFontHeightF = fht
     sres@tmXTLabelFontHeightF = fht
     sres@tmYLLabelFontHeightF = fht
     sres@tiMainFontHeightF = 0.025
     ; sres@tiMainFont = "helvetica"
     sres@tiMainString = ""
     sres@gsnRightString = ""
     sres@tiYAxisString = ""
     sres@gsnLeftString = "";; "trd of NAJ lat"
     ;sres@tmXBLabelAngleF   = 20 
     sres@tmXTOn = False
     sres@tmYROn = False   
     sres@xyLineColor = "gray15"
     sres@xyLineThicknesses = 1.
     sres@gsnXYBarChartColors2 = (/"black","green4","red","cornflowerblue","purple"/)    
     xx = ispan(1,5,1) * 1.
     plot1 = gsn_csm_xy(wks,xx,yval(:,2),sres)     ; draw each time series  his

;;11.erro bar 
     lnres = True
     lnres@gsLineColor = "gray75"
     lnres@gsLineThicknessF = 1.5
     lnres@gsLineDashPattern = 0.
     lnres@tfPolyDrawOrder   = "PreDraw"
     lnres@gsLineOpacityF   = (/0.7,0.7,0.7,0.7,0.7/)
color = sres@gsnXYBarChartColors2

     dum_ver = new(5, graphic)
     dum_up = new(5, graphic)
     dum_dw = new(5, graphic)
do iexp = 0,4
         yyy = (/yval(iexp,0),yval(iexp,4)/)   ;;65  (/40.,50,60/)
         xxx = iexp + 1.
         lnres@gsLineColor = color(iexp)
         dum_ver(iexp) = gsn_add_polyline(wks,plot1, (/xxx,xxx/), (/yval(iexp,0),yval(iexp,4)/) , lnres)
         dum_up(iexp) = gsn_add_polyline(wks, plot1, (/xxx-0.3*barwidth,xxx+0.3*barwidth/), (/yval(iexp,0),yval(iexp,0)/), lnres)
         dum_dw(iexp) = gsn_add_polyline(wks, plot1, (/xxx-0.3*barwidth,xxx+0.3*barwidth/), (/yval(iexp,4),yval(iexp,4)/), lnres)
end do 


  txres  = True
  txres@txFontHeightF = 0.017
  txres@txFontThicknessF = 0.03
  txres@txFontColor = "black"
  Yposition   = -1.4 
  ypct_string =  tostring_with_format(yval(:,2), "%3.2f") +"~S~o~N~ 10a~S~-1~N~"       ;;tostring(sprintf("%3.2f", yval(:,2)))
  xaixs  = ispan(1, 5, 1)  
  text  = new(5, graphic)
  do i = 0, 4
     txres@txFontColor = color(i)
     text(i)  = gsn_add_text(wks, plot1, ypct_string(i), xaixs(i), Yposition, txres)
  end do
  
draw(plot1)   
     
frame(wks)  
end