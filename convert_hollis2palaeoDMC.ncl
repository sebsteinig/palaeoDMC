; get path to working directory
dir = getenv("dir")

f_proxy_data_raw  = "full_deepMIP_temperature_compilation_1.5_subsampled.csv"

cut_off           = tofloat(getenv("CUT_OFF"))
ref_frame         = getenv("REF_FRAME") 

proxy_data_lines  = asciiread(dir+"/data/"+f_proxy_data_raw,-1,"string")

num_lines         = dimsizes(proxy_data_lines)-1
delim             = ","

proxy_data        = new((/num_lines,4/),"float")
time_slice        = str_get_field(proxy_data_lines(1:num_lines),2,delim)           ;--- time slice
realm             = str_get_field(proxy_data_lines(1:num_lines),3,delim)           ;--- SST or LAT measurements
proxy             = str_get_field(proxy_data_lines(1:num_lines),4,delim)           ;--- proxy technique
;load subsampling flags
flag_all          = str_get_field(proxy_data_lines(1:num_lines),21,delim)
flag_NoFrosty     = str_get_field(proxy_data_lines(1:num_lines),22,delim)
flag_NoTEX        = str_get_field(proxy_data_lines(1:num_lines),23,delim)
flag_NoMBT        = str_get_field(proxy_data_lines(1:num_lines),24,delim)
flag_NoMammal     = str_get_field(proxy_data_lines(1:num_lines),25,delim)
flag_NoSWPacific  = str_get_field(proxy_data_lines(1:num_lines),26,delim)

realm_flag  = where(realm .eq. "sst",1,0)

if (ref_frame .eq. "pmag") then
  proxy_data(:,0)   = tofloat(str_get_field(proxy_data_lines(1:num_lines),10,delim)) ;--- paleo-longitude
  proxy_data(:,1)   = tofloat(str_get_field(proxy_data_lines(1:num_lines),9,delim))  ;--- paleo-latitude
else if (ref_frame .eq. "mantle") then
  proxy_data(:,0)   = tofloat(str_get_field(proxy_data_lines(1:num_lines),8,delim)) ;--- paleo-longitude
  proxy_data(:,1)   = tofloat(str_get_field(proxy_data_lines(1:num_lines),7,delim))  ;--- paleo-latitude
end if
end if

proxy_data(:,2)   = tofloat(str_get_field(proxy_data_lines(1:num_lines),12,delim)) ;--- mean T (50 percentile of temperature estimates)
upper_bound       = tofloat(str_get_field(proxy_data_lines(1:num_lines),18,delim))
lower_bound       = tofloat(str_get_field(proxy_data_lines(1:num_lines),17,delim))
upper_bound       = where(upper_bound .eq. -999.9,upper_bound@_FillValue,upper_bound)
lower_bound       = where(lower_bound .eq. -999.9,lower_bound@_FillValue,lower_bound)

proxy_data(:,3)   = where(proxy .eq. "tex86" .or. proxy .eq. "mgca",(((upper_bound-proxy_data(:,2)) / 1.64 + (proxy_data(:,2)-lower_bound) / 1.64))/2.,((upper_bound-proxy_data(:,2)) + (proxy_data(:,2)-lower_bound)) / 2.)
sd_orig           = ((upper_bound-proxy_data(:,2)) + (proxy_data(:,2)-lower_bound)) / 2. ;save original SD definition for comparison

;--- 
;--- find sites with a reported error of 0 and set error to missing
proxy_data(:,3) = where(proxy_data(:,3) .eq. 0,proxy_data@_FillValue,proxy_data(:,3))
;--- find sites with a reported error lower than the cut-off (assumed as unrealistic) and set error to the chosen cut-off value
proxy_data(:,3) = where(proxy_data(:,3) .lt. cut_off,cut_off,proxy_data(:,3))
;--- find sites with no reported error and set error to average error of proxy method
do nn=0,num_lines-1
  if (ismissing(proxy_data(nn,3)) .eq. True) then
    ind_proxy := ind(proxy .eq. proxy(nn)) ;find other records of proxy method
    proxy_data(nn,3) = avg(proxy_data(ind_proxy,3)) ;compute average error of proxy method
  end if
end do

sd_new  = proxy_data(:,3)

;--- plot overview of SD/errors for different proxies
  proxy_labels        = get_unique_values(proxy)
  proxy_nums          = new(dimsizes(proxy_labels),"integer")
  proxy_nums_lp       = new(dimsizes(proxy_labels),"integer")
  proxy_nums_petm     = new(dimsizes(proxy_labels),"integer")
  proxy_nums_eeco     = new(dimsizes(proxy_labels),"integer")
 
  proxy_err_lp       = new((/dimsizes(proxy_labels),2,3/),"float")
  proxy_err_petm     = new((/dimsizes(proxy_labels),2,3/),"float")
  proxy_err_eeco     = new((/dimsizes(proxy_labels),2,3/),"float")

  do nn=0,dimsizes(proxy_labels)-1   
    proxy_nums(nn)  = num(ind(proxy .eq. proxy_labels(nn)))
  end do

  ip  = dim_pqsort(proxy_nums,-1)
  proxy_labels_sorted = proxy_labels(ip)
  proxy_labels_sorted(2)  = "mgca"
  proxy_labels_sorted(3)  = "del47"
  proxy_labels_sorted(4)  = "mbtcbt"
  proxy_labels_sorted(5)  = "lma"

  do nn=0,dimsizes(proxy_labels_sorted)-1   
    ind_lp              := ind(proxy .eq. proxy_labels_sorted(nn) .and. time_slice .eq. "lp")
    ind_petm            := ind(proxy .eq. proxy_labels_sorted(nn) .and. time_slice .eq. "petm")
    ind_eeco            := ind(proxy .eq. proxy_labels_sorted(nn) .and. time_slice .eq. "eeco")

    proxy_nums_lp(nn)    = num(.not. ismissing(ind_lp))
    proxy_nums_petm(nn)  = num(.not. ismissing(ind_petm))
    proxy_nums_eeco(nn)  = num(.not. ismissing(ind_eeco))

    if (proxy_nums_lp(nn) .gt. 0) then
      proxy_err_lp(nn,0,0)   = avg(sd_orig(ind_lp))
      proxy_err_lp(nn,0,1)   = min(sd_orig(ind_lp))
      proxy_err_lp(nn,0,2)   = max(sd_orig(ind_lp))
      proxy_err_lp(nn,1,0)   = avg(sd_new(ind_lp))
      proxy_err_lp(nn,1,1)   = min(sd_new(ind_lp))
      proxy_err_lp(nn,1,2)   = max(sd_new(ind_lp))
    end if
    if (proxy_nums_petm(nn) .gt. 0) then
      proxy_err_petm(nn,0,0)   = avg(sd_orig(ind_petm))
      proxy_err_petm(nn,0,1)   = min(sd_orig(ind_petm))
      proxy_err_petm(nn,0,2)   = max(sd_orig(ind_petm))
      proxy_err_petm(nn,1,0)   = avg(sd_new(ind_petm))
      proxy_err_petm(nn,1,1)   = min(sd_new(ind_petm))
      proxy_err_petm(nn,1,2)   = max(sd_new(ind_petm))
    end if
    if (proxy_nums_eeco(nn) .gt. 0) then
      proxy_err_eeco(nn,0,0)   = avg(sd_orig(ind_eeco))
      proxy_err_eeco(nn,0,1)   = min(sd_orig(ind_eeco))
      proxy_err_eeco(nn,0,2)   = max(sd_orig(ind_eeco))
      proxy_err_eeco(nn,1,0)   = avg(sd_new(ind_eeco))
      proxy_err_eeco(nn,1,1)   = min(sd_new(ind_eeco))
      proxy_err_eeco(nn,1,2)   = max(sd_new(ind_eeco))
    end if

  end do  
  proxy_nums_all =  proxy_nums_lp + proxy_nums_petm + proxy_nums_eeco
  
  wks = gsn_open_wks("pdf",dir+"/results/plots/DeepMIP_proxy_SD")

  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@tiXAxisFontHeightF  = 0.02
  res@tiYAxisFontHeightF  = 0.02
  res@trYMinF             = 0.5
  res@trYMaxF             = tofloat(dimsizes(proxy_labels_sorted))+.5
  res@xyMarkLineModes     = "Markers"
  res@xyMarkerSizeF       = 0.01
  res@xyMarkers           = 16.
  res@trYReverse          = True
  res@vpHeightF           = 0.9
  res@vpWidthF            = 0.3
  res@tmYLMode            = "Explicit"	
  res@tmYLValues          = ispan(1,dimsizes(proxy_labels_sorted),1)
  res@tmYLLabels          = proxy_labels_sorted
  res@tmYLLabelAngleF     = -45
  res@tmYLLabelJust       = "CenterCenter"
  res@tmXTOn              = True
  res@tmXTLabelsOn        = True
  res@tiXAxisSide         = "Top"

  pgres                  = True
  pgres@gsEdgesOn        = True
  pgres@gsLineThicknessF  = 3.
  pgres@tfPolyDrawOrder   = "PreDraw"

  pmres                  = True
  pmres@gsMarkerIndex    = 4
  pmres@gsMarkerSizeF    = 0.008
  pmres@gsMarkerColor    = "black"
  pmres@gsMarkerThicknessF  = 2.5
;  pmres@tfPolyDrawOrder   = "PostDraw"

  pmres_ind              = pmres
  pmres_ind@gsMarkerThicknessF  = 2.
  pmres_ind@gsMarkerOpacityF  = 0.6

  resP              = True
  resP@gsnFrame     = False
;  resP@gsnMaximize  = True
  resP@gsnPanelFigureStringsFontHeightF = 0.015
  resP@amJust         = "BottomRight"
;  resP@gsnPanelFigureStrings = (/"available data","original errors","adjusted errors"/)

  plot1 = new(3,graphic)
  plot2 = new(3,graphic)
  plot3 = new(3,graphic)
  mark1 = new(3,graphic)
  mark2 = new(3,graphic)
  mark3 = new(3,graphic)

  poly1 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  poly2 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  poly3 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  circ1 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  circ2 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  circ3 = new(dimsizes(proxy_labels_sorted)*3,graphic)

  shade1 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  shade2 = new(dimsizes(proxy_labels_sorted)*3,graphic)
  shade3 = new(dimsizes(proxy_labels_sorted)*3,graphic)

  width = 0.25

  res@trXMinF             = 0
  res@trXMaxF             = 20
  res@tiXAxisString       = "# locations"
  res@tmXBMode            = "Explicit"	
  res@tmXBValues          = (/0,4,8,12,16/)
  res@tmXBLabels          = tostring(res@tmXBValues)
  res@tmXBMinorValues     = ispan(0,20,1)

  txres               = True                            
  txres@txFontHeightF = 0.015             ; Set the font height
  txres@txAngleF      = -90
  txres@txJust        = "TopCenter"
  dummy = proxy_err_lp(:,0,0)
  dummy = -999.

  res@xyMarkerColor       = (/27,158,119/)/256.
  yspan                   = fspan(1-width, tofloat(dimsizes(proxy_labels_sorted))-width,dimsizes(proxy_labels))
  plot1(0) = gsn_csm_xy(wks,dummy,yspan,res)
    
  pgres@gsEdgesOn        = False
  pgres@gsFillColor     := "gray90"

  yspan                   = fspan(1-width, tofloat(dimsizes(proxy_labels_sorted))-width,dimsizes(proxy_labels))
  do pp=0,dimsizes(proxy_labels)-1
    if (pp/2 .ne. pp/2.) then
      shade1(pp)  = gsn_add_polygon(wks,plot1(0),(/0,100,100,0,0/),(/yspan(pp)+3*width,yspan(pp)+3*width,yspan(pp)-width,yspan(pp)-width,yspan(pp)+3*width/),pgres)
    end if
  end do 

  pgres@gsEdgesOn         = True
  pgres@gsFillColor      := (/27,158,119/)/256.
  do pp=0,dimsizes(proxy_labels)-1
    poly1(pp)  = gsn_add_polygon(wks,plot1(0),(/0,proxy_nums_lp(pp),proxy_nums_lp(pp),0,0/),(/yspan(pp)+width/2.,yspan(pp)+width/2.,yspan(pp)-width/2.,yspan(pp)-width/2.,yspan(pp)+width/2./),pgres)
  end do

  pgres@gsFillColor      = (/217,95,2/)/256.
  yspan                   = fspan(1., tofloat(dimsizes(proxy_labels_sorted)),dimsizes(proxy_labels_sorted))
  do pp=0,dimsizes(proxy_labels)-1
    poly2(pp)  = gsn_add_polygon(wks,plot1(0),(/0,proxy_nums_petm(pp),proxy_nums_petm(pp),0,0/),(/yspan(pp)+width/2.,yspan(pp)+width/2.,yspan(pp)-width/2.,yspan(pp)-width/2.,yspan(pp)+width/2./),pgres)
  end do

  pgres@gsFillColor      = (/117,112,179/)/256.
  yspan                   = fspan(1+width, tofloat(dimsizes(proxy_labels_sorted))+width,dimsizes(proxy_labels_sorted))
  do pp=0,dimsizes(proxy_labels)-1
    poly3(pp)  = gsn_add_polygon(wks,plot1(0),(/0,proxy_nums_eeco(pp),proxy_nums_eeco(pp),0,0/),(/yspan(pp)+width/2.,yspan(pp)+width/2.,yspan(pp)-width/2.,yspan(pp)-width/2.,yspan(pp)+width/2./),pgres)
  end do

  pgres@gsFillColor      := "white"
  pgres@gsEdgesOn         = True
  legbox  = gsn_add_polygon(wks,plot1(0),(/13.3,18.7,18.7,13.3,13.3/),(/12.7,12.7,16.8,16.8,12.7/)+2.,pgres)

  pgres@gsEdgesOn         = False
  pgres@gsFillColor      := (/27,158,119/)/256.
  leg1  = gsn_add_polygon(wks,plot1(0),(/17,18,18,17,17/),(/13,13,15,15,13/)+2.,pgres)
  pgres@gsFillColor      = (/217,95,2/)/256.
  leg2  = gsn_add_polygon(wks,plot1(0),(/15.5,16.5,16.5,15.5,15.5/),(/13,13,15,15,13/)+2.,pgres)
  pgres@gsFillColor      = (/117,112,179/)/256.
  leg3  = gsn_add_polygon(wks,plot1(0),(/14,15,15,14,14/),(/13,13,15,15,13/)+2.,pgres)
  txt1  = gsn_add_text(wks,plot1(0),(/"LP","PETM","EECO"/),(/17.5,16.,14.5/),(/15.3,15.3,15.3/)+2.,txres)
  res@trXMinF             = 0

  pgres@gsEdgesOn        = False
  pgres@gsFillColor     := "gray90"

  res@tmYLLabelsOn       = False

  res1                     = True
  res2                     = True
  res1@gsnMaximize         = True
 
  do nn=0,1
    if (nn .eq. 0) then
        res@tiXAxisString       = "DeepMIP error"
        res@trXMaxF            := 10.
        res@tmXBValues         := (/0,2,4,6,8/)
        res@tmXBLabels         := tostring(res@tmXBValues)
        res@tmXBMinorValues    := ispan(0,100,5)/10.

    else
        res@tiXAxisString       = "adjusted SD [~S~o~N~C]"
        res@trXMinF            := 0.
        res@trXMaxF            := 8.
    end if

    res@xyMarkerColor       = (/27,158,119/)/256.
    yspan                   = fspan(1-width, tofloat(dimsizes(proxy_labels_sorted))-width,dimsizes(proxy_labels_sorted))
    plot1(nn+1) = gsn_csm_xy(wks,proxy_err_lp(:,nn,0),yspan,res)
  
    do pp=0,dimsizes(proxy_labels)-1
      if (pp/2 .ne. pp/2.) then
        shade2((nn+1)*dimsizes(proxy_labels) + pp)  = gsn_add_polygon(wks,plot1(nn+1),(/0,100,100,0,0/),(/yspan(pp)+3*width,yspan(pp)+3*width,yspan(pp)-width,yspan(pp)-width,yspan(pp)+3*width/),pgres)
      end if
    end do 

    pmres_ind@gsMarkerColor      := (/27,158,119/)/256.
    do pp=0,dimsizes(proxy_labels)-1
      indices :=  ind(proxy .eq. proxy_labels_sorted(pp) .and. time_slice .eq. "lp")
      if (ismissing(indices(0)) .eq. False) then
        yplot := new(dimsizes(indices),"float")
        yplot = yspan(pp)
        if (nn .eq. 0) then
          circ1((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot1(nn+1),sd_orig(indices),yplot,pmres_ind)
        else
          circ1((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot1(nn+1),sd_new(indices),yplot,pmres_ind)
        end if
      end if
    end do



    res@xyMarkerColor       = (/217,95,2/)/256.
    yspan                   = fspan(1., tofloat(dimsizes(proxy_labels_sorted)),dimsizes(proxy_labels_sorted))
    plot2(nn+1) = gsn_csm_xy(wks,proxy_err_petm(:,nn,0),yspan,res)

    pmres_ind@gsMarkerColor      = (/217,95,2/)/256.
    do pp=0,dimsizes(proxy_labels)-1
      indices :=  ind(proxy .eq. proxy_labels_sorted(pp) .and. time_slice .eq. "petm")
      if (ismissing(indices(0)) .eq. False) then
        yplot := new(dimsizes(indices),"float")
        yplot = yspan(pp)
        if (nn .eq. 0) then
          circ2((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot2(nn+1),sd_orig(indices),yplot,pmres_ind)
        else
          circ2((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot2(nn+1),sd_new(indices),yplot,pmres_ind)
        end if
      end if
    end do


    res@xyMarkerColor       = (/117,112,179/)/256.
    yspan                   = fspan(1+width, tofloat(dimsizes(proxy_labels_sorted))+width,dimsizes(proxy_labels_sorted))
    plot3(nn+1) = gsn_csm_xy(wks,proxy_err_eeco(:,nn,0),yspan,res)

    pmres_ind@gsMarkerColor      = (/117,112,179/)/256.
    do pp=0,dimsizes(proxy_labels)-1
      indices :=  ind(proxy .eq. proxy_labels_sorted(pp) .and. time_slice .eq. "eeco")
      if (ismissing(indices(0)) .eq. False) then
        yplot := new(dimsizes(indices),"float")
        yplot = yspan(pp)
        if (nn .eq. 0) then
          circ3((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot3(nn+1),sd_orig(indices),yplot,pmres_ind)
        else
          circ3((nn+1)*dimsizes(proxy_labels) + pp) = gsn_add_polymarker(wks,plot3(nn+1),sd_new(indices),yplot,pmres_ind)
        end if
      end if
    end do

    yspan                   = fspan(1-width, tofloat(dimsizes(proxy_labels_sorted))-width,dimsizes(proxy_labels_sorted))
    mark1(nn+1) = gsn_add_polymarker(wks,plot3(nn+1),proxy_err_lp(:,nn,0),yspan,pmres)
    yspan                   = fspan(1., tofloat(dimsizes(proxy_labels_sorted)),dimsizes(proxy_labels_sorted))
    mark2(nn+1) = gsn_add_polymarker(wks,plot3(nn+1),proxy_err_petm(:,nn,0),yspan,pmres)
    yspan                   = fspan(1+width, tofloat(dimsizes(proxy_labels_sorted))+width,dimsizes(proxy_labels_sorted))
    mark3(nn+1) = gsn_add_polymarker(wks,plot3(nn+1),proxy_err_eeco(:,nn,0),yspan,pmres)

    overlay(plot1(nn+1),plot2(nn+1))
    overlay(plot1(nn+1),plot3(nn+1))
  end do

  anno1 = gsn_attach_plots(plot1(1),plot1(2),res1,res2)
  anno2 = gsn_attach_plots(plot1(0),plot1(1),res1,res2)

  draw(plot1(0))


;--- convert absolute proxy temperatures to anomaly with respect to present-day reanalysis
  ;--- load reanalysis data
  f_HadISST   = addfile(dir+"/reanalysis_reference/HadISST_sst_1979-2018.timmean.zonmean.nc","r")
  f_ERA5      = addfile(dir+"/reanalysis_reference/ERA5_t2m_1979-2018.timmean.zonmean.nc","r")
  SST_zm      = f_HadISST->sst(0,:,0)
  GAT_zm      = short2flt(f_ERA5->t2m(0,:,0))
  GAT_zm      = GAT_zm-273.15

  ;--- calculate respective anomaly
 do nn=0,num_lines-1
   if (realm(nn) .eq. "sst") then ;--- anomaly to zonal mean SST
     proxy_data(nn,2) = proxy_data(nn,2) - SST_zm({lat|proxy_data(nn,1)})
   else if (realm(nn) .eq. "lat") then ;--- anomaly to zonal mean GAT
     proxy_data(nn,2) = proxy_data(nn,2) - GAT_zm({lat|proxy_data(nn,1)})
   end if
   end if
 end do

;--- write selected data to individual ASCII files
time_slice_names  = (/"petm","lp","eeco"/)
do tt=0,dimsizes(time_slice_names)-1
  system("mkdir -p " + dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt))
  system("rm -f " + dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_*.text")
  do rr=0,2
    if (rr .eq. 0) then
      filename_ctrl           = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_ctrl_all.txt"
      filename_WithFrosty     = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_WithFrosty_all.txt"
      filename_NoTEX          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoTEX_all.txt"
      filename_NoMBT          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMBT_all.txt"
      filename_NoMammal       = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMammal_all.txt"
;      filename_NoSWPacific    = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoSWPacific.txt"
    else if (rr .eq. 1) then
      filename_ctrl            = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_ctrl_sst.txt"
      filename_WithFrosty       = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_WithFrosty_sst.txt"
      filename_NoTEX          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoTEX_sst.txt"
      filename_NoMBT          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMBT_sst.txt"
      filename_NoMammal       = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMammal_sst.txt"
;      filename_NoSWPacific    = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoSWPacific_sst.txt"
    else if (rr .eq. 2) then
      filename_ctrl            = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_ctrl_lat.txt"
      filename_WithFrosty       = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_WithFrosty_lat.txt"
      filename_NoTEX          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoTEX_lat.txt"
      filename_NoMBT          = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMBT_lat.txt"
      filename_NoMammal       = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoMammal_lat.txt"
;      filename_NoSWPacific    = dir+"/Observation_Data/DeepMIP_"+time_slice_names(tt)+"/DeepMIP_"+ref_frame+"_"+time_slice_names(tt)+"_NoSWPacific_lat.txt"
    end if
    end if
    end if
    system("rm -f " + filename_ctrl)
    system("rm -f " + filename_WithFrosty)
    system("rm -f " + filename_NoTEX)
    system("rm -f " + filename_NoMBT)
    system("rm -f " + filename_NoMammal)
;    system("rm -f " + filename_NoSWPacific)

    header = "x y z std sst"
    print("creating file "+filename_ctrl)
    print("creating file "+filename_WithFrosty)
    print("creating file "+filename_NoTEX)
    print("creating file "+filename_NoMBT)
    print("creating file "+filename_NoMammal)
;    print("creating file "+filename_NoSWPacific)

    write_table(filename_ctrl, "w", [/header/], "%s")    ; Use "w" to create file
    write_table(filename_WithFrosty, "w", [/header/], "%s")    ; Use "w" to create file
    write_table(filename_NoTEX, "w", [/header/], "%s")    ; Use "w" to create file
    write_table(filename_NoMBT, "w", [/header/], "%s")    ; Use "w" to create file
    write_table(filename_NoMammal, "w", [/header/], "%s")    ; Use "w" to create file
;    write_table(filename_NoSWPacific, "w", [/header/], "%s")    ; Use "w" to create file

    count = 0
    do nn=0,num_lines-1
      if (time_slice(nn) .eq. time_slice_names(tt)) then
        if (rr .eq. 1 .and. realm(nn) .ne. "sst") then
          continue
        else if (rr .eq. 2 .and. realm(nn) .ne. "lat") then
          continue
        end if
        end if
        a :=  [/proxy_data(nn,0),proxy_data(nn,1),proxy_data(nn,2),proxy_data(nn,3),realm_flag(nn)/]
        ;a :=  [/proxy_data(nn,0),proxy_data(nn,1),proxy_data(nn,2),proxy_data(nn,3)/]
        format  = "%5.2f%5.2f%5.2f%5.2f%i"
        ;format  = "%5.2f%5.2f%5.2f%5.2f"

        if (flag_all(nn) .eq. 1) then
          write_table(filename_WithFrosty, "a", a, format)
        end if
        if (flag_NoFrosty(nn) .eq. 1) then
          write_table(filename_ctrl, "a", a, format)
        end if
        if (flag_NoTEX(nn) .eq. 1 .and. flag_NoFrosty(nn) .eq. 1) then
          write_table(filename_NoTEX, "a", a, format)
        end if
        if (flag_NoMBT(nn) .eq. 1 .and. flag_NoFrosty(nn) .eq. 1) then
          write_table(filename_NoMBT, "a", a, format)
        end if
        if (flag_NoMammal(nn) .eq. 1 .and. flag_NoFrosty(nn) .eq. 1) then
          write_table(filename_NoMammal, "a", a, format)
        end if
 ;       if (flag_NoSWPacific(nn) .eq. 1 .and. flag_NoFrosty(nn) .eq. 1) then
 ;         write_table(filename_NoSWPacific, "a", a, format)
 ;       end if
      end if
    end do
  end do
end do
