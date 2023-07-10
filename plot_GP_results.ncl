dir = getenv("dir")
time_slices = (/"lp","petm","eeco"/)
areas       = (/"_sst","_sst","_lat"/)
;experiments = (/"ctrl","WithFrosty","NoTEX","NoMBT","NoMammal","NoSWPacific"/)
;labels      = (/"Default","WithFrosty","NoTEX","NoMBT","NoMammal","NoSWPacific"/)
experiments = (/"ctrl","WithFrosty","NoTEX","NoMBT","NoMammal"/)
labels      = (/"Default","WithFrosty","NoTEX","NoMBT","NoMammal"/)

plot_time_slice_maps      = 1
plot_sst_lat_maps         = 0
plot_subsampling_maps     = 0
calc_gm                   = 0
plot_gm                   = 0
plot_zm                   = 0

f_ERA5      = addfile(dir+"/reanalysis_reference/ERA5_t2m_1979-2018.timmean.zonmean.1deg.nc","r")
OBS_zm      = short2flt(f_ERA5->t2m(0,:,0))
OBS_zm      = OBS_zm-273.15

;f_mask      = addfile(dir+"/Model_data/no_model/mask.nc","r")
;f_mask      = addfile(dir+"/mask_36x18.nc","r")
; f_mask      = addfile(dir+"/mask_144x72.nc","r")
f_mask      = addfile(dir+"/mask_72x36.nc","r")
ls_mask     = f_mask->mask
lat         = f_mask->lat
lon         = f_mask->lon

f_mask_outline  = addfile(dir+"/mask.nc","r")
ls_mask_outline  = f_mask_outline->mask
paleo_outline (ls_mask_outline,ls_mask_outline&lat,ls_mask_outline&lon,1,"herold_eocene")

cut_off     = getenv("CUT_OFF")    ; return string "1901"
ref_frame   = getenv("REF_FRAME") 

GP_mean            = new((/dimsizes(time_slices),dimsizes(areas),dimsizes(experiments),dimsizes(lat),dimsizes(lon)/),"float")
GP_SD              = new((/dimsizes(time_slices),dimsizes(areas),dimsizes(experiments),dimsizes(lat),dimsizes(lon)/),"float")
GP_mean_sample     = new((/dimsizes(time_slices),dimsizes(areas),dimsizes(experiments),20,dimsizes(lat),dimsizes(lon)/),"float")

do tt=0,dimsizes(time_slices)-1
  do rr=0,dimsizes(areas)-1
    do nn=0,dimsizes(experiments)-1
      f                  = addfile(dir+"/results/netcdf/GP_predict_DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_"+experiments(nn)+areas(rr)+".nc","r")
      GP_mean(tt,rr,nn,:,:)           = f->GP_mean
      GP_SD(tt,rr,nn,:,:)             = f->GP_SD
      GP_mean_sample(tt,rr,nn,:,:,:)  = f->Ysample_masked_reshape(0:19,:,:)
      if (areas(rr) .eq. "_sst") then ; mask land points
        GP_mean(tt,rr,nn,:,:) = where(ls_mask .eq. 0,GP_mean(tt,rr,nn,:,:),GP_mean@_FillValue)
        GP_SD(tt,rr,nn,:,:)   = where(ls_mask .eq. 0,GP_SD(tt,rr,nn,:,:),GP_SD@_FillValue)
      end if
      if (areas(rr) .eq. "_lat") then ; mask ocean points
        GP_mean(tt,rr,nn,:,:) = where(ls_mask .eq. 1,GP_mean(tt,rr,nn,:,:),GP_mean@_FillValue)
        GP_SD(tt,rr,nn,:,:)   = where(ls_mask .eq. 1,GP_SD(tt,rr,nn,:,:),GP_SD@_FillValue)
      end if
      y_obs_weighted_mean = f->y_obs_weighted_mean
     ; proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(tt)+"/DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_"+experiments(nn)+areas(rr)+".txt",5,"float",1)
      GP_mean(tt,rr,nn,:,:)          = GP_mean(tt,rr,nn,:,:) + y_obs_weighted_mean
      GP_mean_sample(tt,rr,nn,:,:,:) = GP_mean_sample(tt,rr,nn,:,:,:) + y_obs_weighted_mean
    end do
  end do
end do

;replace missing values with zeros and add individual SST and LAT predictions for global field
SST_mean            = where(ismissing(GP_mean(:,1,:,:,:)) .eq. True,0,GP_mean(:,1,:,:,:))
LAT_mean            = where(ismissing(GP_mean(:,2,:,:,:)) .eq. True,0,GP_mean(:,2,:,:,:))
SST_SD              = where(ismissing(GP_SD(:,1,:,:,:)) .eq. True,0,GP_SD(:,1,:,:,:))
LAT_SD              = where(ismissing(GP_SD(:,2,:,:,:)) .eq. True,0,GP_SD(:,2,:,:,:))

SST_mean_sample     = where(ismissing(GP_mean_sample(:,1,:,:,:,:)) .eq. True,0,GP_mean_sample(:,1,:,:,:,:))
LAT_mean_sample     = where(ismissing(GP_mean_sample(:,2,:,:,:,:)) .eq. True,0,GP_mean_sample(:,2,:,:,:,:))

GP_mean(:,0,:,:,:)  = SST_mean + LAT_mean
GP_SD(:,0,:,:,:)    = SST_SD   + LAT_SD
GP_mean_sample(:,0,:,:,:,:) = SST_mean_sample + LAT_mean_sample
if (calc_gm .eq. 1) then
  
  ;---part 1: calculate GMST of mean GP fields 
    ;calculate zonal mean SST+LAT individually
    SST_zm    = dim_avg_n_Wrap(GP_mean(:,1,:,:,:),3)
    LAT_zm    = dim_avg_n_Wrap(GP_mean(:,2,:,:,:),3)
  
    ;convert relative anomalies (wrt to present-day reanalysis) to absolute temperatures 
      ;load reanalysis data
      f_HadISST   = addfile(dir+"/reanalysis_reference/HadISST_sst_1979-2018.timmean.zonmean.r72.nc","r")
      f_ERA5      = addfile(dir+"/reanalysis_reference/ERA5_t2m_1979-2018.timmean.zonmean.r72.nc","r")
      HAD_zm      = f_HadISST->sst(0,:,0)
      ERA_zm      = short2flt(f_ERA5->t2m(0,:,0))
      ERA_zm      = ERA_zm-273.15
      ;replace high-latitude missing values in modern geography
      HAD_zm(4)   = HAD_zm(5)
      HAD_zm(71)  = HAD_zm(70)
      ;add modern reference temperatures to get absolute temperatures
      SST_zm      = SST_zm + conform_dims(dimsizes(SST_zm),HAD_zm,2)
      LAT_zm      = LAT_zm + conform_dims(dimsizes(LAT_zm),ERA_zm,2)
      ;replace missing values with zeroes
      SST_zm_fill = where(ismissing(SST_zm) .eq. True,0,SST_zm)
      LAT_zm_fill = where(ismissing(LAT_zm) .eq. True,0,LAT_zm)

    ;weight zonal mean surface temperature by respective land fraction at each latitude
    land_frac = dim_avg_n_Wrap(ls_mask,1)
    T_zm      = LAT_zm_fill*conform_dims(dimsizes(LAT_zm_fill),land_frac,2) + SST_zm_fill*(1.-conform_dims(dimsizes(SST_zm_fill),land_frac,2))

    rad       = 4.0*atan(1.0)/180.0
    clat      = cos(lat*rad)
    clat!0    = "lat"
    clat&lat  = lat

    T_gm      = dim_avg_wgt_n(T_zm,clat,1,2)

  ;---part 2: calculate GMST for N=10000 subsamples from posterior
    N                 = 20
    GP_sample_mean    = new((/dimsizes(time_slices),dimsizes(areas),dimsizes(experiments),N,dimsizes(lat)/),"float")

    do tt=0,dimsizes(time_slices)-1
      do rr=0,dimsizes(areas)-1
        do nn=0,dimsizes(experiments)-1
          f                             = addfile(dir+"/results/netcdf/GP_predict_DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_"+experiments(nn)+areas(rr)+".nc","r")
          GP_sample_mean(tt,rr,nn,:,:)  = f->Ysample_masked_reshape_zm
          y_obs_weighted_mean = f->y_obs_weighted_mean
        ;  proxy                        := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(tt)+"/DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_"+experiments(nn)+areas(rr)+".txt",5,"float",1)
          GP_sample_mean(tt,rr,nn,:,:)  = GP_sample_mean(tt,rr,nn,:,:) + y_obs_weighted_mean
        end do
      end do
    end do

    ;add modern reference temperatures to get absolute temperatures
    GP_sample_mean(:,1,:,:,:)      = GP_sample_mean(:,1,:,:,:) + conform_dims(dimsizes(GP_sample_mean(:,1,:,:,:)),HAD_zm,3)
    GP_sample_mean(:,2,:,:,:)      = GP_sample_mean(:,2,:,:,:) + conform_dims(dimsizes(GP_sample_mean(:,2,:,:,:)),ERA_zm,3)

    ;replace missing values with zeroes
    GP_sample_mean(:,1,:,:,:) = where(ismissing(GP_sample_mean(:,1,:,:,:)) .eq. True,0,GP_sample_mean(:,1,:,:,:))
    GP_sample_mean(:,2,:,:,:) = where(ismissing(GP_sample_mean(:,2,:,:,:)) .eq. True,0,GP_sample_mean(:,2,:,:,:))

    T_sample_zm      = GP_sample_mean(:,2,:,:,:)*conform_dims(dimsizes(GP_sample_mean(:,2,:,:,:)),land_frac,3) + GP_sample_mean(:,1,:,:,:)*(1.-conform_dims(dimsizes(GP_sample_mean(:,1,:,:,:)),land_frac,3))
    T_sample_gm      = dim_avg_wgt_n(T_sample_zm,clat,1,3)
    T_sample_gm_sd   = dim_stddev_n(T_sample_gm,2)
    
    header  = (/"----------------------------------------------------------------", \
               "GMST estimates (+/- 1 SD) from the Gaussian Process", \ 
               "regression using the "+ref_frame+" reference frame ", \
               "----------------------------------------------------------------"/)

    hlist   = [/header/]
    cnames  = [/"Experiment","   |       LP","      |      PETM ","    |     EECO   "/]
    
    rule    = [/(/"----------------------------------------------------------------"/)/]
  
    labels_table      = (/"Default     ","Frosty      ","NoTEX       ","NoMBT       ","NoMammal    "/)
    
    fname   = dir+"/results/GP_results_"+ref_frame+"_"+cut_off+"_degC_SD_cut_off.txt"
    write_table(fname, "w", hlist, "%s")
    write_table(fname, "a", cnames, "%s %s %s %s ")
    write_table(fname, "a", rule, "%s")
    
    SD  = 1.0
    
    do nn=0,4
      results    = [/(/labels_table(nn)/),T_gm(0,nn),T_sample_gm_sd(0,nn),T_gm(1,nn),T_sample_gm_sd(1,nn),T_gm(2,nn),T_sample_gm_sd(2,nn)/]
      write_table(fname, "a", results, "%s  | %5.1f (+/-%2.1f) | %5.1f (+/-%2.1f) | %5.1f (+/-%2.1f)")
    end do
    
end if

if (plot_time_slice_maps .eq. 1) then
  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_time_slice_maps_"+ref_frame+"_new")

  res                    := True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@tiXAxisFontHeightF  = 0.02
  res@tiYAxisFontHeightF  = 0.02
  res@mpDataBaseVersion   = "Ncarg4_1"
  res@mpDataSetName       = "herold_eocene"

  res@tmXBMode            = "Explicit"	
  res@tmXBValues          = (/-180,-120,-60,0,60,120,180/)
  res@tmXBLabels          = (/"180~S~o~N~W","120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~E"/)
  res@tmYLMode            = "Explicit"	
  res@tmYLValues          = (/-90,-60,-30,0,30,60,90/)
  res@tmYLLabels          = (/"90~S~o~N~S","60~S~o~N~S","30~S~o~N~S","0~S~o~N~","30~S~o~N~N","60~S~o~N~N","90~S~o~N~N"/)
  res@cnFillOn            = True
  res@cnLinesOn           = False
  res@mpFillOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  cn_levels_mean         := ispan(-28,28,4)
  res@cnLevels            = cn_levels_mean
  res@cnFillPalette       = "amwg_blueyellowred"
  res@mpGeophysicalLineThicknessF = 1.5       
  res@lbLabelBarOn	              = False
  res@lbOrientation               = "horizontal"
  res@lbLabelFontHeightF	        =	0.02
  res@lbTitleString               = "~F8~D~F21~GMST (~S~o~N~C)"
  res@lbTitleFontHeightF          = 0.02
;  res@lbLabelStride               = 2
  res@pmLabelBarWidthF            = 0.39
  res@pmLabelBarHeightF           = 0.1
  res@lbTitlePosition             = "Bottom"
  res@mpOutlineDrawOrder          = "Draw"
  res@cnFillDrawOrder             = "PreDraw"
  
  res@gsnStringFontHeightF  = 0.02
;  res@gsnLeftStringOrthogonalPosF = 0.02
  res@tmXTOn                    = False
  res@mpCenterLonF          = -90.
  res@mpGridAndLimbOn       = True
  res@mpGridLineColor       = "gray50"
  res@mpGridLineDashPattern = 16
  res@mpLimbLineThicknessF  = 1.5
  res@mpGridLineThicknessF  = .1

  res@mpGridLatSpacingF     = 30. 
  res@mpGridLonSpacingF     = 60. 
  res@mpProjection = "Robinson"
  res@mpPerimOn     = False

  left_strings  = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)
  center_strings= (/"D~B~surf~N~-3-default LP","D~B~surf~N~-3-default LP","D~B~surf~N~-3-default PETM","D~B~surf~N~-3-default PETM","D~B~surf~N~-3-default EECO","D~B~surf~N~-3-default EECO"/)


  res_sd                  = res
;  cn_levels_sd           := (/3.5,4.,4.5,5.,5.5,6,7,8,10,12/)
  cn_levels_sd           := ispan(4,13,1)

  res_sd@cnLevels        := cn_levels_sd
  cmap_sd                 = "WhiteBlueGreenYellowRed"
  res_sd@cnFillPalette    = cmap_sd
  res_sd@lbTitleString               = "Standard Deviation (~S~o~N~C)"

  resP              = True
  resP@gsnFrame     = False
;  resP@gsnMaximize  = True
  resP@gsnPanelFigureStringsFontHeightF = 0.016
  resP@amOrthogonalPosF = -.62
  resP@amJust         = "TopLeft"
  resP@gsnPanelYWhiteSpacePercent = 2.
  resP@gsnPanelBottom = .1
  resP@gsnPanelFigureStrings = left_strings

  pmres                  := True
  pmres@gsMarkerColor    = "black"
  pmres@gsMarkerThicknessF  = 1.
  pmres@gsMarkerSizeF     = 0.01

;  pmres@tfPolyDrawOrder   = "PostDraw"

  mstring = "u"
  fontnum = 34
  xoffset = 0.0
  yoffset = 0.0
  ratio   = 1.2
  size    = 1.0
  angle   = 0.0
  new_index = NhlNewMarker(wks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)
  
  pfres                   = pmres


  txres                   = True
  txres@txFontHeightF     = 0.02

  plot           :=  new(6,graphic)
  poly_outline1  :=  new(200*3,graphic)
  poly_fill1     :=  new(200*3,graphic)
  poly_outline2  :=  new(200*3,graphic)
  poly_fill2     :=  new(200*3,graphic)

  do tt=0,2
    if (tt .eq. 2) then
      res@lbLabelBarOn	              = True
      res_sd@lbLabelBarOn	              = True
    end if

;    res@gsnLeftString    = left_strings(tt*2)
;    res_sd@gsnLeftString = left_strings(tt*2+1)
    res@gsnCenterString    = center_strings(tt*2)
    res_sd@gsnCenterString = center_strings(tt*2+1)

    plot(tt*2)    = gsn_csm_contour_map(wks,GP_mean(tt,0,0,:,:),res)
    plot(tt*2+1)  = gsn_csm_contour_map(wks,GP_SD(tt,0,0,:,:),res_sd)

    ;load and plot proxy data
    proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(tt)+"/DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_ctrl_all.txt",5,"float",1)
    dims_proxy         := dimsizes(proxy)
    colors_mean        := get_color_rgba("amwg_blueyellowred",cn_levels_mean,proxy(:,2))
    colors_sd          := get_color_rgba(cmap_sd,cn_levels_sd,proxy(:,3))
    lon_orig           := proxy(:,0)
    lat_orig           := proxy(:,1)
    lon_new            := lon_orig
    lat_new            := lat_orig
    unique_lons        := get_unique_values(lon_orig)
    unique_lons_flag   := toint(unique_lons)
    unique_lons_flag    = 0
    
    do nn=0,dims_proxy(0)-1
      index = ind(proxy(nn,0) .eq. unique_lons)
      if (unique_lons_flag(index) .gt. 1) then
        if (unique_lons_flag(index)/2 .eq. unique_lons_flag(index)/2.) then
          lat_new(nn) = lat_new(nn)+unique_lons_flag(index)*1
        else
          lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*1
        end if
      else if (unique_lons_flag(index) .eq. 1) then
        lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*2
      end if
      end if

      unique_lons_flag(index) = unique_lons_flag(index) + 1
      if (proxy(nn,4) .eq. 1) then ; SST location
          pfres@gsMarkerIndex    = 16
          pmres@gsMarkerIndex    = 4
          pfres@gsMarkerSizeF    = 0.01
      else
          pfres@gsMarkerIndex    = new_index
          pmres@gsMarkerIndex    = 7
          pfres@gsMarkerSizeF    = 0.011
      end if      

      pfres@gsMarkerColor     := colors_mean(nn,:)
      poly_fill1(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt*2),proxy(nn,0),lat_new(nn),pfres)
      poly_outline1(tt*200+nn) = gsn_add_polymarker(wks,plot(tt*2),proxy(nn,0),lat_new(nn),pmres)
      pfres@gsMarkerColor     := colors_sd(nn,:)
      poly_fill2(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt*2+1),proxy(nn,0),lat_new(nn),pfres)
      poly_outline2(tt*200+nn) = gsn_add_polymarker(wks,plot(tt*2+1),proxy(nn,0),lat_new(nn),pmres) 
    end do

  end do

;  gsn_text_ndc(wks,(/"temperature anomaly","standard deviation"/),(/0.26,0.765/),(/0.935,0.935/),txres)
  gsn_panel(wks,plot,(/3,2/),resP)
  
  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_time_slice_maps_means_"+ref_frame)

  left_strings  := (/"(a)","(b)","(c)"/)
  center_strings:= (/"D~B~surf~N~-3-default LP","D~B~surf~N~-3-default PETM","D~B~surf~N~-3-default EECO"/)

  resP@gsnPanelYWhiteSpacePercent = 2.
  resP@gsnPanelFigureStrings := left_strings
  resP@gsnPanelBottom = .12

  mstring = "u"
  fontnum = 34
  xoffset = 0.0
  yoffset = 0.0
  ratio   = 1.2
  size    = 1.0
  angle   = 0.0
  new_index = NhlNewMarker(wks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)
  
  plot           :=  new(3,graphic)
  poly_outline1  :=  new(200*3,graphic)
  poly_fill1     :=  new(200*3,graphic)
  poly_outline2  :=  new(200*3,graphic)
  poly_fill2     :=  new(200*3,graphic)

  do tt=0,2
    if (tt .eq. 2) then
      res@lbLabelBarOn	              = True
    else
      res@lbLabelBarOn	              = False
    end if

    res@gsnCenterString    = center_strings(tt)
    plot(tt)    = gsn_csm_contour_map(wks,GP_mean(tt,0,0,:,:),res)

    ;load and plot proxy data
    proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(tt)+"/DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_ctrl_all.txt",5,"float",1)
    dims_proxy         := dimsizes(proxy)
    colors_mean        := get_color_rgba("amwg_blueyellowred",cn_levels_mean,proxy(:,2))
    colors_sd          := get_color_rgba(cmap_sd,cn_levels_sd,proxy(:,3))
    lon_orig           := proxy(:,0)
    lat_orig           := proxy(:,1)
    lon_new            := lon_orig
    lat_new            := lat_orig
    unique_lons        := get_unique_values(lon_orig)
    unique_lons_flag   := toint(unique_lons)
    unique_lons_flag    = 0
    
    do nn=0,dims_proxy(0)-1
      index = ind(proxy(nn,0) .eq. unique_lons)
      if (unique_lons_flag(index) .gt. 1) then
        if (unique_lons_flag(index)/2 .eq. unique_lons_flag(index)/2.) then
          lat_new(nn) = lat_new(nn)+unique_lons_flag(index)*1
        else
          lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*1
        end if
      else if (unique_lons_flag(index) .eq. 1) then
        lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*2
      end if
      end if

      unique_lons_flag(index) = unique_lons_flag(index) + 1
      if (proxy(nn,4) .eq. 1) then ; SST location
          pfres@gsMarkerIndex    = 16
          pmres@gsMarkerIndex    = 4
          pfres@gsMarkerSizeF    = 0.01
      else
          pfres@gsMarkerIndex    = new_index
          pmres@gsMarkerIndex    = 7
          pfres@gsMarkerSizeF    = 0.011
      end if      

      pfres@gsMarkerColor     := colors_mean(nn,:)
      poly_fill1(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt),proxy(nn,0),lat_new(nn),pfres)
      poly_outline1(tt*200+nn) = gsn_add_polymarker(wks,plot(tt),proxy(nn,0),lat_new(nn),pmres)
    end do

  end do
  gsn_panel(wks,plot,(/3,1/),resP)
  

  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_eeco_uncertainty_"+ref_frame)

  left_strings  := (/"(a)","(b)","(c)"/)
;  center_strings:= (/"D~B~surf~N~-3-default LP","D~B~surf~N~-3-default PETM","D~B~surf~N~-3-default EECO"/)

;=  resP@gsnPanelYWhiteSpacePercent = 2.
;  resP@gsnPanelFigureStrings := left_strings
;  resP@gsnPanelBottom = .12

  mstring = "u"
  fontnum = 34
  xoffset = 0.0
  yoffset = 0.0
  ratio   = 1.2
  size    = 1.0
  angle   = 0.0
  new_index = NhlNewMarker(wks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)
  
  plot           :=  new(20,graphic)
  poly_outline1  :=  new(200*20,graphic)
  poly_fill1     :=  new(200*20,graphic)
  poly_outline2  :=  new(200*20,graphic)
  poly_fill2     :=  new(200*20,graphic)

  do tt=2,2
    res@lbLabelBarOn	              = False
    ;load and plot proxy data
    proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(tt)+"/DeepMIP_"+ref_frame+"_"+time_slices(tt)+"_ctrl_all.txt",5,"float",1)
    dims_proxy         := dimsizes(proxy)
    colors_mean        := get_color_rgba("amwg_blueyellowred",cn_levels_mean,proxy(:,2))
    colors_sd          := get_color_rgba(cmap_sd,cn_levels_sd,proxy(:,3))
    lon_orig           := proxy(:,0)
    lat_orig           := proxy(:,1)
    lon_new            := lon_orig
    lat_new            := lat_orig
    unique_lons        := get_unique_values(lon_orig)
    unique_lons_flag   := toint(unique_lons)
    unique_lons_flag    = 0

    do ss=0,19
      if (ss .eq. 18) then
        res@lbLabelBarOn	              = True
      end if
      
      if (ss .eq. 0) then
        res@gsnCenterString    = "GP mean field EECO"
        plot(ss)    = gsn_csm_contour_map(wks,GP_mean(tt,0,0,:,:),res)
      else
        res@gsnCenterString    = "posterior supsampling #"+tostring(ss)
        plot(ss)    = gsn_csm_contour_map(wks,GP_mean_sample(tt,0,0,ss,:,:),res)
      end if
    
      do nn=0,dims_proxy(0)-1
        index = ind(proxy(nn,0) .eq. unique_lons)
        if (unique_lons_flag(index) .gt. 1) then
          if (unique_lons_flag(index)/2 .eq. unique_lons_flag(index)/2.) then
            lat_new(nn) = lat_new(nn)+unique_lons_flag(index)*1
          else
            lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*1
          end if
        else if (unique_lons_flag(index) .eq. 1) then
          lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*2
        end if
        end if

        unique_lons_flag(index) = unique_lons_flag(index) + 1
        if (proxy(nn,4) .eq. 1) then ; SST location
            pfres@gsMarkerIndex    = 16
            pmres@gsMarkerIndex    = 4
            pfres@gsMarkerSizeF    = 0.01
        else
            pfres@gsMarkerIndex    = new_index
            pmres@gsMarkerIndex    = 7
            pfres@gsMarkerSizeF    = 0.011
        end if      

        pfres@gsMarkerColor     := colors_mean(nn,:)
        poly_fill1(ss*200+nn)    = gsn_add_polymarker(wks,plot(ss),proxy(nn,0),lat_new(nn),pfres)
        poly_outline1(ss*200+nn) = gsn_add_polymarker(wks,plot(ss),proxy(nn,0),lat_new(nn),pmres)
      end do
    end do

  end do
  
  delete(resP@gsnPanelBottom)
  delete(resP@gsnPanelFigureStrings)
  gsn_panel(wks,plot,(/6,3/),resP)

end if

if (plot_sst_lat_maps .eq. 1) then
  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_sst_lat_maps_"+ref_frame)

  res                    := True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@tiXAxisFontHeightF  = 0.02
  res@tiYAxisFontHeightF  = 0.02
  res@mpDataBaseVersion   = "Ncarg4_1"
  res@mpDataSetName       = "herold_eocene"

  res@tmXBMode            = "Explicit"	
  res@tmXBValues          = (/-180,-120,-60,0,60,120,180/)
  res@tmXBLabels          = (/"180~S~o~N~W","120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~E"/)
  res@tmYLMode            = "Explicit"	
  res@tmYLValues          = (/-90,-60,-30,0,30,60,90/)
  res@tmYLLabels          = (/"90~S~o~N~S","60~S~o~N~S","30~S~o~N~S","0~S~o~N~","30~S~o~N~N","60~S~o~N~N","90~S~o~N~N"/)
  res@cnFillOn            = True
  res@cnLinesOn           = False
  res@mpFillOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  cn_levels_mean         := ispan(-28,28,4)
  res@cnLevels            = cn_levels_mean
  res@cnFillPalette       = "amwg_blueyellowred"
  res@mpGeophysicalLineThicknessF = 1.5       
  res@lbLabelBarOn	              = False
  res@lbOrientation               = "horizontal"
  res@lbLabelFontHeightF	        =	0.025
  res@lbTitleString               = "~F8~D~F21~GMST (~S~o~N~C)"
  res@lbTitleFontHeightF          = 0.025
  res@pmLabelBarWidthF            = 0.6
  res@pmLabelBarOrthogonalPosF           = 0.12
  res@lbTitlePosition             = "Bottom"
  res@mpOutlineDrawOrder          = "Draw"
  res@cnFillDrawOrder             = "PreDraw"
  res@cnMissingValFillColor       = "gray95"
  res_sd                  = res

  res_sd@cnLevels        := cn_levels_sd
  cmap_sd                 = "WhiteBlueGreenYellowRed"
  res_sd@cnFillPalette    = cmap_sd

  resP              = True
  resP@gsnFrame     = False
;  resP@gsnMaximize  = True
  resP@gsnPanelFigureStringsFontHeightF = 0.013
  resP@amJust         = "TopLeft"
  resP@gsnPanelFigureStrings  := (/"(a) SST","(b) SST","(c) LAT","(d) LAT"/)

  mstring = "u"
  fontnum = 34
  xoffset = 0.0
  yoffset = 0.0
  ratio   = 1.2
  size    = 1.0
  angle   = 0.0
  new_index = NhlNewMarker(wks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)

  txres                   = True
  txres@txFontHeightF     = 0.02

  plot           :=  new(4,graphic)
  poly_outline1  :=  new(200*2,graphic)
  poly_fill1     :=  new(200*2,graphic)
  poly_outline2  :=  new(200*2,graphic)
  poly_fill2     :=  new(200*2,graphic)

  do tt=0,1
    if (tt .eq. 1) then
      res@lbLabelBarOn	              = True
      res_sd@lbLabelBarOn	              = True
    end if

    plot(tt*2)    = gsn_csm_contour_map_ce(wks,GP_mean(2,1+tt,0,:,:),res)
    plot(tt*2+1)  = gsn_csm_contour_map_ce(wks,GP_SD(2,1+tt,0,:,:),res_sd)

    ;load and plot proxy data
    proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_eeco/DeepMIP_"+ref_frame+"_eeco_ctrl"+areas(1+tt)+"_all.txt",5,"float",1)
    dims_proxy         := dimsizes(proxy)
    colors_mean        := get_color_rgba("amwg_blueyellowred",cn_levels_mean,proxy(:,2))
    colors_sd          := get_color_rgba(cmap_sd,cn_levels_sd,proxy(:,3))
    lon_orig           := proxy(:,0)
    lat_orig           := proxy(:,1)
    lon_new            := lon_orig
    lat_new            := lat_orig
    unique_lons        := get_unique_values(lon_orig)
    unique_lons_flag   := toint(unique_lons)
    unique_lons_flag    = 0
    

    do nn=0,dims_proxy(0)-1
      index = ind(proxy(nn,0) .eq. unique_lons)
      if (unique_lons_flag(index) .gt. 1) then
        if (unique_lons_flag(index)/2 .eq. unique_lons_flag(index)/2.) then
          lat_new(nn) = lat_new(nn)+unique_lons_flag(index)*1
        else
          lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*1
        end if
      else if (unique_lons_flag(index) .eq. 1) then
        lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*2
      end if
      end if

      if (proxy(nn,4) .eq. 1) then ; SST location
          pfres@gsMarkerIndex    = 16
          pmres@gsMarkerIndex    = 4
          pfres@gsMarkerSizeF    = 0.01
      else
          pfres@gsMarkerIndex    = new_index
          pmres@gsMarkerIndex    = 7
          pfres@gsMarkerSizeF    = 0.011
      end if      

      unique_lons_flag(index) = unique_lons_flag(index) + 1

      pfres@gsMarkerColor     := colors_mean(nn,:)
      poly_fill1(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt*2),proxy(nn,0),lat_new(nn),pfres)
      poly_outline1(tt*200+nn) = gsn_add_polymarker(wks,plot(tt*2),proxy(nn,0),lat_new(nn),pmres)
      pfres@gsMarkerColor     := colors_sd(nn,:)
      poly_fill2(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt*2+1),proxy(nn,0),lat_new(nn),pfres)
      poly_outline2(tt*200+nn) = gsn_add_polymarker(wks,plot(tt*2+1),proxy(nn,0),lat_new(nn),pmres) 
    end do

  end do

  gsn_text_ndc(wks,(/"temperature anomaly","standard deviation"/),(/0.26,0.765/),(/0.775,0.775/),txres)
  gsn_panel(wks,plot,(/2,2/),resP)

end if

if (plot_subsampling_maps .eq. 1) then

  do pp=0,2

    wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_subsampling_maps_"+time_slices(pp)+"_"+ref_frame)

  ;  resP@gsnPanelFigureStrings  := (/"(a) "+labels(0),"(b) "+labels(1),"(c) "+labels(2),"(d) "+labels(3),"(e) "+labels(4),"(f) "+labels(5)/)
  ;  resP@gsnPanelFigureStrings  := (/"(a)","(b)","(c)","(d)"/)
  ;  resP@gsnPanelFigureStrings  := (/"(a) "+labels(1),"(b) "+labels(2),"(c) "+labels(3),"(d) "+labels(4)/)

    mstring = "u"
    fontnum = 34
    xoffset = 0.0
    yoffset = 0.0
    ratio   = 1.2
    size    = 1.0
    angle   = 0.0
    new_index = NhlNewMarker(wks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)

    left_strings  := (/"(a)","(b)","(c)","(d)"/)
    center_strings  := (/"D~B~surf~N~-3-Frosty","D~B~surf~N~-3-NoTEX","D~B~surf~N~-3-NoMBT","D~B~surf~N~-3-NoMammal"/)

    res@pmLabelBarParallelPosF  = -0.00
    resP@gsnPanelFigureStrings  := left_strings
    
    txres                   = True
    txres@txFontHeightF     = 0.02

    plot           :=  new(4,graphic)
    poly_outline1  :=  new(200*6,graphic)
    poly_fill1     :=  new(200*6,graphic)

    do tt=0,3
      if (tt .eq. 3) then
        res@lbLabelBarOn	              = True
      else
        res@lbLabelBarOn	              = False
      end if
    
      res@gsnCenterString = center_strings(tt)
      plot(tt)    = gsn_csm_contour_map(wks,GP_mean(pp,0,tt+1,:,:),res)

      ;load and plot proxy data
      proxy              := readAsciiTable(dir+"/Observation_Data/DeepMIP_"+time_slices(pp)+"/DeepMIP_"+ref_frame+"_"+time_slices(pp)+"_"+experiments(tt+1)+"_all.txt",5,"float",1)

      dims_proxy         := dimsizes(proxy)
      colors_mean        := get_color_rgba("amwg_blueyellowred",cn_levels_mean,proxy(:,2))
      lon_orig           := proxy(:,0)
      lat_orig           := proxy(:,1)
      lon_new            := lon_orig
      lat_new            := lat_orig
      unique_lons        := get_unique_values(lon_orig)
      unique_lons_flag   := toint(unique_lons)
      unique_lons_flag    = 0
    

      do nn=0,dims_proxy(0)-1
        index = ind(proxy(nn,0) .eq. unique_lons)
        if (unique_lons_flag(index) .gt. 1) then
          if (unique_lons_flag(index)/2 .eq. unique_lons_flag(index)/2.) then
            lat_new(nn) = lat_new(nn)+unique_lons_flag(index)*1
          else
            lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*1
          end if
        else if (unique_lons_flag(index) .eq. 1) then
          lat_new(nn) = lat_new(nn)-unique_lons_flag(index)*2
        end if
        end if

        unique_lons_flag(index) = unique_lons_flag(index) + 1
      if (proxy(nn,4) .eq. 1) then ; SST location
          pfres@gsMarkerIndex    = 16
          pmres@gsMarkerIndex    = 4
          pfres@gsMarkerSizeF    = 0.009
      else
          pfres@gsMarkerIndex    = new_index
          pmres@gsMarkerIndex    = 7
          pfres@gsMarkerSizeF    = 0.01
      end if      

        pfres@gsMarkerColor     := colors_mean(nn,:)
        poly_fill1(tt*200+nn)    = gsn_add_polymarker(wks,plot(tt),proxy(nn,0),lat_new(nn),pfres)
        poly_outline1(tt*200+nn) = gsn_add_polymarker(wks,plot(tt),proxy(nn,0),lat_new(nn),pmres)
      end do

    end do

    gsn_panel(wks,plot,(/2,2/),resP)
    
  end do
  
end if

if (plot_gm .eq. 1) then
  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_global_means_"+ref_frame)

  res                    := True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@tiXAxisFontHeightF  = 0.02
  res@tiYAxisFontHeightF  = 0.01
  res@trXMinF             = 0.5
  res@trXMaxF             = tofloat(dimsizes(experiments))+.5
  res@xyMarkLineModes     = "Markers"
  res@xyMarkerSizeF       = 0.008
  res@xyMarkers           = 16.
  res@vpHeightF           = 0.2
  res@vpWidthF            = 0.25
  res@tmXBMode            = "Explicit"	
  res@tmXBValues          = ispan(1,dimsizes(experiments),1)
  res@tmXBLabels          = labels
  res@tmXBLabelAngleF     = -45
  res@tmXBLabelJust       = "CenterCenter"
  res@tmYROn              = True
  res@tmYRLabelsOn        = True
  res@tiYAxisSide         = "Left"
  res@trYMinF             = -2
  res@trYMaxF             = 26
  res@xyCurveDrawOrder    = "PostDraw"
  res@tmYLMode            = "Explicit"	
  res@tmYLValues          = (/0,4,8,12,16,20,24/)
  res@tmYLLabels          = tostring(res@tmYLValues)
  res@tmYLMinorValues     = ispan(-10,30,1)

  
  res_hollow                    = res
  res_hollow@xyMarkers          = 4.
  res_hollow@xyMarkerColor      = "black"
  res_hollow@xyMarkerThicknessF = 1.5

  pmres                  := True
  pmres@gsMarkerIndex    = 4
  pmres@gsMarkerSizeF    = 0.008
  pmres@gsMarkerColor    = "black"
  pmres@gsMarkerThicknessF  = 1.5
  pmres@tfPolyDrawOrder   = "PostDraw"

  pgres                  := True
  pgres@gsLineThicknessF  = 1.5
  pgres@tfPolyDrawOrder   = "PostDraw"
  pgres@gsLineColor       = "gray50"

  polyres                := True
  polyres@gsEdgesOn      = False
  polyres@gsFillColor   := "gray90"
  polyres@tfPolyDrawOrder = "PreDraw"

  resP              = True
  resP@gsnFrame     = False
;  resP@gsnMaximize  = True
  resP@gsnPanelFigureStringsFontHeightF = 0.025
  resP@amJust         = "TopLeft"
;  resP@gsnPanelFigureStrings  := (/"(a) SST+LAT","(b) SST","(c) LAT"/)
  resP@gsnPanelFigureStrings  := ""

  plot1 :=  new(3,graphic)
  plot2 :=  new(3,graphic)
  plot3 :=  new(3,graphic)
  mark1 :=  new(3,graphic)
  mark2 :=  new(3,graphic)
  mark3 :=  new(3,graphic)

  poly1  = new(dimsizes(experiments)*3,graphic)
  poly2  = new(dimsizes(experiments)*3,graphic)
  poly3  = new(dimsizes(experiments)*3,graphic)
  poly1a = new(dimsizes(experiments)*3,graphic)
  poly2a = new(dimsizes(experiments)*3,graphic)
  poly3a = new(dimsizes(experiments)*3,graphic)
  poly1b = new(dimsizes(experiments)*3,graphic)
  poly2b = new(dimsizes(experiments)*3,graphic)
  poly3b = new(dimsizes(experiments)*3,graphic)

  bpoly1  = new(dimsizes(experiments)*3,graphic)
  bpoly2  = new(dimsizes(experiments)*3,graphic)
  bpoly3  = new(dimsizes(experiments)*3,graphic)
  bpoly1a = new(dimsizes(experiments)*3,graphic)
  bpoly2a = new(dimsizes(experiments)*3,graphic)
  bpoly3a = new(dimsizes(experiments)*3,graphic)
  bpoly1b = new(dimsizes(experiments)*3,graphic)
  bpoly2b = new(dimsizes(experiments)*3,graphic)
  bpoly3b = new(dimsizes(experiments)*3,graphic)

  circ1 = new(dimsizes(experiments)*3,graphic)
  circ2 = new(dimsizes(experiments)*3,graphic)
  circ3 = new(dimsizes(experiments)*3,graphic)

  shade1 = new(dimsizes(experiments)*3,graphic)
  shade2 = new(dimsizes(experiments)*3,graphic)

  col1  :=  (/27,158,119/)/256.
  col2  :=  (/217,95,2/)/256.
  col3  :=  (/117,112,179/)/256.

  width = 0.2
  ystrings  = (/"GMST increase [~S~o~N~C]","GMSST increase [~S~o~N~C]"/)

  do rr=0,0
    res@tiYAxisString       = ystrings(rr)

    yspan  = fspan(1-width,tofloat(dimsizes(experiments))-width,dimsizes(experiments))
    res@xyMarkerColor    = col1
    plot1(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(0,rr,:),res)
    mark1(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(0,rr,:),res_hollow)

    yspan_shade         = fspan(1., tofloat(dimsizes(experiments)),dimsizes(experiments))
    do pp=0,dimsizes(experiments)-1
      if (pp/2 .ne. pp/2.) then
        shade1(rr*dimsizes(experiments)+pp)  = gsn_add_polygon(wks,plot1(rr),(/yspan_shade(pp)-2.5*width,yspan_shade(pp)+2.5*width,yspan_shade(pp)+2.5*width,yspan_shade(pp)-2.5*width,yspan_shade(pp)-2.5*width/),(/-30,-30,30,30,-30/),polyres)
      end if
    end do 

    do pp=0,dimsizes(experiments)-1
      poly1(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot1(rr),(/yspan(pp),yspan(pp)/),(/GP_mean_gm(0,rr,pp)-GP_SD_gm(0,rr,pp),GP_mean_gm(0,rr,pp)+GP_SD_gm(0,rr,pp)/),pgres)
      poly1a(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot1(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(0,rr,pp)-GP_SD_gm(0,rr,pp),GP_mean_gm(0,rr,pp)-GP_SD_gm(0,rr,pp)/),pgres)
      poly1b(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot1(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(0,rr,pp)+GP_SD_gm(0,rr,pp),GP_mean_gm(0,rr,pp)+GP_SD_gm(0,rr,pp)/),pgres)
    end do

    yspan               = fspan(1., tofloat(dimsizes(experiments)),dimsizes(experiments))
    res@xyMarkerColor   = col2
    plot2(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(1,rr,:),res)
    mark2(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(1,rr,:),res_hollow)


    do pp=0,dimsizes(experiments)-1
      poly2(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot2(rr),(/yspan(pp),yspan(pp)/),(/GP_mean_gm(1,rr,pp)-GP_SD_gm(1,rr,pp),GP_mean_gm(1,rr,pp)+GP_SD_gm(1,rr,pp)/),pgres)
      poly2a(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot2(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(1,rr,pp)-GP_SD_gm(1,rr,pp),GP_mean_gm(1,rr,pp)-GP_SD_gm(1,rr,pp)/),pgres)
      poly2b(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot2(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(1,rr,pp)+GP_SD_gm(1,rr,pp),GP_mean_gm(1,rr,pp)+GP_SD_gm(1,rr,pp)/),pgres)
    end do

    yspan               = fspan(1+width, tofloat(dimsizes(experiments))+width,dimsizes(experiments))
    res@xyMarkerColor   = col3
    plot3(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(2,rr,:),res)
    mark3(rr)            = gsn_csm_xy(wks,yspan,GP_mean_gm(2,rr,:),res_hollow)

    do pp=0,dimsizes(experiments)-1
      poly3(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot3(rr),(/yspan(pp),yspan(pp)/),(/GP_mean_gm(2,rr,pp)-GP_SD_gm(2,rr,pp),GP_mean_gm(2,rr,pp)+GP_SD_gm(2,rr,pp)/),pgres)
      poly3a(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot3(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(2,rr,pp)-GP_SD_gm(2,rr,pp),GP_mean_gm(2,rr,pp)-GP_SD_gm(2,rr,pp)/),pgres)
      poly3b(rr*dimsizes(experiments)+pp)  = gsn_add_polyline(wks,plot3(rr),(/yspan(pp)-(width/3),yspan(pp)+(width/3)/),(/GP_mean_gm(2,rr,pp)+GP_SD_gm(2,rr,pp),GP_mean_gm(2,rr,pp)+GP_SD_gm(2,rr,pp)/),pgres)
    end do

    overlay(plot1(rr),mark1(rr))
    overlay(plot2(rr),mark2(rr))
    overlay(plot3(rr),mark3(rr))
    overlay(plot1(rr),plot2(rr))
    overlay(plot1(rr),plot3(rr))
  end do

  res1                     = True
  res2                     = True
;  res2@gsnAttachPlotsXAxis = True
;  anno1 = gsn_attach_plots(plot1(0),plot1(1),res1,res2)
;  anno2 = gsn_attach_plots(plot1(0),plot1(1),res1,res2)

  gsn_panel(wks,plot1,(/1,1/),resP)

;  draw(plot1)

end if

if (plot_zm .eq. 1) then
  wks = gsn_open_wks("pdf",dir+"/results/plots/GP_regression_zonal_means_"+ref_frame)

  printVarSummary(T_sample_zm)
  
  res_zm               = True
  res_zm@gsnDraw       = False
  res_zm@gsnFrame      = False
  res_zm@tiXAxisString       = "zonal mean Tsurf"
  res_zm@tiXAxisFontHeightF  = 0.02
  res_zm@tiYAxisFontHeightF  = 0.02
  res_zm@trYMinF             = -90
  res_zm@trYMaxF             = 90
  res_zm@trXMinF             = -60
  res_zm@trXMaxF             = 60
  res_zm@gsnXRefLine         = 0.              ; reference line   
  res_zm@gsnXRefLineDashPattern  = 2
  res_zm@vpWidthF  = 0.3
  res_zm@vpHeightF = 0.8
  res_zm@tmYLMode        = "Explicit"	
  res_zm@tmYLValues      = (/-90,-60,-30,0,30,60,90/)
  res_zm@tmYLLabels      = (/"90~S~o~N~S","60~S~o~N~S","30~S~o~N~S","0~S~o~N~","30~S~o~N~N","60~S~o~N~N","90~S~o~N~N"/) 
  res_zm@tmXBMode        = "Explicit"	
  res_zm@tmXBValues      = (/-40,-20,0,20,40,60/)
  res_zm@tmXBLabels      = "" + res_zm@tmXBValues
  res_zm@pmLegendDisplayMode  = "NoCreate"
  
  res_zm_mean     = res_zm
  
  cmap            := new((/14,4/),"float")
  cmap(0,:)        = namedcolor2rgba("gray")
  cmap(1,:)        = namedcolor2rgba("black")
  cmap(2:13,:)     = read_colormap_file("GMT_paired")

;  res_zm@xyLineColors     := cmap((/1,3,5,7,9,11/),:)

;  res_zm@xyLineThicknessF       = .1
;  res_zm@xyLineOpacityF         = .1

  res_zm@xyDashPattern         := 0

  resP_zm              = True
  resP_zm@gsnFrame     = False
  resP_zm@gsnMaximize  = True
  resP_zm@gsnPanelFigureStringsFontHeightF = 0.015
  resP_zm@amJust         = "TopLeft"
  resP_zm@gsnPanelFigureStrings = (/"LP","PETM","EECO"/)
  resP_zm@pmLegendDisplayMode  = "NoCreate"

  plot1  = new(3,graphic)
  plot2  = new(3,graphic)

  ;res_zm@xyExplicitLegendLabels := " gm: "+sprintf("%4.2f", gm(tt,:,0))+" ("+sprintf("%4.2f", gm(tt,:,1))+"-"+sprintf("%4.2f", gm(tt,:,2))+") "+experiments

  lat(0)  = -90
  lat(71) = 90
 
  res_zm@xyLineThicknessF       = .1
  res_zm@xyLineColor           := "gray"
;  res_zm_mean@xyLineColors     := (/cmap(0,:),cmap(2,:),cmap(4,:),cmap(6,:),cmap(8,:)/)
  res_zm_mean@xyLineColors     := (/cmap(0,:)/)
  res_zm_mean@xyLineThicknessF  = 3.
 
  do tt=0,dimsizes(time_slices)-1
    plot1(tt) = gsn_csm_xy(wks,T_sample_zm(tt,0,:,:),lat,res_zm)
    plot2(tt) = gsn_csm_xy(wks,T_zm(tt,0,:),lat,res_zm_mean)
    overlay(plot1(tt),plot2(tt))


  end do

  gsn_panel(wks,plot1,(/1,3/),resP_zm)
end if