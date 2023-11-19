## make spg fingerprint (note cdo weights averages by gridcell area in fldmean)

cdo sub -fldmean -sellonlatbox,-55,-20,46,61 -selvar,${var} ${ifile} -fldmean -selvar,${var} ${ifile} ${ofile}

## yearly average weighted by the number of days per month

cdo -L divdpy -yearsum -muldpm ${ifile} ${ofile}

## yearly average of uncertainty file

cdo sqrt -divdpy -yearsum -muldpm -sqr ${ifile} ${ofile}

## select just the Atlantic ocean

cdo sellonlatbox,-90,29,-90,90 -selvar,${var} ${ifile} ${ofile}

## mean of top 300m in the salinity dataset (cdo weights average by level height)

cdo -L vertmean -sellevidx,1/19 ${ifile} ${ofile}

## mean of top 300m for salinity uncertainty (in reality the uncertainties are not independent, but the salinity uncertainty is only used in Figure S5 and not in any analysis)

cdo -L sqrt -vertmean -sqr sellevidx,1/19 ${ifile} ${ofile}

## get just SST from HadCRUT files

cdo -expr,'topo = ((topo<0.0)) ? 1.0 : topo/0.0' -remapcon,${ifile} -topo ${land_mask}
cdo -mul ${ifile} ${land_mask} ${ofile}


