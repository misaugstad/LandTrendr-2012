PRO segment_pixels, start_x, end_x, offset, skipfactor, cld_img, img, n_fix_doy_effect, pval, $
        minneeded, seed, minimum_number_years_needed, desawtooth_val, $
        background_val, recovery_threshold, max_segments, distweightfactor, $
        vertexcountovershoot, bestmodelproportion, modifier, mask_img, $
        fitted_image, source_image, vertyear_image, vertvals_image, juldays, $
        segmse_image, segmean_image, mag_image, stats_image, dur_image, $
        idx_img, x_axis, distrec_image

    sz = size(img, /dim)
    ; for x = offset, sz[0]-(offset+1), skipfactor do begin
    ;     for y = offset, sz[1]-(offset+1), skipfactor do begin

    for x = start_x, end_x + 1, skipfactor do begin
        for y = offset, sz[1]-(offset+1), skipfactor do begin

            ;check on the mask image to see if we should run this pixel
            ; oh no load balancing
            if mask_img[x,y] eq 1 then begin
                ;check for clouds
                chunk = img[x-offset:x+offset, y-offset:y+offset, *]
                usable = cld_img[x-offset:x+offset, y-offset:y+offset, *] eq 0

                slice = total(chunk*usable,1)
                slice_usable = total(usable, 1)

                vals = total(slice,1)/total(slice_usable, 1)

                goods= where(cld_img[x,y,*] ne 1, ngds)
                if ngds gt minimum_number_years_needed then begin


                    ;first check to see if fix the doy effect
                    if n_elements(fix_doy_effect) ne 0 then begin

                        idxs = idx_img[x,y,*]

                        uniques = fast_unique(juldays[idxs[goods]])
                        if n_elements(uniques) gt 4 then begin
                            r = poly_fit(juldays[idxs[goods]], vals[goods],2, chisq=chisq,yfit = yfit)
                            m = mean(yfit)
                            zzz = calc_fitting_stats3(vals[goods], yfit, 3, resid=resid)
                            if zzz.p_of_f lt pval then outvals = m+resid else $
                                outvals = vals[goods]
                        end else outvals = vals[goods]
                    end else outvals = vals[goods]
        ;feb 2013
                    ; dampen the first and last years
                    dampen = 0.7
                    
                    diffend = outvals[ngds-1]-outvals[ngds-2]
                    outvals[ngds-1] = outvals[ngds-2]+ ((1-dampen)*diffend)
                    
                    diffbeg = outvals[0]-outvals[1]
                    outvals[0]=outvals[1]+((1-dampen)*diffbeg)
            

                    ok=fit_trajectory_v2(x_axis,goods, outvals, $
                        minneeded, background_val, $
                        modifier, seed, $
                        desawtooth_val, pval, $
                        max_segments, recovery_threshold, $
                        distweightfactor,  vertexcountovershoot, $
                        bestmodelproportion)


                    if ok.ok eq 1 then begin
                        ;take out the bad year
                        ;fitted_image[x,y,*] = round(ok.best_model.yfit[uniq(info.year)])   ;all years, including masked out, will get fittedvals
                        fitted_image[x,y,*] = round(ok.best_model.yfit) ;prior version's line (commented, above) resulted in flatline at end
                                                ;because of illogic -- all of the yfits belong in the fitted image
;plot, x_axis[goods], outvals, psym = 4
;oplot, x_axis, round(ok.best_model.yfit), color = '444499'xl
;qqw = get_kbrd()


                        source_image[x,y,goods] = outvals   ;



                        vertyear_image[x,y,*] = ok.best_model.vertices      ;these are true years
                        vertvals_image[x,y,*] = ok.best_model.vertvals      ;in the units fed to fit_trajectory_v1
                        segmse_image[x,y,*] = ok.best_model.segment_mse     ;mse of each segment
                        for ss = 0, ok.best_model.n_segments-1 do $         ;mean of each segment
                                segmean_image[x,y,ss] = (vertvals_image[x,y,ss]+vertvals_image[x,y,ss+1])/2.

                        ;get the magnitudes and the proportions
                        temp = shift(ok.best_model.vertvals, -1) - ok.best_model.vertvals
                        mag_image[x,y,0:ok.best_model.n_segments-1] = temp[0:ok.best_model.n_segments-1]

                        maxdist = max(mag_image[x,y,0:ok.best_model.n_segments-1], min=maxrec)
                        distrec_image[x,y, 0]=max([maxdist,0])
                        distrec_image[x,y, 1]=min([maxrec, 0])



                        totalmag = total(abs(mag_image[x,y, *]))    ;the total distance traversed, up or down
                        summag = float(total(mag_image[x,y, *]))            ;the actual value with pluses and minuses
                        if totalmag eq 0 then distrec_image[x,y, 2] = (-1500) else $
                            distrec_image[x,y, 2] = (summag/totalmag)*1000  ;will be -1000 if all rec, + 1000 if all dist

                        ;get the durations
                        temp = shift(ok.best_model.vertices, -1) - ok.best_model.vertices

                        dur_image[x,y,0:ok.best_model.n_segments-1] = temp[0:ok.best_model.n_segments-1]


                        ;                   mag_image[x,y,*] =
                        ;                    = intarr(sz[0], sz[1], max_segments)
                        ;           distrec_image = intarr(sz[0], sz[1], max_segments)




                        if ok.best_model.f_stat gt 300 then ok.best_model.f_stat = 300

                        stats_image[x,y,0] = round(ok.best_model.p_of_f*100)
                        stats_image[x,y,1] = round(ok.best_model.f_stat*100)
                        stats_image[x,y,2] = round(ok.best_model.ms_regr/10.)
                        stats_image[x,y,3] = round(ok.best_model.ms_resid/10.)
                        stats_image[x,y,4] = 1          ;directly run?
                        stats_image[x,y,5] = ok.best_model.n_segments
                        stats_image[x,y,6] = x_axis[goods[0]]   ;set to minimum usable year
                        stats_image[x,y,7] = n_elements(goods)  ;number of usable years
                        ;stats_image layers 8 and 9 are used only for interpolation, giving the offset of the pixel
                    end
                end
            end
        end ;y
    end ;x
    ; now return all the data that was modified
end