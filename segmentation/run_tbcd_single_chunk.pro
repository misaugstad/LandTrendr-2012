;**************************************************************************** 
;Copyright © 2008-2011 Oregon State University                                
;All Rights Reserved.                                                         
;                                                                             
;                                                                             
;Permission to use, copy, modify, and distribute this software and its        
;documentation for educational, research and non-profit purposes, without     
;fee, and without a written agreement is hereby granted, provided that the    
;above copyright notice, this paragraph and the following three paragraphs    
;appear in all copies.                                                        
;                                                                             
;                                                                             
;Permission to incorporate this software into commercial products may be      
;obtained by contacting Oregon State University Office of Technology Transfer.
;                                                                             
;                                                                             
;This software program and documentation are copyrighted by Oregon State      
;University. The software program and documentation are supplied "as is",     
;without any accompanying services from Oregon State University. OSU does not 
;warrant that the operation of the program will be uninterrupted or           
;error-free. The end-user understands that the program was developed for      
;research purposes and is advised not to rely exclusively on the program for  
;any reason.                                                                  
;                                                                             
;                                                                             
;IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT, 
;INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST      
;PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN 
;IF OREGON STATE UNIVERSITYHAS BEEN ADVISED OF THE POSSIBILITY OF SUCH        
;DAMAGE. OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES,       
;INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
;FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY WARRANTY OF               
;NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,    
;AND OREGON STATE UNIVERSITY HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,       
;SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.                            
;                                                                             
;**************************************************************************** 

function run_tbcd_single_chunk, info,  $
		subset, index, mask_image, output_image_group, $
		within_layer_offset, layersize, kernelsize, $
		background_val, $
		skipfactor, desawtooth_val, $
		pval, max_segments, normalize , $
		fix_doy_effect, divisor , $
		recovery_threshold, minneeded, $
		distweightfactor,  vertexcountovershoot, $
    	bestmodelproportion, progressbaryesno

	;differs from s3 in that it writes out images according to
	;  abrupt vs. slow disturbance.
	;differs from s4 in that this one eventually calls find_segments6 instead
	;  of find_segments7.  that one was an aborted line in evolution that
	;  used the non-linear fitting.

	;eventually calls find_segments6, which uses quick approach to
	;  id segments.

	;this is intended for raw bands (either 5 or 7) for the sample
	; scenes, but is founded from the biomass fitting programs
	;  that use the segmentation rather than the iterative non-linear
	;  curvefitting.  Also, this will allow two disturbance/recovery sections
	;  subject to some constraints to prevent overfitting.

	;version 3 differs from priors in that here we use only
	;   the hypothesized models, not all possible combos of models.

	;version 22 differs from 2 in that this one uses groups of output iamges
	;   to keep filesize manageable.
	;biomass2 version uses segmentation -- totally different strategy --
	;  that results in different output layers.


	minimum_number_years_needed = minneeded    ;if we have fewer years than this, we can't do it.

	;6/20/08 the number of years is not the number of info items, because now we allow
	;  multiple images per year.  So check on the number of unique years in the
	;  year array

	years = fast_unique(info.year)
	years = years[sort(years)]
	n_yrs = n_elements(years)

	n_images = n_elements(info)	;because we need to go through each image, even if multiple per year

	if n_yrs lt minimum_number_years_needed then begin
		print, 'run_tbcd_single_chunk:  there are fewer than the minimum'
		print, 'number of years available for disturbance/recovery extraction'
		print, 'the minimum is: '+string(minimum_number_years_needed)
		print, 'the number of files given to extract_disturbance_recovery4.pro: '+string(n_yrs)
		print, 'confirm that the information from find_image_stack_files is correct'
		return, {ok:0}
	end

;set up progress variables 

progressval = 0 ;set to record progress
  progressinc = 10  ;for the increment



	;check on the mask image

	if file_exists(mask_image) eq 0 then begin
		print, 'run_tbcd_single_chunk.pro needs to have a mask image.'
		print, 'the mask image should be 0s and 1s, with 1s indicating where
		print, ' to run the curve fitting.
		return, {ok:0}
	end


	;for the first year, just get the full subset, then use that
	;as the template.

	if n_elements(subset) eq 0 then begin
		print, 'run_tbcd_single_chun k needs to have a "subset" keyword set'
		return, {ok:0}
	end



	;check on the mask image

tempsubset = subset
	zot_img, mask_image, mask_hdr, mask_img, subset=tempsubset


	if max(mask_img) gt 1 then begin
		print, 'The mask image must have a maximum value of 1, to indicate
		print, '  where to run the curve-fitting.  The mask image
		print, mask_image
		print, '   has a maximum of '+string(max(mask_img))
		return, {ok:0}
	end


	;  image_file:'', $
	;   			image_path:'', $
	;   			type:0, $				;1 mtbs 2 nonmtbs 3 reference year mtbs
	;   					nbr_file:'', $
	;   					tc_file:'', $
	;   					b6_file:'', $
	;   					year:0, $
	;   					julday:0, $
	;   					unique_year:0, $    ;flagged 1 if only 1 image in this year
	;   					n_in_year:0, $		;number of images in this year
	;   					image_priority:0, $	;priority of picking if more than one image per year
	;   					cloudyear:0, $
	;   					cloud_diff_file:'', $
	;   					shadow_diff_file:'', $
	;   					tc_cloud_diff_file:'', $
	;   					cloud_file:'none', $
	;   					subset:[ [0.d, 0.d],[0.d, 0.d]], $
	;   					useareafile: ''}


	;First, build an image to hold the different years, then read them in

	;use first year as a template
tempsubset=subset

	zot_img, info[0].image_file, hdr, img, subset=tempsubset, layer=[1], /hdronly

	;lcount = n_elements(layer)

	;if lcount eq 0 then layer = 1
	;zot_img, path+file_list[0], hdr, img, subset=subset, layer=layer, /hdronly
	;       if lcount eq 2 then begin
	;           layer2 = layer[1]
	;           layer1 = layer[0]
	;        end else layer1 = layer



	if hdr.pixeltype ne 6 and hdr.pixeltype ne 3 and hdr.pixeltype ne 5 then begin
		print, 'run_tbcd_single_chunk expects the image to be of integer type'
		print, 'this one is type number '+string(hdr.pixeltype)

		return, {ok:0}
	end

	;make up a new image with the right dimensions.
	;the new image could potentially have multiple values for a given year,
	;which will be handled by the cloud mask.
	img = intarr(hdr.filesize[0], hdr.filesize[1], n_yrs)
	cld_img = bytarr(hdr.filesize[0], hdr.filesize[1], n_yrs)		;added v4
	usedmask = intarr(hdr.filesize[0], hdr.filesize[1]) ;valide values for years with multiple image

	;which image was used
	idx_img = bytarr(hdr.filesize[0], hdr.filesize[1], n_yrs)

	;now go through and build it.
	k = 0
	for i = 0, n_yrs-1 do begin
	;for i = 0, n_images-1 do begin
		; zot_img, info[i].image_file, hdr, img1, subset = subset, layer=layer1
		fileid = i + k
		this = where(info.year eq years[i], n)

		;current year
		cur_mask = bytarr(hdr.filesize[0], hdr.filesize[1])
		cur_img = intarr(hdr.filesize[0], hdr.filesize[1])
		cur_idx = bytarr(hdr.filesize[0], hdr.filesize[1])

		;YEARS WITH SINGLE IMAGE

		if n eq 1 then begin
			tempsubset = subset
			landtrendr_image_read, info[fileid], hdr, img1, tempsubset, index, modifier, background_val
			sz = size(img1, /dim)

			;now check vs. background. If so, then assign to the cloud
			;   image, since that's what I check before calling the
			;   fitting algorithm.
			bads = where(img1 eq background_val, n_bads)
			if n_bads ne 0 then cld_img[*,*,i] = (cld_img[*,*,i]+ (img1 eq background_val)) ne 0 	;ne 0 needed incase cloud image and background val!

			if n_elements(sz) gt 2 then begin
				print, 'run_tbcd_single_chunks: each image layer must have a single layer'
				print, 'image '+file_list[fileid]+ 'has more than 1 layer'
				return, {ok:0}
			end

			img[*,*,i] = img1/divisor	;added 2/7/08 this will scale to max of 1000

			idx_img[*,*,i] = replicate(fileid, size(img1, /dim))

			;now read the cloud mask
			; if there is no cloud mask, then just skip this
			if info[fileid].cloud_file ne 'none' and info[fileid].cloud_file ne '' then begin
				tempsubset=subset
				if info[fileid].cloud_file eq 'band8' then $
					zot_img, info[fileid].image_file, clhdr, mimg, layers=[8], subset=tempsubset else $
					zot_img, info[fileid].cloud_file, clhdr, mimg, subset=tempsubset
				cld_img[*,*,i] = (cld_img[*,*,i] + (mimg eq 0)) ne 0
			;  cld_img[*,*,i] = (cld_img[*,*,i] + (img1 ne 0)) ne 0 ;0 is no-cloud in cheng's cloudmasks
			end
		end

		;MULTIPLE IMAGES PER YEAR
		;if multiple image exists for this year, select one and make the others masked out
		if n gt 1 then begin
			victims = info[this]
			;sort by priority
			vicorder = sort(victims.image_priority)
			victims = victims[vicorder]
			;read in the cloud image
			for j = 0, n-1 do begin
				tempsubset=subset

				landtrendr_image_read, victims[j], hdr, img1, tempsubset, index, modifier, background_val

				sz = size(img1, /dim)
				if n_elements(sz) gt 2 then begin
					print, 'run_tbcd_single_chunks: each image layer must have a single layer'
					print, 'image '+file_list[fileid+j]+ 'has more than 1 layer'
					return, {ok:0}
				end

				;now read the cloud mask
				; if there is no cloud mask, then just skip this
;print, victims[j].cloud_file

				mimg = replicate(0, size(img1, /dim))
				if victims[j].cloud_file ne 'none' and victims[j].cloud_file ne '' then begin
						tempsubset=subset
						if victims[j].cloud_file eq 'band8' then $
						zot_img, victims[j].image_file, clhdr, mimg, layers=[8], subset=tempsubset else $
						zot_img, victims[j].cloud_file, clhdr, mimg, subset=tempsubset
						;cld_img[*,*,fileid+j] = (cld_img[*,*,fildid+j] + (mimg gt 2300)) ne 0
				end


				;identify pixels that are not background, that haven't been picked by the higher
				;   priority image, and that are not in the cloud mask

				valid = where(img1 ne background_val and cur_mask eq 0 and mimg eq 1, n_valid)
				if n_valid ne 0 then begin
					cur_img[valid] = img1[valid]
					cur_mask[valid] = 1					;mask gets set to 1 if the pixel is chosen
					cur_idx[valid] = replicate(this[vicorder[j]], n_valid)
				end
			endfor
			k = k + n - 1
			img[*,*,i] = cur_img/divisor
			cld_img[*,*,i] = cur_mask ne 1		;any cur_mask pixels still remaining 0 were not chosen
			idx_img[*,*,i] = cur_idx
		end
	end

	img1 = 0 ;reset to save space
	cur_img = 0
	cur_mask = 0
	;
	;    ;write out the stacked image
	;     outfile = strmid(write_file, 0, strlen(write_file)-4)+'_stack.bsq'
	;
	;     openw, un, outfile, /get_lun
	;     writeu, un, img
	;     free_lun, un
	;     hdr.n_layers = n_yrs
	;     write_im_hdr, outfile, hdr
	;

	;observe year stuff

	sz = size(img, /dim)
	;      n_yrs = sz[2]
	;        x_axis = indgen(n_yrs)

	;v4 has the actual years, offset by the min

	min_year = min(info.year)
	;  x_axis = info.year-min_year	+ 1		;set up by the year, but set to 1 to be consistent with
	;fitting functions


	x_axis = years		;these were "uniqued" early on, so should be okay.
	;x_axis = info.year

	;x_axis = info.year		;changed on 9/5 to match with


	;set up the progress bar:
	progressbaryesno = 0 ; when parallelized, we are not using prgoress bar

	if progressbaryesno eq 1 then begin
		progressBar = Obj_New("PROGRESSBAR", /fast_loop, title = 'Curve-fitting:  percent done')
		progressBar -> Start
	end

	; Instead of initializing here, we are initializing as shared memory
	; vertyear_image = intarr(sz[0], sz[1], output_image_group[0].n_layers)
	; vertvals_image = intarr(sz[0], sz[1], output_image_group[1].n_layers)
	; mag_image = intarr(sz[0], sz[1],output_image_group[2].n_layers)
	; dur_image = intarr(sz[0], sz[1], output_image_group[3].n_layers)
	; distrec_image = intarr(sz[0], sz[1], output_image_group[4].n_layers)

	; fitted_image = intarr(sz[0], sz[1], output_image_group[5].n_layers)
	; stats_image = intarr(sz[0], sz[1], output_image_group[6].n_layers)

	; segmse_image = intarr(sz[0], sz[1], output_image_group[7].n_layers)
	; source_image = intarr(sz[0], sz[1], output_image_group[8].n_layers)
	; segmean_image = intarr(sz[0], sz[1], output_image_group[9].n_layers)

	totalcount = float(sz[0]*sz[1])

	ksq=kernelsize^2
	seed= randomseed()
	if n_elements(skipfactor) eq 0 then skipfactor = 3

	offset = (kernelsize-1)/2

	; If our current chunk contains anything in the usearea mask, then we split
	; the data by rows and run in parallel in child processes. Otherwise, we
	; skip this process altogether.
	; The progress bar feature is currently not implemented in parallel version.
	if TOTAL(mask_img) > 0 then begin
		nCPUs = !CPU.HW_NCPU-1

		; map some anonymous shared memory into into our address space
		SHMMAP, 'seg_cld_img', /INTEGER, DIMENSION=[SIZE(cld_img, /DIMENSIONS)]
		SHMMAP, 'seg_img', /INTEGER, DIMENSION=SIZE(img, /DIMENSIONS)
		SHMMAP, 'seg_mask_img', /INTEGER, DIMENSION=SIZE(mask_img, /DIMENSIONS)
		SHMMAP, 'seg_fitted_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[5].n_layers]
		SHMMAP, 'seg_source_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[8].n_layers]
		SHMMAP, 'seg_vertyear_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[0].n_layers]
		SHMMAP, 'seg_vertvals_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[1].n_layers]
		SHMMAP, 'seg_segmse_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[7].n_layers]
		SHMMAP, 'seg_segmean_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[9].n_layers]
		SHMMAP, 'seg_mag_image', /INTEGER, DIMENSION=[sz[0], sz[1],output_image_group[2].n_layers]
		SHMMAP, 'seg_stats_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[6].n_layers]
		SHMMAP, 'seg_dur_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[3].n_layers]
		SHMMAP, 'seg_distrec_image', /INTEGER, DIMENSION=[sz[0], sz[1], output_image_group[4].n_layers]

		; make variables pointing to our different segments in shared memory
		temp_cld_img = SHMVAR('seg_cld_img')
		temp_img = SHMVAR('seg_img')
		temp_mask_img = SHMVAR('seg_mask_img')
		temp_fitted_image = SHMVAR('seg_fitted_image')
		temp_source_image = SHMVAR('seg_source_image')
		temp_vertyear_image = SHMVAR('seg_vertyear_image')
		temp_vertvals_image = SHMVAR('seg_vertvals_image')
		temp_segmse_image = SHMVAR('seg_segmse_image')
		temp_segmean_image = SHMVAR('seg_segmean_image')
		temp_mag_image = SHMVAR('seg_mag_image')
		temp_stats_image = SHMVAR('seg_stats_image')
		temp_dur_image = SHMVAR('seg_dur_image')
		temp_distrec_image = SHMVAR('seg_distrec_image')

		; initialize the shared memory with the images that we have in local mem
		temp_cld_img[0, 0, 0] = cld_img
		temp_img[0, 0, 0] = img
		temp_mask_img[0, 0, 0] = mask_img
		; initialize the rest of the shared memory variables w/ default values
		temp_fitted_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[5].n_layers)
		temp_source_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[8].n_layers)
		temp_vertyear_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[0].n_layers)
		temp_vertvals_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[1].n_layers)
		temp_segmse_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[7].n_layers)
		temp_segmean_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[9].n_layers)
		temp_mag_image[0, 0, 0] = intarr(sz[0], sz[1],output_image_group[2].n_layers)
		temp_stats_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[6].n_layers)
		temp_dur_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[3].n_layers)
		temp_distrec_image[0, 0, 0] = intarr(sz[0], sz[1], output_image_group[4].n_layers)

		; This is the procedure call that does the work inside the loop. It is
		; in a string so the command can be sent to the child processes.
		proc_str = 'segment_pixels, start_x, end_x, offset, skipfactor, cld_img, img, ' $
			+ 'n_fix_doy_effect, pval, minneeded, seed, ' $
			+ 'minimum_number_years_needed, desawtooth_val, background_val, ' $
			+ 'recovery_threshold, max_segments, distweightfactor, ' $
			+ 'vertexcountovershoot, bestmodelproportion, modifier, mask_img, ' $
			+ 'fitted_image, source_image, vertyear_image, vertvals_image, ' $
			+ 'juldays, segmse_image, segmean_image, mag_image, stats_image, ' $
			+ 'dur_image, idx_img, x_axis, distrec_image'

		; count how many indices are computed and how many for each CPU.
		if sz[0] > 2*offset then xstodo = 1 + (sz[0] - 2*offset - 1)/skipfactor else xstodo = 0
		xs_per_cpu = xstodo / nCPUs

		oBridge = objarr(nCPUs-1)
		for i=0, nCPUs-1 do begin
			; Compute start and end rows for each process, parent process gets
			; all the leftovers from integer division.
			start_x = i * xs_per_cpu * skipfactor + offset
			end_x = start_x + (xs_per_cpu - 1) * skipfactor
	        if i EQ n_elements(oBridge) then end_x = (sz[0]-(offset+2))

	        ; If the parent process, no need to transfer data, just do the work.
	        if i EQ nCPUs-1 then begin
	            segment_pixels, start_x, end_x, offset, skipfactor, temp_cld_img, temp_img, $
	            	n_elements(fix_doy_effect), pval, minneeded, seed, $
	            	minimum_number_years_needed, desawtooth_val, background_val, $
	            	recovery_threshold, max_segments, distweightfactor, vertexcountovershoot, $
	            	bestmodelproportion, modifier, temp_mask_img, temp_fitted_image, $
	            	temp_source_image, temp_vertyear_image, temp_vertvals_image, info.julday, $
	            	temp_segmse_image, temp_segmean_image, temp_mag_image, temp_stats_image, $
	            	temp_dur_image, idx_image, x_axis, temp_distrec_image

	        ; if child process, transfer data and start the work
	        endif else begin
	            oBridge[i] = obj_new('IDL_IDLBridge') ; start up the child process

	            ; pass data (none in this first set are modified)
	            oBridge[i]->SetVar, 'start_x', start_x
	            oBridge[i]->SetVar, 'end_x', end_x
	            oBridge[i]->SetVar, 'offset', offset
	            oBridge[i]->SetVar, 'skipfactor', skipfactor
	            oBridge[i]->SetVar, 'n_fix_doy_effect', n_elements(fix_doy_effect)
	            oBridge[i]->SetVar, 'pval', pval
	            oBridge[i]->SetVar, 'minneeded', minneeded
	            oBridge[i]->SetVar, 'seed', seed
	            oBridge[i]->SetVar, 'minimum_number_years_needed', minimum_number_years_needed
	            oBridge[i]->SetVar, 'desawtooth_val', desawtooth_val
	            oBridge[i]->SetVar, 'background_val', background_val
	            oBridge[i]->SetVar, 'recovery_threshold', recovery_threshold
	            oBridge[i]->SetVar, 'max_segments', max_segments
	            oBridge[i]->SetVar, 'distweightfactor', distweightfactor
	            oBridge[i]->SetVar, 'vertexcountovershoot', vertexcountovershoot
	            oBridge[i]->SetVar, 'bestmodelproportion', bestmodelproportion
	            oBridge[i]->SetVar, 'modifier', modifier
	            oBridge[i]->SetVar, 'x_axis', x_axis
	            oBridge[i]->SetVar, 'juldays', info.julday
	            oBridge[i]->SetVar, 'idx_img', idx_img

	            ; pass dimensions of images in shared memory
	            oBridge[i]->SetVar, 'fitted_dim', SIZE(temp_fitted_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'source_dim', SIZE(temp_source_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'vertyear_dim', SIZE(temp_vertyear_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'vertvals_dim', SIZE(temp_vertvals_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'segmse_dim', SIZE(temp_segmse_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'segmean_dim', SIZE(temp_segmean_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'mag_dim', SIZE(temp_mag_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'stats_dim', SIZE(temp_stats_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'distrec_dim', SIZE(temp_distrec_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'dur_dim', SIZE(temp_dur_image, /DIMENSIONS)
	            oBridge[i]->SetVar, 'cloud_dim', SIZE(temp_cld_img, /DIMENSIONS)
	            oBridge[i]->SetVar, 'img_dim', SIZE(temp_img, /DIMENSIONS)
	            oBridge[i]->SetVar, 'mask_dim', SIZE(temp_mask_img, /DIMENSIONS)

	            ; map addresses of images in anon shared memory into each child's address space.
	            oBridge[i]->Execute, "SHMMAP, 'seg_cld_img', /INTEGER, DIMENSION=cloud_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_img', /INTEGER, DIMENSION=img_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_mask_img', /INTEGER, DIMENSION=mask_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_fitted_image', /INTEGER, DIMENSION=fitted_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_source_image', /INTEGER, DIMENSION=source_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_vertyear_image', /INTEGER, DIMENSION=vertyear_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_vertvals_image', /INTEGER, DIMENSION=vertvals_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_segmse_image', /INTEGER, DIMENSION=segmse_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_segmean_image', /INTEGER, DIMENSION=segmean_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_mag_image', /INTEGER, DIMENSION=mag_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_stats_image', /INTEGER, DIMENSION=stats_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_distrec_image', /INTEGER, DIMENSION=distrec_dim"
	            oBridge[i]->Execute, "SHMMAP, 'seg_dur_image', /INTEGER, DIMENSION=dur_dim"

	            ; Create local variables in each child process pointing to  images in shared memory.
	            oBridge[i]->Execute, "cld_img = SHMVAR('seg_cld_img')"
	            oBridge[i]->Execute, "img = SHMVAR('seg_img')"
	            oBridge[i]->Execute, "mask_img = SHMVAR('seg_mask_img')"
	            oBridge[i]->Execute, "fitted_image = SHMVAR('seg_fitted_image')"
	            oBridge[i]->Execute, "source_image = SHMVAR('seg_source_image')"
	            oBridge[i]->Execute, "vertyear_image = SHMVAR('seg_vertyear_image')"
	            oBridge[i]->Execute, "vertvals_image = SHMVAR('seg_vertvals_image')"
	            oBridge[i]->Execute, "segmse_image = SHMVAR('seg_segmse_image')"
	            oBridge[i]->Execute, "segmean_image = SHMVAR('seg_segmean_image')"
	            oBridge[i]->Execute, "mag_image = SHMVAR('seg_mag_image')"
	            oBridge[i]->Execute, "stats_image = SHMVAR('seg_stats_image')"
	            oBridge[i]->Execute, "distrec_image = SHMVAR('seg_distrec_image')"
	            oBridge[i]->Execute, "dur_image = SHMVAR('seg_dur_image')"

	            ; Compile the function that processes pixels and run it. The nowait keyword means
	            ; this child process doesn't have to finish executing before moving on, so work can
	            ; be done in parallel.
	            oBridge[i]->Execute, ".r segment_pixels.pro"
	            oBridge[i]->Execute, proc_str, /nowait
	        endelse
	    endfor

	    ; Wait for all child processes to finish executing...
	    notdone = 1
	    while notdone do begin
	        done=0
	        for i=0, n_elements(oBridge)-1 do $
	            done = done+oBridge[i]->Status()
	        if done EQ 0 then notdone=done
	    endwhile

	    ; Once everyone is done, stop the child processes.
	    for i=0, n_elements(oBridge)-1 do begin
	        obj_destroy, oBridge[n_elements(oBridge)-1-i]
	    endfor

	    ; Copy data from shared memory over into local memory so it can be released.
	    ; TODO Remove this copying. It is not necessary. Right now we are doing it because string
	    ;      literals are being passed to SHMMAP, and if we don't copy the data over, then the
	    ;      mapped memory's name clashes when you get to a 2nd chunk. We can make it more general
	    ;      where SHMMAP creates a unique name for us, which we can save in a variable, instead
	    ;      of using the same string literal every time.
		vertyear_image = temp_vertyear_image
		vertvals_image = temp_vertvals_image
		mag_image = temp_mag_image
		dur_image = temp_dur_image
		distrec_image = temp_distrec_image
		fitted_image = temp_fitted_image
		stats_image = temp_stats_image
		segmse_image = temp_segmse_image
		source_image = temp_source_image
		segmean_image = temp_segmean_image

	    ; Then make the addresses no longer shared.
		SHMUNMAP, 'seg_cld_img'
		SHMUNMAP, 'seg_img'
		SHMUNMAP, 'seg_mask_img'
		SHMUNMAP, 'seg_fitted_image'
		SHMUNMAP, 'seg_source_image'
		SHMUNMAP, 'seg_vertyear_image'
		SHMUNMAP, 'seg_vertvals_image'
		SHMUNMAP, 'seg_segmse_image'
		SHMUNMAP, 'seg_segmean_image'
		SHMUNMAP, 'seg_mag_image'
		SHMUNMAP, 'seg_stats_image'
		SHMUNMAP, 'seg_dur_image'
		SHMUNMAP, 'seg_distrec_image'

	endif else begin
		; If we don't spawn other processes, at least initialize in private memory.
		vertyear_image = intarr(sz[0], sz[1], output_image_group[0].n_layers)
		vertvals_image = intarr(sz[0], sz[1], output_image_group[1].n_layers)
		mag_image = intarr(sz[0], sz[1],output_image_group[2].n_layers)
		dur_image = intarr(sz[0], sz[1], output_image_group[3].n_layers)
		distrec_image = intarr(sz[0], sz[1], output_image_group[4].n_layers)

		fitted_image = intarr(sz[0], sz[1], output_image_group[5].n_layers)
		stats_image = intarr(sz[0], sz[1], output_image_group[6].n_layers)

		segmse_image = intarr(sz[0], sz[1], output_image_group[7].n_layers)
		source_image = intarr(sz[0], sz[1], output_image_group[8].n_layers)
		segmean_image = intarr(sz[0], sz[1], output_image_group[9].n_layers)
	endelse


	if progressbaryesno eq 1 then begin 
	progressBar -> Destroy
	progressBar = Obj_New("PROGRESSBAR", /fast_loop, title = 'Interpolating:  percent done')
	progressBar -> Start
	end else print
;set up progress variables
  progressval = 0 ;set to record progress
  progressinc = 10  ;for the increment


	;now interpolate
	desired_kernel_size = 15
	ks = min([desired_kernel_size, sz[0], sz[1]])		;make sure that the kernel size is not
	desired_kernel_size = ks
	;bigger than the size of the chunk

	;make a distance matrix that can be clipped to all
		ok = make_distance_grid([(ks * 2)+1, (ks * 2)+1], start = [desired_kernel_size, desired_kernel_size])
					;this is twice as big as anything we'll get, so we can then subset it.

 	 	master_geo_dist = (ok.matrix / (desired_kernel_size/2)) + 1	;scale so penalty is 2 for things kernel_size dist away.



	dmat = fltarr(ks,ks)
	halfval = (ks-1)/2
	for x = 0, sz[0]-1 do begin			;start and end one pixel in,
		for y = 0, sz[1]-1 do begin
			checkval = (stats_image[x,y,4] ne 1) + $
				(mask_img[x,y] eq 1)


			if checkval eq 2 then begin

				;first, get cloud info for the desired pixel
				goods_pix = cld_img[x,y,*] ne 1

				;then calc the range of neighborhood pixels

				start_x = max([x-halfval, 0])		;first set up starting point, make sure in image
				start_y = max([y-halfval, 0])
				end_x = min([start_x+ks-1, sz[0]-1])  ;if we're near the other side, bump
				end_y = min([start_y+ks-1, sz[1]-1])
				start_x = end_x-ks+1								 ;jiggle start if we had to bump
				start_y = end_y-ks+1							 ;won't affect anything if we're not near the end

				;set the source image to img[x,y,combined_goods].
				;  this is not interpolated.

				wh_goods_pix = where(goods_pix eq 1, n_goods_pix)
				if n_goods_pix ne 0 then source_image[x,y,wh_goods_pix] = img[x,y,wh_goods_pix]

;w, 0

				;make a penalty score matrix from the geo distance
				geo_offset_x = 	desired_kernel_size-(x-start_x)
				geo_offset_y = desired_kernel_size-(y-start_y)

				penalty_geo_dist = master_geo_dist[ $
						geo_offset_x:geo_offset_x+desired_kernel_size-1, $
						geo_offset_y:geo_offset_y+desired_kernel_size-1]



		count=0		;then go through and get distance in spectral/temporal space

				for i = start_x, end_x do begin
					for j= start_y, end_y do begin



						;if this is a real curve-fitted pixel, then see
						;  how far away

						if stats_image[i,j,4] eq 1 then begin
							goods_test = cld_img[i,j,*] ne 1
							combined_goods = where(goods_test*goods_pix eq 1, ngds)
							if ngds ne 0 then $
								;						  	      dmat[i-start_x, j-start_y] = $
								;						  					sqrt(total(img[i,j,combined_goods] - $
								;						  					img[x,y,combined_goods])^2) else dmat[i-start_x, j-start_y] = 2e32
								dmat[i-start_x, j-start_y] = $
								total(abs(img[i,j,combined_goods] - $
								img[x,y,combined_goods])) else dmat[i-start_x, j-start_y] = 2e32

;							if ngds ne 0 and x gt 17 and y gt 2 then begin
;							    ;w, count, 400,400
;							    if count eq 0 then plot, img[x,y, combined_goods]
;							    oplot, img[i,j, combined_goods] , color = i*12455 + j*945454
;							    xyouts, .1, .7, string(dmat[i-start_x, j-start_y]), /norm
;							    count=count+1
;							    ;a = get_kbrd()
;							end

;plot, img[i,j, combined_goods], psym = 4
;oplot, img[x,y, combined_goods], color = '444499'xl
;print, dmat[i-start_x, j-start_y]
;qqw = get_kbrd()



						end else dmat[i-start_x, j-start_y] = 2e32		;just set to an absurdly high number

					end
				end


				;penalize appropriately
				decents = where(dmat ne 2e32, n_decents)
				if n_decents gt 0 then 	dmat[decents] = dmat[decents] * penalty_geo_dist[decents]

				;then pick closest one

				closest = where(dmat eq min(dmat))
				closest = closest[0]
				pos = getxy(closest, ks, ks)

;stop
				xoffset = (pos[0]+start_x)-x
				yoffset = (pos[1]+start_y)-y
				match_pos_x = x+xoffset
			    match_pos_y = y+yoffset



;if start_x gt 0 and start_y gt 0 then stop


				if min(dmat) ne 2e32 then begin
;					w, count+1
;					plot, img[x,y,*], thick = 2
;					oplot, img[x,y,*], thick = 4, color = '00ff00'xl
;					oplot, img[match_pos_x, match_pos_y, *], thick = 4
;					xyouts, .6, .2, string(dmat(closest)), /norm
;if  xoffset eq 3 and yoffset eq 4 and x gt 17 and y gt 2 then stop
;
;;						a = get_kbrd()

					vertyear_image[x,y,*] = vertyear_image[match_pos_x,match_pos_y, *]
					vertvals_image[x,y,*] = vertvals_image[match_pos_x,match_pos_y, *]
					segmse_image[x,y,*] = segmse_image[match_pos_x,match_pos_y, *]
					segmean_image[x,y,*] = segmean_image[match_pos_x,match_pos_y, *]

					mag_image[x,y,*] = mag_image[match_pos_x,match_pos_y, *]
					distrec_image[x,y,*] = distrec_image[match_pos_x,match_pos_y, *]
					dur_image[x,y,*] = dur_image[match_pos_x,match_pos_y, *]

					fitted_image[x,y,*] = fitted_image[match_pos_x,match_pos_y, *]
					stats_image[x,y,0:3] = stats_image[match_pos_x,match_pos_y,0:3]
					stats_image[x,y,4] = 2	;interpolated
					stats_image[x,y,5] = stats_image[match_pos_x,match_pos_y, 5]
   					stats_image[x,y,6] = stats_image[match_pos_x,match_pos_y, 6]
   					stats_image[x,y,7] = stats_image[match_pos_x,match_pos_y, 7]
   					stats_image[x,y,8] = xoffset
   					stats_image[x,y,9] = yoffset


;					vertyear_image[x,y,*] = vertyear_image[pos[0]+start_x, pos[1]+start_y, *]
;					vertvals_image[x,y,*] = vertvals_image[pos[0]+start_x, pos[1]+start_y, *]
;					segmse_image[x,y,*] = segmse_image[pos[0]+start_x, pos[1]+start_y, *]
;					segmean_image[x,y,*] = segmean_image[pos[0]+start_x, pos[1]+start_y, *]
;
;					mag_image[x,y,*] = mag_image[pos[0]+start_x, pos[1]+start_y, *]
;					distrec_image[x,y,*] = distrec_image[pos[0]+start_x, pos[1]+start_y, *]
;					dur_image[x,y,*] = dur_image[pos[0]+start_x, pos[1]+start_y, *]
;
;					fitted_image[x,y,*] = fitted_image[pos[0]+start_x, pos[1]+start_y, *]
;					stats_image[x,y,0:3] = stats_image[pos[0]+start_x,pos[1]+start_y,0:3]
;					stats_image[x,y,4] = 2	;interpolated
;					stats_image[x,y,5] = stats_image[pos[0]+start_x, pos[1]+start_y, 5]
;   					stats_image[x,y,6] = stats_image[pos[0]+start_x, pos[1]+start_y, 6]
;   					stats_image[x,y,7] = stats_image[pos[0]+start_x, pos[1]+start_y, 7]


				end else begin 			;if no good vals
					vertyear_image[x,y,*] = -1
					vertvals_image[x,y,*] = -1
					segmse_image[x,y,*] = -1
					segmean_image[x,y,*] = -1
					mag_image[x,y,*] = -1
					distrec_image[x,y,*] = -1
					dur_image[x,y,*] = -1


					fitted_image[x,y,*] = -1
					stats_image[x,y,*] = [1.0, 0, 0, 0, 2, -1, 0, 0, 0, 0]	;set p to 1.0, and set to interpolated
				end

			end	;checkval okay
		end	;y


		percent_done = ((float(x)*y)/ totalcount)*100
		test = round((percent_done) / progressinc)
    if test gt progressval then begin
      print, progressval*progressinc  ;only print if we've bumped to next increment
      progressval = test
    end

		if progressbaryesno eq 1 then 	progressBar -> Update, percent_done
		
	
	end	;x
	if progressbaryesno eq 1 then 	progressBar->Destroy

	;write 'em out

	;vertices


	openu, un, output_image_group[0].filename, /get_lun
	for layercount = 0ll, output_image_group[0].n_layers-1 do begin
		point_lun, un, (output_image_group[0].layersize * $
			layercount)+within_layer_offset
		writeu, un, vertyear_image[*,*,layercount]
	end
	free_lun, un

	;vertvals

	openu, un, output_image_group[1].filename, /get_lun
	for layercount = 0ll, output_image_group[1].n_layers-1 do begin
		point_lun, un, (output_image_group[1].layersize * $
			layercount)+within_layer_offset
		writeu, un, vertvals_image[*,*,layercount]*modifier		;added modifier july 9 2008 so values make sense
	end
	free_lun, un


	;mag image

	openu, un, output_image_group[2].filename, /get_lun
	for layercount = 0ll, output_image_group[2].n_layers-1 do begin
		point_lun, un, (output_image_group[2].layersize * $
			layercount)+within_layer_offset
		writeu, un, mag_image[*,*,layercount]
	end
	free_lun, un

	;duration image


	openu, un, output_image_group[3].filename, /get_lun
	for layercount = 0ll, output_image_group[3].n_layers-1 do begin
		point_lun, un, (output_image_group[3].layersize * $
			layercount)+within_layer_offset
		writeu, un, dur_image[*,*,layercount]
	end
	free_lun, un

	;distrec image


	openu, un, output_image_group[4].filename, /get_lun
	for layercount = 0ll, output_image_group[4].n_layers-1 do begin
		point_lun, un, (output_image_group[4].layersize * $
			layercount)+within_layer_offset
		writeu, un, distrec_image[*,*,layercount]
	end
	free_lun, un
	;fitted

	openu, un, output_image_group[5].filename, /get_lun
	for layercount = 0ll, output_image_group[5].n_layers-1 do begin
		point_lun, un, ulong64(output_image_group[5].layersize) * $
			layercount+within_layer_offset
		writeu, un, fitted_image[*,*,layercount]
	end
	free_lun, un

	;stats

	openu, un, output_image_group[6].filename, /get_lun
	for layercount = 0ll, output_image_group[6].n_layers-1 do begin
		point_lun, un, (output_image_group[6].layersize * $
			layercount)+within_layer_offset
		writeu, un, stats_image[*,*,layercount]
	end
	free_lun, un

	;segment mse

	openu, un, output_image_group[7].filename, /get_lun
	for layercount = 0ll, output_image_group[7].n_layers-1 do begin
		point_lun, un, (output_image_group[7].layersize * $
			layercount)+within_layer_offset
		writeu, un, segmse_image[*,*,layercount]
	end
	free_lun, un

	;source image

	openu, un, output_image_group[8].filename, /get_lun
	for layercount = 0ll, output_image_group[8].n_layers-1 do begin
		point_lun, un, (output_image_group[8].layersize * $
			layercount)+within_layer_offset
		writeu, un, source_image[*,*,layercount]
	end
	free_lun, un

	;segmean image

	openu, un, output_image_group[9].filename, /get_lun
	for layercount = 0ll, output_image_group[9].n_layers-1 do begin
		point_lun, un, (output_image_group[9].layersize * $
			layercount)+within_layer_offset
		writeu, un, segmean_image[*,*,layercount]
	end
	free_lun, un




	return, {ok:1}

end

