pro astrometric

;---------------------------------------------
;This is the control panel, where we toggle things on and off

;You will need to make these directories if they do not exist:  'astro_out' and 'atmospheres_new'

start_image=1     							   ;The starting number of images to run on
end_image=1									;The ending number of images to run on
start_iterate=0									;The starting number of iterations of an image
end_iterate=9									;The ending number of iterations of an image
the_directory=''								;Declares the next value a string
the_directory='examples/'				;The directory where the input catalogs are
prefix='rotate_'										;The prefix of the image name, before the numbers
bin_increment=50								;Size of steps to make in bins of distance


;---------------------------------------------

number_of_columns=(end_iterate-start_iterate)+1
;This code finds the distances between stars in a set of 10 exposures

;LOOP THROUGH ALL SETS OF INPUT IMAGES
for p=start_image, end_image do begin

	;LOOP THROUGH THE 10 EXPOSURES
	for z=start_iterate, end_iterate do begin

	namer=the_directory+prefix+strcompress(p,/remove_all)+'_'+strcompress(z,/remove_all)+'.cat'
	
	;DELETE THE RESULTS FILE SINCE WE USE "APPEND"
	result=file_search(the_directory+'astrometry.'+strcompress(p,/remove_all)+'_'+strcompress(z,/remove_all)+'.zero', count=out1)
	IF out1 GT 0 then begin
	file_delete, the_directory+'astrometry.'+strcompress(p,/remove_all)+'_'+strcompress(z,/remove_all)+'.zero'
	endif

		IF z EQ 0 then begin 
		rdfloat, namer, id, x1, y1, flags, columns=[1,7,8,22], /DOUBLE, skipline=37 
		forprint, id, x1, y1, flags, textout=the_directory+'standard.0', format='(i0,1x,D0.0,1x,D0.0,1x,i0)', /NOCOMMENT
		print, 'Working on file: ', namer
		x=x1
		y=y1
		endif

		IF z GT 0 then begin
		rdfloat, namer, id, x1, y1, flags, columns=[1,7,8,22], /DOUBLE, skipline=37 
		rdfloat, the_directory+'standard.0', id_stan, x_stan, y_stan, flag_stan, /DOUBLE
		print, 'Working on file: ', namer

		x=findgen(N_ELEMENTS(x_stan))
		y=findgen(N_ELEMENTS(x_stan))

			for zango=0, N_ELEMENTS(x_stan)-1 do begin
			diffy_x=(abs(x_stan[zango]-x1))
			diffy_y=(abs(y_stan[zango]-y1))
			diffy=SQRT((diffy_x^2)+(diffy_y^2))
			fin=where(diffy EQ min(diffy))
			x[zango]=x1[fin]
			y[zango]=y1[fin]
			endfor

		endif	

	;SORT THE STARS BY RA SO THEY STAY IN THE SAME ORDER

	dippie=SQRT((x^2)+(y^2))
	w=fsort(dippie, sorted)
	x_ok=x[w]
	y_ok=y[w]
	id_ok=findgen(N_ELEMENTS(x))+1


		;RUN THROUGH EACH STAR, LETTING EACH ONE BE THE ORIGIN
		nardo=N_ELEMENTS(x)
		for i=0, N_ELEMENTS(x)-1 do begin
		openw, lun, the_directory+'astrometry.'+strcompress(p,/remove_all)+'_'+strcompress(z,/remove_all)+'.zero', /get_lun, /append
		
		x0=x_ok[i]
		y0=y_ok[i]
		front=id[i]

		;MEASURE DISTANCE
		dx = double(x_ok-x0)
		dy = double(y_ok-y0)
		dist=SQRT((dx^2)+(dy^2))*0.2
		num1=round(front)
		num2=round(id_ok)
	
		;DO NOT ALLOW REPETITION IN STARS
		fine=where(num1 NE num2 AND num2 GT num1, zippy)

		dist_fine=dist[fine]
		numby=strcompress(round(num1[fine]),/remove_all)+'_'+strcompress(round(num2[fine]),/remove_all)

			IF zippy NE 0 then begin
			;PRINT THE SET OF STAR DISTANCES EACH TIME.  THE IF STATEMENT IS TO DELETE PRINTING THE FINAL DIAGONAL ELEMENT
				for j=0, N_ELEMENTS(dist_fine)-1 do begin
				printf, lun, numby[j], dist_fine[j], format='(a0, 1x, D0.0)' 
				endfor	
			endif 
		close, lun
		free_lun, lun
		endfor
	endfor

;---------------------------------------------------------------

totes=file_lines(the_directory+'astrometry.'+strcompress(p,/remove_all)+'_0.zero')

A=findgen(number_of_columns,totes)

for z=0, number_of_columns-1 do begin

readcol, the_directory+'astrometry.'+strcompress(p,/remove_all)+'_'+strcompress(z,/remove_all)+'.zero', id, dist, format='(a0,D0.0)'

		for i=0, N_ELEMENTS(dist)-1 do begin
		A[z,i]=dist[i]
		endfor

endfor


distance=findgen(totes)
dispersion=findgen(totes)

	for j=0, totes-1 do begin
	distance[j]=MEAN(A[start_iterate:end_iterate,j])
	dispersion[j]=STDDEV(A[start_iterate:end_iterate,j])
	endfor

mad_max=floor(1000./bin_increment)	

meanie=findgen(mad_max)
dispie=findgen(mad_max)

for k=0, mad_max-1 do begin

fine=where(distance GE (k*bin_increment) AND distance LT ((k+1)*bin_increment) AND dispersion LT 0.100)

meanie[k]=MEAN(distance[fine])
dispie[k]=MEAN(dispersion[fine])*1000.


endfor

bad=where(dispersion GT 0.100)

print, 'THIS MANY BAD:', N_ELEMENTS(dispersion[bad]), ' out of ', N_ELEMENTS(dispersion)

forprint, meanie, dispie, textout=the_directory+prefix+'atmosphere_new.'+strcompress(p,/remove_all)+'.zero', /NOCOMMENT, format='(D0.0,1x,D0.0)'

endfor

astrometric_plot

end 
