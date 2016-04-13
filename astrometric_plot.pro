pro astrometric_plot

;SETUP STUFF
!p.multi=[0,1,1,0,0]
A = FINDGEN(17) * (!PI*2/16.)  
USERSYM, COS(A), SIN(A), /FILL  
set_plot, 'PS'
simpctable

DEVICE, /encapsul, /color, /landscape, FILENAME='examples/output.ps'

for j=1, 1 do begin

	rdfloat, 'examples/rotate_atmosphere_new.'+strcompress(j,/remove_all)+'.zero', meanie0, dispie0, /DOUBLE

	If j EQ 1 then begin
	plot, meanie0, dispie0, psym=8, xtitle='Mean separation (arcsec)', ytitle='Dispersion (mas)', title='Astrometric Error vs. Separation', yrange=[0,50], ystyle=1, xrange=[0,800]
	oplot, meanie0, dispie0, linestyle=0
	endif else begin
	oplot, meanie0, dispie0, psym=8, color=!red
	oplot, meanie0, dispie0, linestyle=0, color=!red
	endelse

endfor

device, /close

spawn, 'ps2pdf -dAutoRotatePages=180 examples/output.ps examples/output.pdf'

spawn, 'open -a Adobe\ Reader examples/output.pdf'


end