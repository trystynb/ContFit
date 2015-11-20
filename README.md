# ContFit
A spectrum continuum fitter GUI based on drawing continuum with a mouse and fitting a cubic spline.

Author: Trystyn Berg

This readme contains a brief outline of the purpose of the code, the input/output requirements and formats, and a quick start guide to using the code. For a tutorial, please see CONTFIT_tutorial.pdf in the tutorial folder.

This GUI was written for astronomers to perform continuum fitting of 1D QSO spectra for quasara absorption lines.
ContFit could be extended to other continuum fitting purposes. The fit is found by the user clicking with their mouse
along where they think the continuum is. All the click locations are recorded and used to fit a cubic spline to the
continuum. The input object and error spectra (FITS formatted; but can be edited to any format you want) are divided
by the generated continuum, and saved as output. The program also saves the continuum as an ASCII file, as well as
a python binary pickle file with the x/y coordinates of the spline (the mouse clicks).

Python requirements:
TB_CONTFIT was developed with:

		-Python Version 2.7
		-Matplotlib Version 1.4.3
		-TKINTER Revision: 81008
		-Numpy version 1.9.2
		-Scipy version 0.16.0
		-Pyfits version 3.3
		Other versions of above packages may work...

Input file requirements:

    Input spectrum as a FITS file, with flux as data extension and headers CRVAL1, CDELT1, NAXIS1
    to define wavelength. Wavelength can be in log units, but the check for this is based on 
    optical spectra in units of angstroms. See code for more details. 
    
    Error spectrum (not necessary). Must be of same format and length as input spectrum.
    
Output files:

    Object/error spectra continuum normalized in FITS format (as above)
    
    TSV file with wavelength and flux of continuum
    
    Python binary file of the x/y coordinates of the spline (mouse clicks), 
    saved using python's pickle module. in format [XSPLINE,YSPLINE]
    
Quick usage

    Select input spectrum, error spectrum, and name your output files. NOTE: If you are
    tweaking the continuum, make sure you output continuum file matches the name of your
    original file as the software looks for a pickle file of the same name +'.p'. Choose 
    the "chunk" size to view the spectrum one chunk at a time. If you are a first time
    user, use the drop-down menu in the top corner to toggle a tutorial mode on.
    
    Once ready to fit conitnuum (and your wrist is relatively free of carpal tunnel pains...),
    click the Fit continuum button. This will start the continuum fitting procedure. Two 
    windows will pop-up: one showing the full spectrum to guide you, while the other shows the
    "chunk" of spectrum you are fitting. Click along the spectrum in the Continuum window. You
    can only click new points rightwards of the last point. Once finished, click left of the
    first point (the big red circle with a red dashed line) to refresh the window to the next
    chunk. Repeat until all chunks are done.
    
    After fitting, a new window will pop up showing two panels. The top displays the original
    spectrum and what the continuum looks like. The bottom panel shows what the 
    continuum-normalized spectrum looks like. The error spectrum is plotted in green as 1+err
    and 1-err to reflect the possible scatter about the continuum solely due to noise. This
    is only to guide the eye, and means nothing else.
    
    If you are unhappy with parts of the fit; instead of refitting the continuum, the drop-down
    menu has a continuum fitting feature. By using this, you can add/remove points from the
    spline and see how that affects the continuum. Once it is tweaked to your satisfaction,
    click "save spline" button and exit. This will keep the spline for later, as well as update
    your output files.
    
