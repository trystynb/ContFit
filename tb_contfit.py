#!/usr/bin/python2.7

""" 
TB_CONTFIT.PY
V1.1 (06/11/2015)
	-Fixed zoom feature when add/remove points in continuum tweaker.

Author: Trystyn Berg (trystynb@uvic.ca)


SUMMARY AND USAGE:
	Stand-alone python application to fit continuum to spectrum.
	
	Load a spectrum from fits file, and fit a cubic spline to the spectrum.
	The spline points are chosen by the user clicking along where they think
	the continuum is (this was inspired by DIPSO continuum fitting). 

	The continuum fitting process is done by splitting the spectrum into
	"spectrum chunks", the size of which (in wavelength space) is defined
	by the user. Each chunk is displayed while the user clicks the along the
	continuum. At the end, the user clicks left of the first point (red circle)
	to move to next "chunk".

	Once finished, the user can inspect the fit via a sequence of plots. In
	addition, a continuum tweaker function can be used to add/remove spline
	points, and see their effect on the continuum.

	Once finished either task, the input spectrum and error spectrum are divided
	by the fitted conitnuum and saved as fits files. The continuum itself is saved
	as an ascii file. The spline points (mouse click locations) are saved in a
	python binary file for internal application and future use.
	

PYTHON PACKAGE REQUIRMENTS:
	TB_CONTFIT was developed with:
		-Python Version 2.7
		-Matplotlib Version 1.4.3
		-TKINTER Revision: 81008
		-Numpy version 1.9.2
		-Scipy version 0.16.0
		-Pyfits version 3.3
	Other versions of above packages may work...

FITS FILE FORMAT (INPUT AND OUTPUT)
	Flux in the first data extension
	Contains 'CRVAL1', 'CDELT1','NAXIS1' headers to generate wavelength
	Can be either in log or linear scale (assuming the wavelength is in angstroms, and is
	not <10\AA... See SPECFITS function for more details.

	Feel free to modify the fits reader/write functions (SPECFITS/SAVEFITS) to your 
	appropriate format.


INPUT FILE REQUIREMENTS
	-FITS file containing spectrum for continuum fitting.
	-FITS file containing error spectrum for continuum fitting.

OUTPUT FILES
	-FITS file containing continuum-fitted spectrum
	-FITS file containing continuum-fitted error spectrum
	-ASCII file containing the continuum generated.
		Each line in file must be in the format (whitespace delimeted):
			WAVELENGTH FLUX
	-Python binary file containing XSPLINE/YSPLINE lists from continuum fit. This
		is a pickle file storing these values so you can regernate the interpolation
		of the spline. When reading in, it is of the format [XSPLINE,YSPLINE].

KNOWN ISSUES:
	-In the continuum tweaking function, sometimes the REMOVE POINT feature is either slow or non-responsive.
		Keep trying it will eventually work. Also, if multiple points exist near the one you are trying to 
		remove, zoom in.

"""

######################
###START OF PROGRAM###
######################

#Import all required packages
import matplotlib.backends.backend_tkagg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import numpy as np
import math
import os,sys
import Tkinter,tkFileDialog, tkMessageBox
import scipy.interpolate
import cPickle
import pyfits

""" 
Miscellanious functions

SORTLIST(xs,ys) takes a list of x/y coordinates (XS/YS) and sorts them in numerical order of x values

SAVESPLINE(xspline,yspline,wvlngth,outfile) interpolates the continuum using XSPLINE/YSPLINE information,
	and generates the continuum for the wavlength points WVLNGTH. The WVLNGTH/continuum flux is 
	save to OUTFILE as a tab-separated ASCII file, while the XSPLINE/YSPLINE lists are "pickled"
	into a python binary file for future use.

SPECFITS(INP, LOG=FALSE) reads the INP fits file (see above for required format).It should
	be smart enough to tell the difference between a log and linear scale; but do be careful
	if you are using non-optical wavelength data... If so, set LOG=TRUE/FALSE accordingly

SAVEFITS(NEWFLUX, INFITS, OUTFITS) takes the flux array DATA, and saves in over the original flux
	values in INFITS. The output fits file is called OUTFITS. This may seem strange, but it 
	preserves the header/wavelength information, while only changing the flux.

"""



def SortList(xs,ys):
	#XS - The x-coordinate list
	#YS - The y-coordinate list

	#Convert XS/YS into arrays to use with numpy (NXS/NYS)
	nxs=np.array(xs)
	nys=np.array(ys)
	#Sort the XS array, and get the corresponding indicies of the sort
	xind=np.argsort(nxs)
	#Sort NXS/NYS using the indicies
	xs=nxs[xind]
	ys=nys[xind]
	#Return sorted  XS/YS as lists
	return list(xs), list(ys)


def SaveSpline(xspline,yspline,wvlngth,outfile):
	#XSPLINE - The x-coordinates for the spline fit of the continuum
	#YSPLINE - The y-coordinates for the spline fit of the continuum
	#WVLNGTH - The wavelength array of the spectrum to generate a continuum
	#OUTFILE - The name of the output ascii file of the continuum

	#Interpolate a cubic spline from XSPLINE/YSPLINE
	spline=scipy.interpolate.interp1d(xspline,yspline,kind='cubic')
	#Make spline "pickle" file name
	splinefile=outfile+'.p'
	#Save the XSPLINE/YSPLINE information to the binary python file SPLINEFILE
	cPickle.dump([xspline,yspline],open(splinefile,'wb'))
	#Generate the continuum from the interpolated spline
	cont=spline(wvlngth)
	#Write it out to a tab-separated file OUTFILE
	f=open(outfile,'w')
	f.write('#Wvlngth\tFlux\n')
	for ii in range(len(wvlngth)):
		f.write('%.7f\t%s\n'%(wvlngth[ii],str(cont[ii]).format(precision=7,type='e')))
	f.close()

#Function to load fits file
def specfits(inp, log=False):
	hdulist=pyfits.open(inp)#Opens 1D spectrum
	spectrum=hdulist[0].data#Read 1D spectrum into array
	wstart=hdulist[0].header['CRVAL1']#Convert starting wavelength into angstroms
	dw=hdulist[0].header['CDELT1']#convert wavelength/pixel scale into angstroms
	npix=hdulist[0].header['NAXIS1']#Get the number of pixels
	#Check to see if pixel scale is log or linear
	#My trick here is that if the first wavleneght value is less than 5 units,
	#it is likely a log scale. HOWEVER, if your spectrum actually is a linear scale,
	#you might want to comment this out.... A warning will be printed if it does this...
	if wstart<5:
		log=True
		print "WARNING: TB_CONTFIT - SPECFITS() set wavelength scale to LOG"

	hdulist.close()#Close fits image
	#Generate the spectrum axis of the 1D spectrum
	wvlngth=np.zeros(npix)
	for i in range(npix):
		if log:	wvlngth[i]=10**(wstart+i*dw)
		else: wvlngth[i]=wstart+i*dw
	#return an array of the wavelength and flux of the spectrum.
	return wvlngth, spectrum

#Function to rewrite a spectrum (to preserve header information, but just change flux values)
#Please be carefule if you use this code for something else...
def savefits(newflux, infits ,outfits):
	#Open original files
	hdulist=pyfits.open(infits)
	#Load the original flux
	fitsdat=hdulist[0].data
	#Replace the flux wit the new flux
	for ii in range(len(newflux)):
		fitsdat[ii]=newflux[ii]
	#Create new file. If it already exists, warn the user.
	hdulist.writeto(outfits, output_verify='warn', clobber=True)
	#Close the fits image
	hdulist.close()



""" 
CLASS: LINEBUILDER

This is the heart of the continuum fitting routine. It takes mouse 
clicks on the matplotlib continuum fitting window and builds the 
spline. It will draw ontop of the appropriate display windows
where you have clicked, and connect the points with a blue line to show the 
rough continuum.


Call linebuilder=LineBuilder(LINE,XMIN,XMAX,SPECPLOT,FULLAX,FULLSPECPLOT):

INPUTS: LINE - a line artist in matplotlib. This contains the spline points that you draw
	XMIN - The minimum wavelength in the continuum fitting window. Used to check if click
		within wavelength range
	XMAX - The maximum wavelength in the continuum fitting window. Used to check if click
		withing wavelength range
	SPECPLOT - The TKINTER matplotlib canvas (the continuum fitting window)
	FULLAX - The matplotlib Axes instance for the entire spectrum plot 
		in FULLSPECPLOT (display purposes)
	FULLSPECPLOT - The TKINTER matplotlib canvas showing the entire spectrum (display purposes)
	

ATTRIBUTES:
	linebuilder.LINE - The LINE input
	linebuilder.XMIN - The supplied XMIN value
	linebuilder.XMAX - The supplied XMAX value
	linebuilder.XS - a list of all the wavelengths of the spline (where the mouse was clicked)
	linebuilder.YS - a list of all the fluxes of the spline (where the mouse was clicked)
	linebuilder.XSTART - The starting wavelength to build the spline drawing;
		reference for exiting class
	linebuilder.SPECPLOT - The SPECPLOT input
	linebuilder.FULLSPECPLOT - The FULLSPECPLOT input
	linebuilder.FULLAX - The FULLAX input
	linebuilder.CID - The Event ID number for drawing continuum
	

NOTES
	The important thing to get from this are the XS and YS for the spline.
	It is a modified version from the matplotlib event handling tutorial.

	To exit LINEBUILDER events, you need to click behind the last point of the
	spline from the previous iteration (i.e. linebuilder.XSTART)

"""

class LineBuilder:
	#Initialize when called
	def __init__(self1, line, xmin,xmax,SpecPlot,fullax,FullSpecPlot):
		self1.line = line #Keeps track of the line artist being drawn
		#Get the X/Y data of the line artist for internal use
		self1.xs = list(line.get_xdata())
		self1.ys = list(line.get_ydata())
		#Get the minimum/maximum wavelength for the display window 
		self1.xmin=xmin
		self1.xmax=xmax
		#The starting point of the line being drawn. It is needed
		#to keep track of when the user is done (click behind line)
		self1.xstart=self1.xs[-1]
		#Need to keep track of the two display windows
		self1.SpecPlot=SpecPlot#display window for drawing the continuum
		#display for seeing the full spectrum while doing a piece by piece fit
		self1.FullSpecPlot=FullSpecPlot
		self1.fullax=fullax#MPL.axes instance of full spectrum window
		#Start up the mouse click event to start drawing the line, and leave running
		self1.cid = self1.line.figure.canvas.mpl_connect('button_press_event', self1)
		self1.line.figure.canvas.start_event_loop(0)
		return

	def __call__(self1, event):
		#Check to make sure that the click is within the axes of the window
        	if event.inaxes!=self1.line.axes: return
		#Get the x/y coordiantes of the mouse click (in wavelength-flux space)
		x=event.xdata
		y=event.ydata
		#Make sure the points are within the minimum/maximum value for recording
		#Also, make sure the point is beyond the last click.
		if x>=self1.xmin and x<=self1.xmax and x>self1.xs[-1]:
			#Append the coordinates to the appropriate lists of clicks
		        self1.xs.append(x)
		        self1.ys.append(y)
			#Update the plotting windows to draw on where the
			#mouse was clicked (full spectrum) and blue line
			#connecting clicks in continuum window.
			self1.fullax.plot(x,y,'or')
			self1.FullSpecPlot.draw()
		        self1.line.set_data(self1.xs, self1.ys)
		        self1.line.figure.canvas.draw()
		#if user has clicked behind the first point of the fit, end the event loop
		elif x<=self1.xstart:
			self1.line.figure.canvas.stop_event_loop()
		return
""" 
CLASS: CONTTWEAKER

This is a feature available if there are bits in the continuum you fit that you aren't happy
with and need to refit/tweak parts without doing the whole continuum again. Here you can
keep track how changing (removing/adding) the spline points affects the continuum, and compare
to your original fit.

Call CT.ContTweaker(IWVLNGTH,ISPECTRUM,ESPECTRUM,OLDXSPLINE,OLDYSPLINE,OUTFILE)

INPUTS:
	IWVLNGTH - The wavelength of the input spectrum (that you fitted)
	ISPECTRUM - The original spectrum (i.e. before continuum fit)
	ESPECTRUM - The original error spectrum (bfore dividing by continuum)
	OLDXSPLINE - The original x-coordinates (wavelength) for the spline of the continuum (from LINEBUILDER)
	OLDYSPLINE - The original y-coordinates (flux) for the spline of the continuum (from LINEBUILDER)
	OUTFILE - The name of the ouput file you are saving the continuum fitted spectrum to.

ATTRIBUTES:
	CT.OUTFILE - Same as input OUTFILE
	CT.IWVLNGTH - Same as input IWVLNGTH
	CT.ISPECTRUM - Same as input ISPECTRUM
	CT.ESPECTRUM - Same as input ESPECTRUM
	CT.XSPLINE - The x-coordinates (wavelength) for genereating the spline
	CT.YSPLINE - The y-coordinates (flux) for genereating the spline
	CT.OLDXSPLINE - The input OLDXSPLINE (for reference)
	CT.OLDYSPLINE - The input OLDYSPLINE (for reference)
	CT.OLDSPLINE - The interpolated spline instance from CT.OLDXSPLINE/CT.OLDYSPLINE
	CT.SPLINE - The interpolated spline from CT.XSPLINE/CT.YSPLINE
	CT.POPUP - The  root (extension) of the buttons used to control the tweaker int he Tkinter widget

	CT.TWEAKFIG - The matplotlib figure to display the spectrum while tweaking
	CT.TFROOT - The Tkinter widget for the CT.TWEAKFIG display
	CT.TFPLOT - The Tkinter canvas to display the CT.TWEAKFIG (embedded in CT.TFROOT)
	CT.AX1 - The MPL axes in CT.TWEAKFIG that will display the original continuum and spectrum
	CT.AX2 - The MPL axes in CT.TWEAKFIG that will display the current spline points/continuum and spectrum
	CT.AX3 - The MPL axes in CT.TWEAKFIG that displays the newest continuum-fitted spectrum
	CT.YMIN - The minimum y-value to set the CT.AX? windows on CT.refresh()
	CT.YMIN - The maximum y-value to set the CT.AX? windows on CT.refresh()

	CT.CIDBUT - The MPL event id for clicking the mouse (button)
	CT.CIDPICK - The MPL event id for picking an artist element in the axes.

INTERNAL FUNCTIONS:
	CT.MAKEPLOT() - Generates the plotting window (CT.TFPLOT) for displaying the plots for tweaking
	CT.REFRESH(YREFRESH=TRUE) - Refreshes the plotting window (CT.TFPLOT) in a particular way
	CT.ADDPOINT() - Run a MPL event to add a spline point, and run CT.CLICKADD()	
	CT.REMOVEPOINT() - Run a MPL event to remove a spline point, and run CT.CLICKREMOVE()
	CT.QUITEDIT(CID) - Quit an MPL event (with id CID)
	CT.SAVE() - Save the spline to a file
	CT.EXIT() - Stop the continuum fitter
	CT.CLICKADD() - Add the x/y coordinates of the mouse event to CT.XSPLINE/CT.YSPLINE
	CT.CLICKREMOVE() - Find the associated spline point with the mouse event and remove from CT.XSPLINE/CT.YSPLINE

Notes:
	-CT.CLICKREMOVE is sometimes a bit buggy. Not 100% what the problem is.
	-You CANNOT remove the first/last spline point. if you want to get rid of it, generate
		a point before/after it.

"""
class ContTweaker:
	def __init__(self,iwvlngth,ispectrum,espectrum,oldxspline,oldyspline,outfile):
		#Initialize input variables for the class
		self.outfile=outfile
		self.iwvlngth=iwvlngth
		self.ispectrum=ispectrum
		self.espectrum=espectrum
		#intialize CT.XSPLINE/YSPLINE to the old values at first (will be replaced)
		self.xspline=oldxspline
		self.yspline=oldyspline
		self.oldxspline=oldxspline
		self.oldyspline=oldyspline
		#Interpolate cubic spline of the continuum for both "NEW" and old spline
		self.oldspline=scipy.interpolate.interp1d(oldxspline,oldyspline,kind='cubic')
		self.spline=scipy.interpolate.interp1d(oldxspline,oldyspline,kind='cubic')
		#Make the MPL figures appear
		self.MakePlot()
		#Add buttons to the bottom of the main TB_CONTFIT widget menu
		self.popup=Tkinter.Frame()
		self.popup.grid()
		But1=Tkinter.Button(self.popup,text='Refresh Plot',command=self.Refresh)
		But1.grid(column=0,row=0)
		But3=Tkinter.Button(self.popup,text='Remove Point',command=self.RemovePoint)
		But3.grid(column=0,row=1)
		But4=Tkinter.Button(self.popup,text='Add Point',command=self.AddPoint)
		But4.grid(column=0,row=2)
		But5=Tkinter.Button(self.popup,text='Reopen Window',command=self.MakePlot)
		But5.grid(column=0,row=4)		
		But6=Tkinter.Button(self.popup,text='Save Spline',command=self.Save)
		But6.grid(column=0,row=5)
		But7=Tkinter.Button(self.popup,text='Exit',command=self.Exit)
		But7.grid(column=0,row=6)


		
		#In development: Continous editing rather than having to click
		""" 
		MODES=[("Single/Stop",False),("Continuous",True)]
		self.useContinuous=Tkinter.BooleanVar()
		self.useContinuous.set(False)
		ii=1
		for text,mode in MODES:
			b=Tkinter.Radiobutton(self.popup,text=text,variable=self.useContinuous,value=mode)
			b.grid(column=ii,row=1)
			ii+=1
			if ii==1: b.select()
		"""

		#Wait until the window is closed (i.e. CT.EXIT() is run by clicking EXIT button)
		self.popup.wait_window()




	def MakePlot(self):
		#This function generates the plot window that will be the vidual interface for how
		#the continuum tweaking is going

		#Make the Spline from the CURRENT spline values
		xcont=self.iwvlngth
		ycont=self.spline(xcont)
		#Make the spline from the OLD spline values
		ospectrum=self.oldspline(xcont)

		#Create a MPL figure
		self.tweakfig=plt.figure(figsize=(8,8))
		#First subplot is for displaying the original fit to the spectrum
		self.ax1=self.tweakfig.add_subplot(3,1,1)
		#Second subplot is for displaying the current fit to the spectrum
		#To help with looking at all regions of the same time, share the X and Y
		#axes to CT.AX1.
		self.ax2=self.tweakfig.add_subplot(3,1,2,sharex=self.ax1,sharey=self.ax1)
		#Third subplot for displaying the continuum-fitted spectrum based on the
		#current spline. Because this is now 0 and 1ish, only share the X axis
		#to CT.AX1
		self.ax3=self.tweakfig.add_subplot(3,1,3,sharex=self.ax1)
		#Create the TK window for embedding the figure, set it up, and draw
		self.TFroot=Tkinter.Tk()
		self.TFroot.wm_title("Continuum Tweaker")
		self.TweakPlot=FigureCanvasTkAgg(self.tweakfig,master=self.TFroot)
		#This toolbar is mainly for zooming. By sharing the x (and y) axes, when zooming in
		#on one axes, all three will zoom appropriately
		NavTweakPlot=NavigationToolbar2TkAgg(self.TweakPlot, self.TFroot)
		self.TweakPlot.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

		#For CT.AX1, plot the original specturm, the error spectrum, and the original continuum fitted
		self.ax1.plot(self.iwvlngth,self.ispectrum,'k',drawstyle='steps',label='Spectrum')
		self.ax1.plot(self.iwvlngth,self.espectrum,'g',drawstyle='steps',label='Error')
		self.ax1.plot(xcont,ycont,'b',label='Continuum')
		self.ax1.plot(self.oldxspline,self.oldyspline,'or',label='Spline Pts.')
		#Make a legend
		self.ax1.legend(ncol=4,frameon=False, loc=9, bbox_to_anchor=(0.5, 1.3))
		self.ymin,self.ymax=self.ax1.get_ylim()
		self.ax1.set_ylabel('Flux\n(original)')

		#For CT.AX2, plot the spectrum, and the current spline points and the resulting continuum
		self.ax2.plot(self.iwvlngth,self.ispectrum,'k',drawstyle='steps')
		self.ax2.plot(xcont,ycont,'b')
		self.ax2.plot(self.xspline,self.yspline,'or',picker=5)#Picker needed to pick which point
		self.ax2.set_ylabel('Flux\n(tweaked fit)')



		#For CT.AX3, plot the spectrum upon dividing by the continuum. Plot the
		#error spectrum as 1-err and 1+err to reflect how much noise one MIGHT expect
		#As a rule of thumb, one should aim to have the continuum fluctuations within the error
		self.ax3.plot(xcont,self.ispectrum/ospectrum,'k',drawstyle='steps')
		self.ax3.plot(xcont,1.0-self.espectrum/ospectrum,'--g',drawstyle='steps')
		self.ax3.plot(xcont,1.0+self.espectrum/ospectrum,'--g',drawstyle='steps')
		self.xmin,self.xmax=self.ax3.get_xlim()
		self.ax3.plot([self.xmin,self.xmax],[1,1],'--r')
		self.ax3.plot([self.xmin,self.xmax],[0,0],'--r')
		#restrict CT.AX3 to only show the normalized range of values
		self.ax3.set_ylim(-1,2)
		self.ax3.set_ylabel('Relative flux')
		self.ax3.set_xlabel('Wavlength')

		#Draw all the new plotted stuff
		self.TweakPlot.draw()



	def Refresh(self,yrefresh=True):
		#This function takes any modifications that have taken place and updates
		#the plots accordingly. If YREFRESH=TRUE, this sets the y-axis to the
		#original scaling. Otherwise keep the current values.

		#Get the current x values (incase the user has zoomed in using toolbar)
		xmin,xmax=self.ax2.get_xlim()
		ymin,ymax=self.ax2.get_ylim()
		#if YREFRESH=TRUE, set to original values (SELF.YMIN/YMAX)
		if yrefresh:
			ymin=self.ymin
			ymax=self.ymax
		#Clear CT.AX2/AX3 of all previous information
		self.ax2.clear()
		self.ax3.clear()
		#Remake the spline based on the new XSPLINE/YSPLINE points, and save
		self.spline=scipy.interpolate.interp1d(self.xspline,self.yspline,kind='cubic')
		#Generate the continuum based on the new spline
		xcont=self.iwvlngth
		ycont=self.spline(xcont)

		#Update CT.AX1, perserve the x-axis bounds, and return to original y-axis
		self.ax1.set_ylim(ymin,ymax)
		self.ax1.set_xlim(xmin,xmax)

		#Update CT.AX2 by plotting new spline, also perserve the x-axis bounds
		self.ax2.plot(self.iwvlngth,self.ispectrum,'k',drawstyle='steps')
		self.ax2.plot(xcont,ycont,'b')
		self.ax2.plot(self.xspline,self.yspline,'or',picker=5)
		self.ax2.set_ylabel('Flux\n(tweaked fit)')
		self.ax2.set_ylim(ymin,ymax)
		self.ax2.set_xlim(xmin,xmax)

		#in CT.AX2, Divide out the spectrum&errorspectrum by the continuum, and plot
		self.ax3.plot(xcont,self.ispectrum/ycont,'k',drawstyle='steps')
		self.ax3.plot(xcont,1.0-self.espectrum/ycont,'--g',drawstyle='steps')
		self.ax3.plot(xcont,1.0+self.espectrum/ycont,'--g',drawstyle='steps')
		self.ax3.set_ylabel('Relative flux')
		self.ax3.set_xlabel('Wavlength')
		self.ax3.set_ylim(-1,2)
		self.ax3.plot([self.xmin,self.xmax],[1,1],'--r')
		self.ax3.plot([self.xmin,self.xmax],[0,0],'--r')
		self.ax3.set_xlim(xmin,xmax)

		#Update plotting window
		self.TweakPlot.draw()


	#Function for closing the current MPL event by it's ID, and stop the event loop
	#It might seem redundant that I have both things, but I intend to have a continous
	#editing method, which would need the looper.
	def QuitEdit(self,cid):
		#Close event ID
		self.TweakPlot.mpl_disconnect(cid)
		#Stop event loop
		self.TweakPlot.stop_event_loop()
	#Function when "Add Point" button is clicked.
	def AddPoint(self):
		#Show Tutorial message for what to do.
		if usetutorial: tkMessageBox.showinfo("Help Message", "Click where to add point.")
		#Start mouse click event, and run CT.ClickAdd
		self.cidbut=self.TweakPlot.mpl_connect('button_press_event',self.ClickAdd)
		self.TweakPlot.start_event_loop(0)
	#Given a mouse event for adding a point...	
	def ClickAdd(self,event):
		#Grab the x/y coordiantes of the click, and add to spline
		self.xspline.append(event.xdata)
		self.yspline.append(event.ydata)
		#Sort the spline data to be in order by wavelength
		self.xspline,self.yspline=SortList(self.xspline,self.yspline)
		#Refresh the plot with new data, but keep y-axis
		self.Refresh(yrefresh=False)
		#Close the MPL event stuff
		self.QuitEdit(self.cidbut)
	#Function ro remove a point when "Remove Point" button pressed
	def RemovePoint(self):
		#Show tutorial message on what to do
		if usetutorial: tkMessageBox.showinfo("Help Message", "Click point to remove.")
		#Start MPL event for picking an MPL artist, and start the loop. Run CT.ClickRemove
		self.cidpick=self.TweakPlot.mpl_connect('pick_event',self.ClickRemove)
		self.TweakPlot.start_event_loop(0)
	#Given a picker event for removing a point...
	def ClickRemove(self,event):
		#Get the spline point that you picked, it's x and y coordinates
		splinepoint = event.artist
		xsplineval=splinepoint.get_xdata()
		ysplineval=splinepoint.get_ydata()
		#Index of the artist
		ind = event.ind
		#Make sure the point is in the spline point lists
		if xsplineval[ind] in self.xspline:
			if ysplineval[ind] in self.yspline:
				#Remove that point from the spline, I think this is where sorting is important...
				self.xspline.pop(ind)
				self.yspline.pop(ind)
		#Refresh the plot with new spline, but keep y-axis
		self.Refresh(yrefresh=False)
		#Close the event and stop the event loop
		self.QuitEdit(self.cidpick)
	#Function saves the spline using SaveSpline function
	def Save(self):
		print "Saving Spline"
		SaveSpline(self.xspline,self.yspline,self.iwvlngth,self.outfile)
	#Destroy CT.POPUP and CT.TFROOT (the masters for the CT buttons and plots) and leave CT.
	def Exit(self):
		self.popup.destroy()
		self.TFroot.destroy()
		return


##############################
###INTIALIZATION OF PROGRAM###
##############################

#If you choose, you can run TB_CONTFIT from the command line with inputs loaded from command line

#Get input from command line.
inputs=sys.argv
#If you want to have popup messages telling you how to use TB_CONTFIT, set to True
usetutorial=False

#Default values for input into TB_CONTFIT widget

#Initial spectrum (to be continuum fitted) path/file
initinspec=''
#Initial output spectrum (that is continuum fitted) path/file
initoutspec=''
#Initial error spectrum path/file
initespec=''
#Initial continuum fitting chunk size
initlen='500'
#initial name of the continuum file generated
initcont='mycont.dat'


#If, from command line, you type:
#	tb_contfit 'initspec' 'initoutspec' 'initespec' 'initlen' 'initcont'
#It will save these files into the default values to load. Otherwsie, you will have to set them manually.
if len(sys.argv)==6:
	initinspec=sys.argv[1]
	initoutspec=sys.argv[2]
	initespec=sys.argv[3]
	initlen=sys.argv[4]
	initcont=sys.argv[5]


""" 
MAIN GUI APPLICATION. This widget houses the main interface for loading the inputouput files,
doing the continuum fit, and allowing the user to tweak the continuum.


ATTRIBUTES:

	ENTRYIN - The string containing the path/filename of the input spectrum for continuum fitting
	ENTRYERR - The string containing the path/filename of the error spectrum
	ENTRYCONT - The string containing the path/filename of the output continuum ASCII file
	ENTRYOUT - The string containing the path/filename of the output continuum-fitted spectrum
	ENTRYXR - The string containing the size of the wavelength chunk to fit the continuum by

	FULLSPECFIG - The MPL figure for the full spectrum plot
	FULLSPECPLOT - The MPL canvas containing FULLSPECFIG
	FFROOT - The TK window to embed FULLSPECPLOT in
	FULLAX - The Axes instance for plotting the full spectrum in FULLSPECFIG

	SPECFIG - The MPL figure for the chunk of the spectrum to be continuum fitted
	SPECPLOT - The MPL canvas containing SPECFIG
	SFROOT - The TK window to embed SPECPLOT in
	AX - The Axes instance for plotting the spectrum chunk SPECFIG
	
	
INTERNAL FUNCTIONS:
	ONEXIT() - Quits TB_CONTFIT
	INITIALIZE() - Initializes the TB_CONTFIT main widget
	GETINFILE() - Use a TK file finder to get the input spectrum file (ENTRYIN) 
	GETERRFILE() - Use a TK file finder to get the error spectrum (ENTRYERR)
	GETCONTFILE() - Use a TK file finder to get the output continuum file (ENTRYCONT)
	GETOUTFILE() - Use a TK file finder to get the output continuum-fitted spectrum (ENTRYOUT)
	ONTUTORIAL() - Function to toggle tutorial message boxes on/off
	ONTWEAK() - Function to startup the Tweaker class and save resulting file
	ONFITCONT() - Function to startup the LineBuilder class and save the resulting continuum files

"""
#Define the main class for the widget
class simpleapp_tk(Tkinter.Tk):# Base for main window
        #This defines the main window
        def __init__(self,parent):
                Tkinter.Tk.__init__(self,parent)
                self.parent =parent
                self.initialize()

	#Function for killing TB_CONTFIT when Exit button clicked in dropdown menu
        def onExit(self):
                print "Quitting..."
                self.quit()
	#Initialize the main TKinter window interface for TB_CONTFIT
        def initialize(self):
		#Display a tutorial message telling the user what to do when first starting.
		if usetutorial: tkMessageBox.showinfo("Help Message", "Select input spectrum, error spectrum, and output continuum files."+\
			"The X-RANGE field provides the size of each wavelength chunk you want to fit.")
		#Set up a grid for the widget.
                self.grid()
                #Set-up drop down menu with buttons to tweak continuum, toggle tutorial mode, and exit TB_CONTFIT
                mb=Tkinter.Menubutton(self,text='Menu',relief="raised")
                mb.grid(column=0,row=0, sticky='EW')
                picks=Tkinter.Menu(mb,tearoff=0)
		#Set to run ONTWEAK for Tweak Continuum button
                picks.add_command(label="Tweak Continuum",command=self.onTweak)
		#Set to run ONTUTORIAL for Tutorial toggle button
		picks.add_command(label="Tutorial mode on/off",command=self.onTutorial)
		#Set to run ONEXIT to kill TB_CONTFIT on exit
                picks.add_command(label="Exit",command=self.onExit)
                mb.config(menu=picks)

                ##############################
                #Define the coordinate fields#
                ##############################

                #The label for the input spectrum field
                labelIn=Tkinter.StringVar()
                labelin=Tkinter.Label(self,textvariable=labelIn,\
                        anchor="w",fg="black")
                labelin.grid(column=0, row=1, sticky='EW')
                labelIn.set(u"Input spectrum:")
                #The field for the input spectrum field (ENTRYIN). Initialize to default INITINSPEC
                self.entryIn=Tkinter.StringVar()
                entryin=Tkinter.Entry(self,textvariable=self.entryIn)
                entryin.grid(column=1,row=1,stick='EW')
                self.entryIn.set(initinspec)
		#Search file button for finding the Input spectrum in a window (runs GETINFILE())
		sinbutton=Tkinter.Button(self,text="Search",command=self.getInFile)
		sinbutton.grid(column=2,row=1, sticky='EW')



                #The label for the input error spectrum field
                labelErr=Tkinter.StringVar()
                labelerr=Tkinter.Label(self,textvariable=labelErr,\
                        anchor="w",fg="black")
                labelerr.grid(column=0, row=2, sticky='EW')
                labelErr.set(u"Error spectrum:")
                #The field for the input error spectrum field (ENTRYERR). Initialize to default INITESPEC
                self.entryErr=Tkinter.StringVar()
                entryerr=Tkinter.Entry(self,textvariable=self.entryErr)
                entryerr.grid(column=1,row=2,stick='EW')
                self.entryErr.set(initespec)
		#Search file button for finding error spectrum in window (GETERRFILE)
		serrbutton=Tkinter.Button(self,text="Search",command=self.getErrFile)
		serrbutton.grid(column=2,row=2, sticky='EW')
		
                #The label for the output spectrum file field
                labelOut=Tkinter.StringVar()
                labelout=Tkinter.Label(self,textvariable=labelOut,\
                        anchor="w",fg="black")
                labelout.grid(column=0, row=3, sticky='EW')
                labelOut.set(u"Output spectrum file (FITS):")
                #The field for the output sprectrum file (ENTRYOUT). Initialize to INITOUTSPEC
                self.entryOut=Tkinter.StringVar()
                entryout=Tkinter.Entry(self,textvariable=self.entryOut)
                entryout.grid(column=1,row=3,stick='EW')
                self.entryOut.set(initoutspec)
		#Search file button for finding directory/file of output spectrum file (runs GETOUTFILE)
		soutbutton=Tkinter.Button(self,text="Search",command=self.getOutFile)
		soutbutton.grid(column=2,row=3, sticky='EW')


                #The label for the output continuum file field
                labelCont=Tkinter.StringVar()
                labelcont=Tkinter.Label(self,textvariable=labelCont,\
                        anchor="w",fg="black")
                labelcont.grid(column=0, row=4, sticky='EW')
                labelCont.set(u"Output Continuum file (ASCII):")
                #The field for the output continuum file (ENTRYCONT). Initialize to INITCONT
                self.entryCont=Tkinter.StringVar()
                entrycont=Tkinter.Entry(self,textvariable=self.entryCont)
                entrycont.grid(column=1,row=4,stick='EW')
                self.entryCont.set(initcont)
		#Search file button for finding directory/file of output continuum file (runs GETCONTFILE)
		sContbutton=Tkinter.Button(self,text="Search",command=self.getContFile)
		sContbutton.grid(column=2,row=4, sticky='EW')


                #The field label for the wavelength chunk for continuum fitting.
                labelXR=Tkinter.StringVar()
                labelxr=Tkinter.Label(self,textvariable=labelXR,\
                        anchor="w",fg="black")
                labelxr.grid(column=0, row=5, sticky='EW')
                labelXR.set(u"Wavelength Chunk Size:")
                #The field for the wavelength chunk (ENTRYXR). Initalize to default INITLEN
                self.entryXR=Tkinter.StringVar()
                entryxr=Tkinter.Entry(self,textvariable=self.entryXR)
                entryxr.grid(column=1,row=5,stick='EW')
                self.entryXR.set(initlen)

                #Define button for Fitting conitnuum 
                plotbutton=Tkinter.Button(self,text=u'Fit Continuum',\
                        command=self.OnFitCont)#Runs ONFITCONT function when clicked
                plotbutton.grid(column=0,row=6,columnspan=2,sticky='EW')

	#TK widget for finding the input spectrum. Save it to ENTRYIN when found
	def getInFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current spectrum file
		browser_opt['initialfile']=self.entryIn.get()
		#Label the window
		browser_opt['title']='Select Input spectrum'
		#Run the dialog, and save the output to the spectrum file name ENTRYIN
		self.entryIn.set(tkFileDialog.askopenfilename(**browser_opt))
		return
	
	#TK widget for finding the error spectrum. Save it to ENTRYERR when found
	def getErrFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current spectrum file
		browser_opt['initialfile']=self.entryErr.get()
		#Label the window
		browser_opt['title']='Select Error spectrum'
		#Run the dialog, and save the output to the error spectrum file name
		self.entryErr.set(tkFileDialog.askopenfilename(**browser_opt))
		return

	#TK widget for finding the output continuum. Save it to ENTRYOUT when found
	def getOutFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current spectrum file
		browser_opt['initialfile']=self.entryOut.get()
		#Label the window
		browser_opt['title']='Name of output continuum file'
		#Run the dialog, and save the output to the output continuum file name
		self.entryOut.set(tkFileDialog.askopenfilename(**browser_opt))
		return

	#TK widget for finding the output continuum-fitted spectrum. Save it to ENTRYCONT when found
	def getContFile(self):
		#Options for the browser dialog window
		browser_opt={}
		#Have access to all files
		browser_opt['filetypes']=[('all files','.*')]
		#Set initial file to current spectrum file
		browser_opt['initialfile']=self.entryCont.get()
		#Label the window
		browser_opt['title']='Name of output continuum-fitted spectrum'
		#Run the dialog, and save the output to the continuum-fitted spectrum file name
		self.entryCont.set(tkFileDialog.askopenfilename(**browser_opt))
		return

	#Function to toggle tutorial mode on/off.
	def onTutorial(self):
		#USETUTORIAL: A global variable to keep track of whether the tutorial mode is on (TRUE) or off (FALSE)
		global usetutorial
		#Turn off if previously on. Notify user
		if usetutorial:
			usetutorial=False
			tkMessageBox.showinfo("Help Message", "Tutorial mode is now off.")
		#Turn on if previously off. Notify user
		else:
			usetutorial=True
			tkMessageBox.showinfo("Help Message", "Tutorial mode is now on.")
		return
	#Function to start-up continuum tweaker, and save output to file.
	def onTweak(self):
		#Let user know what to do in tutorial mode.
		if usetutorial: tkMessageBox.showinfo("Help Message", "Use buttons in main widget to add/remove spline points. Once happy, save before exiting.")
		#Get the names of the input spectrum, error spectrum, and output files
		inspec=self.entryIn.get()
		outspec=self.entryOut.get()
		errspec=self.entryErr.get()
		fcont=self.entryCont.get()

		print "Out", outspec
		print "Cont", fcont
		print "ERR", errspec

		#Load the input and error spectrum. ?WVLNGTH and ?SPECTRUM are the wavelength and flux arrays for the two files
		iwvlngth, ispectrum=specfits(inspec)
		ewvlngth, espectrum=specfits(errspec)
		#Define the ouput spline python binary file by adding .p to the end of the output continuum file 
		splinefile=fcont+'.p'
		#Load the spline from the previously-generated splinefile.
		[xspline,yspline]=cPickle.load(open(splinefile,'rb'))
		#Run the continuum tweaker
		CT=ContTweaker(iwvlngth,ispectrum,espectrum,xspline,yspline,fcont)

		#Save final tweaked spectrum. This takes the most up-to-date pickled SPLINE file
		#And uses this to generate spectrum. This way if you DO NOT save in the Continuum fitter
		#Nothing new will happen.
		print "Updating output files..."
		#Load the spline from the pickle file
		[xspline,yspline]=cPickle.load(open(splinefile,'rb'))
		#Interpolate the continuum from the cubic spline.
		spline=scipy.interpolate.interp1d(xspline,yspline,kind='cubic')
		#Generate the continuum
		cont=spline(iwvlngth)
		#Define name of output continuum-fitted error spectrum
		outespec=outspec.split('.fit')[0]+'_err.fits'
		#Save to continuum-fitted spectrum, as well as continuum-fitted errorspecturm
		savefits(ispectrum/cont,inspec,outspec)
		if os.path.isfile(errspec): savefits(espectrum/cont,errspec,outespec)

		return
	#Function for starting up the continuum fitting routine
	def OnFitCont(self):
		#Tell user how to use continuum fitter 
		if usetutorial: tkMessageBox.showinfo("Help Message", "Click along the continuum. When finished, click left of the first red point to refresh window.")
		#Get file names/wavelength chunk range from fields in widget
		inspec=self.entryIn.get()
		outspec=self.entryOut.get()
		xr=float(self.entryXR.get())
		#Load input spectrum to be fitted
		wvlngth, spectrum=specfits(inspec)
		#Set the minimum x-value for the continuum fitting window to be the smallest wavelength - 5% of the chunk range
		xmin=min(wvlngth)-0.05*xr
		#Set the maximum x-value to be the size of the wavelength chunk + 5%. 
		xmax=xmin+1.10*xr
		#Calculate the number of iterations required to loop through the entire
		#spectrum based on the wavelength chunk size
		numiters=int(math.ceil((max(wvlngth)-xmin)/xr))+1
		#Set the first value of the spline to the first entry of the spectrum.
		#WARNING: Be careful here as it might be a bad pixel, but we need to start somewhere...
		xspline=[wvlngth[0]]
		yspline=[spectrum[0]]

		#Create matplotlib figure to display the full spectrum
		#Use this to help the user to see how the fit is going as well as seeing what is coming
		#up next in the spectrum for fitting.
		self.fullspecfig=plt.figure()#MPL figure for full spectrum
		self.FFroot=Tkinter.Tk()#TK widget window to embed full spectrum figure in
		self.FFroot.wm_title("Full spectrum")#Window title
		#Set up the MPL canvas to embed the figure in, which is embedded in the TK window
		self.FullSpecPlot=FigureCanvasTkAgg(self.fullspecfig,master=self.FFroot)
		self.FullSpecPlot.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
		#Set up MPL toolbar for zooming, panning, etc.
		NavFullSpecPlot=NavigationToolbar2TkAgg(self.FullSpecPlot, self.FFroot)
		#Define the axes for the plot, and plot the spectrum. Label appropriately
		self.fullax=self.fullspecfig.add_subplot(1,1,1)
		self.fullax.plot(wvlngth,spectrum,'k',drawstyle='steps')
		self.fullax.set_ylabel('Flux')
		self.fullax.set_xlabel('Wavelength')
		self.fullax.set_title(inspec)
		#Draw the spectrum in the window.
		self.FullSpecPlot.draw()

		#Create matplotlib figure window for the continuum-fitting figure. Displays spectrum chunk
		self.specfig=plt.figure()#MPL figure displaying the spectrum chunk
		self.SFroot=Tkinter.Tk()#Window for spectrum chunk display
		self.SFroot.wm_title("Continuum fitting display")
		#MPL Canvas for spectrum chunk within window.
		self.SpecPlot=FigureCanvasTkAgg(self.specfig,master=self.SFroot)
		self.SpecPlot.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
		NavSpecPlot=NavigationToolbar2TkAgg(self.SpecPlot, self.SFroot)
		#Create axes instance for plotting the spectrum chunk to be fitted.
		self.ax=self.specfig.add_subplot(1,1,1)
		#Loop through each wavelength chunk and ... 
		for ii in range(numiters):
			#OFF is an offest to add
			#off=0.0
			#if ii==0:
			#	off=0.95*xr
			#xmax=xmax-off
			#Clear the chunk axes isntance from previous run.
			self.ax.clear()
			#Do I need this???
			#self.SpecPlot.draw()
			#self.specfig.set_canvas(self.SpecPlot)
			#self.ax=self.specfig.add_subplot(1,1,1)

			#Find the index of WVLNGTH where it is equal to XMIN (or as cloase as possible)
			minind=np.where(wvlngth>=xmin)[0][0]
			#Find the index of WVLNGTH where it is equal to XMAX (or as close as possible)
			maxind=len(wvlngth)-1#Initally set to the largest wavelength possible
			#Want to add a little bit more than the chunk size to display to guide the user's eye...
			if xmax+0.05*xr<=max(wvlngth):
				maxind=np.where(wvlngth>=xmax+0.05*xr)[0][0]
			#If somehow the indicies are too small/large, set to the smallest/largest they can be
			if minind<0: minind=0
			if maxind>len(wvlngth)-1: maxind=len(wvlngth)-1
			#With the indicies defining the next spectrum chunk, get the spectrum to be plotted
			tmpwvl=wvlngth[minind:maxind]
			tmpspec=spectrum[minind:maxind]
			
			#Plot the spectrum
			self.ax.plot(tmpwvl,tmpspec,'k',drawstyle='steps')
			#Set the plot limits to +/-10% of the chunk range, so you can see all the data.
			self.ax.set_xlim(xmin-0.10*xr,xmax+0.10*xr)
			#To set the y-limits of the window, select the upper/lower 2.5 percentiles of the flux
			pymin=np.percentile(tmpspec,'2.5')
			pymax=np.percentile(tmpspec,'97.5')
			#Define the y-range as the difference in these percentiles
			yr=pymax-pymin
			#If the last spline point is not within this range, adjust the range appropriately
			if yspline[-1]>(pymax): pymax=yspline[-1]
			if yspline[-1]<(pymin): pymin=yspline[-1]
			#set the y-limits of the plot to the minimum/maximum values defined above.
			#Add 10% of the range to make sure that any emission/absorption lines aren't
			#missed in the display
			self.ax.set_ylim(pymin-0.10*yr, pymax+0.10*yr)
			#Plot the last point of the continuum fit spline for reference
			self.ax.plot(xspline[-1],yspline[-1], 'or')
			#Get the y-limits of the plot to plot vertical lines
			#indicating the user where to start/stop continuum fitting.
			ymin,ymax=self.ax.get_ylim()
			self.ax.plot([xmax,xmax],[ymin,ymax],':r')
			self.ax.plot([xspline[-1],xspline[-1]],[ymin,ymax],':r')
			#define the LINE artist that will display the spline of the continuum fit.
			line,=self.ax.plot(xspline[-1],yspline[-1])
			self.ax.set_title(inspec)
			#update the plot
			self.SpecPlot.draw()
			#Run the continuum fitter (LINEBUILDER CLASS)
			linebuilder= LineBuilder(line, xmin,xmax,self.SpecPlot,self.fullax,self.FullSpecPlot)
			#grab the modifed full-spectrum axes/Canvas from LINEBUILDER to keep for next chunk
			self.fullax=linebuilder.fullax
			self.FullSpecPlot=linebuilder.FullSpecPlot
			#Grab the x/y coordiantes of the continuum from the user clicking
			xcoos=linebuilder.xs
			ycoos=linebuilder.ys
			#Add the coordinates to the spline lists
			for ii in range(len(xcoos)):
				#ignore the first point because it was the last point in the fit from the previous chunk
				if ii>0:
					xspline.append(xcoos[ii])
					yspline.append(ycoos[ii])
			#reset the xmin/xmax to the last pline point -5%of the chunk size, and add a new chunk size
			xmin=xspline[-1]-0.05*xr
			xmax+=xr
			#Stop the plot if the extent is further than the spectrum
			if xmax>max(wvlngth): xmax=max(wvlngth)
		#once finished, close all the plots!
		plt.close('all')
		self.SFroot.destroy()
		self.FFroot.destroy()


		#This block of code does some quick tweaking to save a couple of issues.
		#First: The first spline point was set to be the first picel in the spectrum. This is wrong.
		#To correct, we will take second spline point, and just extend a flat line it back to the first wavelength.
		yspline[0]=yspline[1]
		#Second: The last spline point will be at the last wavelength of the spectrum. Let's set the
		#y-value to the last yspline value (again, a flat line)
		xspline.append(wvlngth[-1])
		yspline.append(yspline[-1])


		#Interpolate a cubic fit to the spline points, and generate the SPLINE instance 
		spline=scipy.interpolate.interp1d(xspline,yspline,kind='cubic')
		#Save the spline to the pickle file, and the continuum to an ascii file
		SaveSpline(xspline,yspline,wvlngth,self.entryCont.get())
		#Get the continuum from the SPLINE instance using the input wavelength
		cont=spline(wvlngth)
		#Get the continuum-fitted spectrum by taking the spectrum and dividing by the continuum
		cont_spec=spectrum/cont
		#Get the error spectrum name, create arrays for the error spectrum
		errspec=self.entryErr.get()
		ewvlngth=np.zeros(len(wvlngth))
		espectrum=np.zeros(len(wvlngth))
		cont_espec=np.zeros(len(wvlngth))
		#If the error spectrum exists, load into these arrays.
		if os.path.isfile(errspec):
			outespec=outspec.split('.fit')[0]+'_err.fits'
			ewvlngth, espectrum=specfits(errspec)
			econt=spline(ewvlngth)
			#Generate the continuum-fitted error spectrum
			cont_espec=espectrum/econt


		#The next step is to display what the resulting continuum fit looks like.
		#If using turorial, inform user of this.
		if usetutorial: tkMessageBox.showinfo("Help Message", "Verify continuum fit looks good. If not, redo fit or use Tweak Continuum feature.")

		# make a figure with two plots. The top panel will be the spectrum witht he continuum overlayed
		#The second panel will be the continuum-fitted spectrum
		plt.figure()
		plt.title('Spline fit')
		plt.subplot(2,1,1)
		plt.plot(wvlngth, spectrum,'k',drawstyle='steps', label='Spectrum')
		#Plot the spectrum and the spline points
		plt.plot(wvlngth,cont,'r',label='Fit')
		plt.plot(xspline,yspline,'ob',label='Spline')
		#Plot the error for reference
		plt.plot(ewvlngth,espectrum,'--g',label='Error')
		plt.ylabel('Flux')
		ymin=np.percentile(spectrum,0.5)
		ymax=np.percentile(spectrum,99.5)
		plt.ylim(ymin,ymax)
		plt.legend(ncol=4,frameon=False, loc=9, bbox_to_anchor=(0.5, 1.3))
		plt.subplot(2,1,2)
		#Plot the continuum-fitted spectrum.
		plt.plot(wvlngth,cont_spec,'r')
		#Here, you may expect to zeroth order the scatter around the continuum should
		#be the same size as the error. Therefore, plot the 1 +/- error spectrum (after cotninuum fit)
		#to demonstrate the expected scatter about the continuum 
		plt.plot(ewvlngth,1.0-cont_espec,'--g')	
		plt.plot(ewvlngth,1.0+cont_espec,'--g')
		xmin,xmax=plt.xlim()
		plt.plot([xmin,xmax],[1,1],'--k')
		plt.ylim(-0.5,1.5)
		plt.ylabel('Relative flux')
		plt.show()
		#Once the user closes the window, save the spectrum and error spectrum 
		savefits(cont_spec,inspec,outspec)		
		if os.path.isfile(errspec): savefits(cont_espec,errspec,outespec)
		print 'Saved to %s'%outspec

#Main program loop.
if __name__=="__main__":
        app=simpleapp_tk(None)
        app.title('TB_CONTFIT')#Name of application
        print "Starting TB's continuum fitter"
        app.mainloop()

