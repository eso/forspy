try:
  import numpy
  from pipeline_display import *
  from pipeline_product import *
  from reflex_plot_widgets import *
  from numpy.polynomial import Polynomial
  numpy.seterr(invalid='ignore')
except ImportError:
  donothing=1


class PlotableReducedArc :
  def __init__(self, fits):
    self.arc     = PipelineProduct(fits)
    self.arcdisp = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    self.arc.readImage()
    self.arc.read2DLinearWCS()

  def plot(self, subplot, title, tooltip):
    self.arcdisp.setLabels('Lambda [Angstrom]', 'Y [pix]')
    self.arcdisp.setZLimits((-100., 1000.))
    self.arcdisp.setXLinearWCSAxis(self.arc.crval1, self.arc.cdelt1, 
                                   self.arc.crpix1)
    self.arcdisp.display(subplot, title, tooltip, self.arc.image)
    
class PlotableFlat(object) :
  def __init__(self, fits_flat, fits_slittrace):
    self.flat      = PipelineProduct(fits_flat)
    self.slittrace = None
    if fits_slittrace is not None:
      self.slittrace = PipelineProduct(fits_slittrace)
    self.flatdisp  = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the flat image
    self.flat.readImage()
    
    #Reading the polinomial traces
    if self.slittrace is not None:
      ndegree = self.slittrace.getTableNcols(1) - 1
      self.nslits  = self.slittrace.getTableNrows(1) / 2
      degreecols = []
      for deg in range(ndegree):
        colname = 'c%d'%deg
        self.slittrace.readTableColumn(1, colname)
        degreecols.append(self.slittrace.column)
    
      top_trace_polynomials = []
      bottom_trace_polynomials = []
      for slit in range(self.nslits) :
        top_trace_coeff = []
        bottom_trace_coeff = []
        for deg in range(ndegree) :
          top_trace_coeff.append(degreecols[deg][2*slit])
          bottom_trace_coeff.append(degreecols[deg][2*slit + 1])
        
        top_trace_pol = Polynomial(top_trace_coeff)  
        bottom_trace_pol = Polynomial(bottom_trace_coeff)
        top_trace_polynomials.append(top_trace_pol) 
        bottom_trace_polynomials.append(bottom_trace_pol) 

      #Creating the points to plot based on the polynomail traces
      self.xpos_traces = []
      self.ypos_top_traces = []
      self.ypos_bottom_traces = []
      for slit in range(self.nslits) :
        ypos_top = []
        ypos_bottom = []
        xpos = []
        for xpix in range(self.flat.image.shape[1]) :
          xpos.append(xpix+1) 
          ypos_top.append(top_trace_polynomials[slit](xpix)+1) 
          ypos_bottom.append(bottom_trace_polynomials[slit](xpix)+1)
        self.xpos_traces.append(xpos)
        self.ypos_top_traces.append(ypos_top) 
        self.ypos_bottom_traces.append(ypos_bottom)

  def plot(self, subplot, title, tooltip):
    self.flatdisp.setLabels('X [pix]', 'Y [pix]')
    self.flatdisp.display(subplot, title, tooltip, self.flat.image)
    
    if self.slittrace is not None:
      subplot.autoscale(enable=False)
      for slit in range(self.nslits) :
        subplot.plot(self.xpos_traces[slit], self.ypos_top_traces[slit],
                     linestyle='solid',color='red')
        subplot.plot(self.xpos_traces[slit], self.ypos_bottom_traces[slit],
                     linestyle='solid',color='darkred')

class PlotableNormFlat (PlotableFlat) :
  def __init__(self, fits_flat, fits_slittrace):
    super(PlotableNormFlat, self).__init__(fits_flat, fits_slittrace)

  def plot(self, subplot, title, tooltip):
    self.flatdisp.setZLimits((0.9, 1.1))
    self.flat.image[self.flat.image > 5.] = 0
    super(PlotableNormFlat, self).plot(subplot, title, tooltip)

class PlotableRawFlat (PlotableFlat) :
  def __init__(self, fits_flat_raw, fits_master_flat, fits_slittrace):
    super(PlotableRawFlat, self).__init__(fits_flat_raw, fits_slittrace)
    if self.slittrace is not None and fits_master_flat is not None :
      master_flat = PipelineProduct(fits_master_flat)
      self.trimm_lly = master_flat.all_hdu[0].header.get('HIERARCH ESO QC TRIMM LLY')

      #Change the traces by the amount of overscan in Y that has been removed
      for ypos_top in self.ypos_top_traces:
        for j, ypos in enumerate(ypos_top):
          ypos_top[j] = ypos + self.trimm_lly - 1
      for ypos_bottom in self.ypos_bottom_traces:
        for j, ypos in enumerate(ypos_bottom):
          ypos_bottom[j] = ypos + self.trimm_lly -1
    

  def plot(self, subplot, title, tooltip):
    super(PlotableRawFlat, self).plot(subplot, title, tooltip)

class PlotableSpatialMap :
  def __init__(self, fits_spatialmap):
    self.spatialmap      = PipelineProduct(fits_spatialmap)
    self.spatialmapdisp  = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the flat image
    self.spatialmap.readImage()

  def plot(self, subplot, title, tooltip):
    self.spatialmapdisp.setLabels('X', 'Y')
    self.spatialmapdisp.setZLimits((0., 100))
    self.spatialmapdisp.display(subplot, title, tooltip, self.spatialmap.image)

class PlotableMappedScience :
  def __init__(self, fits_mappedscience, fits_objecttable):
    self.mappedscience      = PipelineProduct(fits_mappedscience)
    self.mappedsciencedisp  = ImageDisplay()
    if fits_objecttable is not None:
      self.objecttable        = PipelineProduct(fits_objecttable)
    else :
      self.objecttable        = None
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the flat image
    self.mappedscience.readImage()
    
    #Reading the object table
    if self.objecttable is not None:
      nslit = self.objecttable.getTableNrows(1)
      maxobjectperslit = (self.objecttable.getTableNcols(1) - 7 ) / 4
      start_extracted_cols = []
      end_extracted_cols = []
      for obj in range(maxobjectperslit):
        colname = 'start_%d'%(obj+1)
        self.objecttable.readTableColumn(1, colname)
        start_extracted_cols.append(self.objecttable.column)
        colname = 'end_%d'%(obj+1)
        self.objecttable.readTableColumn(1, colname)
        end_extracted_cols.append(self.objecttable.column)

      self.ybottom_obj_extract = []
      self.ytop_obj_extract = []
      for slit in range(nslit) :
        for obj in range(maxobjectperslit) :
          ybottom = start_extracted_cols[obj][slit]
          ytop    = end_extracted_cols[obj][slit]
          if ybottom != -1 :
            self.ybottom_obj_extract.append(ybottom)
            self.ytop_obj_extract.append(ytop)
      self.nobjects = len(self.ybottom_obj_extract)

  def plot(self, subplot, title, tooltip):
    self.mappedsciencedisp.setLabels('X [pix]', 'Y [pix]')
    self.mappedsciencedisp.setZLimits((0., 0.9))
    self.mappedsciencedisp.display(subplot, title, tooltip, self.mappedscience.image)
    if self.objecttable is not None:
      subplot.autoscale(enable=False)
      for obj in range(self.nobjects) :
        subplot.axhline(self.ytop_obj_extract[obj], linestyle='solid',color='red')
        subplot.axhline(self.ybottom_obj_extract[obj], linestyle='solid',color='yellow')
        
  def getObjectInPosition(self, ypos) :
    for obj in range(self.nobjects) :
      if ypos > self.ybottom_obj_extract[obj] and \
         ypos < self.ytop_obj_extract[obj] :
        return self.nobjects - obj
    return -1

class PlotableDispResiduals :
  def __init__(self, fits_dispresiduals):
    self.dispresiduals = PipelineProduct(fits_dispresiduals)
    self.resdisplay  = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the residuals table
    self.dispresiduals.readTableColumn(1, 'wavelength')
    self.wave = self.dispresiduals.column
    nwave    = self.dispresiduals.getTableNrows(1)
    ncolumns = self.dispresiduals.getTableNcols(1)
    nselectedrows = (ncolumns - 1) // 3
    self.residuals = []
    self.allwave = []
    self.allypos = []
    self.allresiduals = []
    for i in range(nselectedrows) :
      #TODO: Currently the residuals are computed every 10 rows. 
      #This is hard-coded in the pipeline. It would be better just to detect the
      #columns whose name start with 'r' 
      colname = 'r%d'%(i*10) 
      self.dispresiduals.readTableColumn(1, colname)
      row_residuals = self.dispresiduals.column
      self.residuals.append(row_residuals)
      self.allwave.extend(self.wave)
      self.allresiduals.extend(row_residuals)
      ypos = i*10.
      self.allypos.extend([ypos] * nwave)

  def plotResVsWave(self, subplot, title, tooltip):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [pix]')
    self.resdisplay.display(subplot, title, tooltip, self.allwave,
                            self.allresiduals)

  def plotResVsY(self, subplot, title, tooltip):
    self.resdisplay.setLabels('Ypos [pix]','Residual [pix]')
    self.resdisplay.display(subplot, title, tooltip, self.allypos,
                            self.allresiduals)
  def getClosestLine(self, wave_selected) :
    
    distance = numpy.fabs(self.wave - wave_selected)
    idx = numpy.nanargmin(distance)
    return self.wave[idx]

class PlotableDetectedLines :
  def __init__(self, fits_detectedlines):
    self.detectedlines = PipelineProduct(fits_detectedlines)
    self.xydisplay     = ScatterDisplay()
    self.resdisplay    = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the residuals table
    try :
      self.detectedlines.readTableColumn(1, 'xpos_rectified')
      self.x_pix = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'ypos_rectified')
      self.y_pix = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'xpos_rectified_iter')
      self.x_pix_iter = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'ypos_rectified_iter')
      self.y_pix_iter = self.detectedlines.column
    except KeyError:
      self.detectedlines.readTableColumn(1, 'xpos')
      self.x_pix = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'ypos')
      self.y_pix = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'xpos_iter')
      self.x_pix_iter = self.detectedlines.column
      self.detectedlines.readTableColumn(1, 'ypos_iter')
      self.y_pix_iter = self.detectedlines.column

      
    self.detectedlines.readTableColumn(1, 'wave_ident')
    self.wave  = self.detectedlines.column
    self.detectedlines.readTableColumn(1, 'wave_ident_iter')
    self.wave_iter  = self.detectedlines.column
    self.detectedlines.readTableColumn(1, 'res_xpos')
    self.res_xpos  = self.detectedlines.column

  def plotXVsY(self, subplot, title, tooltip):
    #We first plot all the detected lines
    self.xydisplay.setLabels('Xpos [pix]','Ypos [pix]')
    self.xydisplay.setColor('black')
    self.xydisplay.display(subplot, title, tooltip, self.x_pix,
                           self.y_pix)
    #We then overplot the identified lines in the second iteration
    self.xydisplay.setColor('lightgreen')
    self.xydisplay.display(subplot, title, tooltip,
                           self.x_pix_iter[numpy.isfinite(self.wave_iter)],
                           self.y_pix_iter[numpy.isfinite(self.wave_iter)])
    #And then we overplot the identified lines in the first iteration
    self.xydisplay.setColor('green')
    self.xydisplay.display(subplot, title, tooltip, 
                           self.x_pix[numpy.isfinite(self.wave)],
                           self.y_pix[numpy.isfinite(self.wave)])

  def plotResVsWave(self, subplot, title, tooltip, excluded_lines = None):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [pix]')
    self.resdisplay.setColor('black')
    self.resdisplay.display(subplot, title, tooltip, 
                            self.wave[numpy.isfinite(self.res_xpos)],
                            self.res_xpos[numpy.isfinite(self.res_xpos)])
    if excluded_lines is not None :
      for line in excluded_lines : 
        subplot.axvline(line, linestyle='solid',color='red')


class PlotableSkylinesOffsets :
  def __init__(self, fits_skylines_off):
    self.skylines_off  = PipelineProduct(fits_skylines_off)
    self.resdisplay    = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the slylines offset table
    nslits = self.skylines_off.getTableNcols(1) - 1

    skylines_wave = self.skylines_off.readTableColumn(1, 'wave')
    self.allskylines_wave = list()
    self.allwave_res = list()
    
    for col in range(nslits) :
      self.allskylines_wave.extend(skylines_wave)
      wave_res = self.skylines_off.readTableColumn(1, col + 1)
      self.allwave_res.extend(wave_res)
      
  def plot(self, subplot, title, tooltip):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [Ang]')
    self.resdisplay.setColor('black')
    self.resdisplay.setPointSize(7)
    self.resdisplay.display(subplot, title, tooltip, 
                            self.allskylines_wave, self.allwave_res)


class PlotableExtractedScience :
  def __init__(self, fits_extractedscience):
    self.obj_id = -1
    self.extractedscience  = PipelineProduct(fits_extractedscience)
    self.spectrumdisplay   = SpectrumDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the flat image
    self.extractedscience.readImage()
    self.nobj = self.extractedscience.image.shape[0]
    self.crpix1  = self.extractedscience.readKeyword('CRPIX1', 0)
    self.crval1  = self.extractedscience.readKeyword('CRVAL1', 0)
    self.cdelt1  = self.extractedscience.readKeyword('CD1_1', 0)
    self.bunit = self.extractedscience.readKeyword('BUNIT', 0)
    self.nwave  =   self.extractedscience.image.shape[1]
    self.wave = numpy.arange(1, self.nwave+1, 1)
    self.wave = (self.wave - self.crpix1) * self.cdelt1 + self.crval1
    if(self.obj_id == -1) : # Select brightest
      self.selectBrightest()
    self.setFluxSelected()

  def selectBrightest(self):
    if self.nobj == 1:
      self.obj_id = 1
    median = 0
    for obj in range(self.nobj) :
      new_median = numpy.median(self.extractedscience.image[obj,:]) 
      if new_median > median :
        median = new_median
        self.obj_id = obj + 1
  
  def setFluxSelected(self) :
    self.flux = self.extractedscience.image[self.obj_id-1,:]

  def selectObject(self, obj_id):
    self.obj_id = obj_id
    self.setFluxSelected()

  def plot(self, subplot, title, tooltip):

    self.spectrumdisplay.setLabels('Lambda', 'Total Flux ['+self.bunit+']')
    self.spectrumdisplay.display(subplot, title, tooltip, self.wave, self.flux,
                                autolimits = True)

class PlotableSpecPhot :
  def __init__(self, fits):
    self.resp     = PipelineProduct(fits)
    self.respdisp = SpectrumDisplay()
    self.tabdisp  = ScatterDisplay()
    self.flat_sed = False
    self.loadFromFits()

  def loadFromFits(self) :
    self.wave         = self.resp.readTableColumn(1, 'WAVE')
    self.wave_obs     = self.resp.readTableColumn(2, 'WAVE')
    self.std_ref_flux = self.resp.readTableColumn(1, 'STD_FLUX')
    self.std_obs_flux = self.resp.readTableColumn(1, 'OBS_FLUX')
    if 'RESPONSE' in self.resp.all_hdu[2].columns.names :
      self.fit_response = self.resp.readTableColumn(2, 'RESPONSE')
      self.raw_response = self.resp.readTableColumn(1, 'RAW_RESPONSE')
    else :
      self.fit_response = self.resp.readTableColumn(2, 'RESPONSE_FFSED')
      self.raw_response = self.resp.readTableColumn(1, 'RAW_RESPONSE_FFSED')
      self.flat_sed = True
    self.used_fit     = self.resp.readTableColumn(1, 'USED_FIT')
    self.raw_response_nonnull = self.raw_response[self.raw_response > 0]
    self.wave_nonnull = self.wave[self.raw_response > 0]
    self.wave_used = self.wave[self.used_fit > 0]
    self.raw_response_used = self.raw_response[self.used_fit > 0]    

  def plotResponse(self, subplot, title, tooltip):
    self.respdisp.setLabels('Angstrom','10^ (-16) erg/(cm^ (-2) e-)')
    self.respdisp.flux_lim = 0., numpy.max(self.raw_response_nonnull) * 1.1
    self.respdisp.display(subplot, title, tooltip, self.wave_obs, self.fit_response, autolimits = False)
    subplot.scatter(self.wave_nonnull, self.raw_response_nonnull, color='darkblue')
    subplot.scatter(self.wave_used, self.raw_response_used, color='lightgreen')

  def plotStdExtracted(self, subplot, title, tooltip):
    self.respdisp.setLabels('Angstrom','e-/ (s Angstrom)')
    std_obs_flux_nonnull = self.std_obs_flux[self.std_obs_flux > 0]
    wave_nonnull = self.wave[self.std_obs_flux > 0]
    self.respdisp.display(subplot, title, tooltip,
                          wave_nonnull, std_obs_flux_nonnull, autolimits = True)

  def plotStdTabulated(self, subplot, title, tooltip):
    self.tabdisp.setLabels('Angstrom','10^ (-16) erg/(cm^ (-2) e-)')
    self.tabdisp.display(subplot, title, tooltip, self.wave, self.std_ref_flux)

class PlotableStdTabRedFlux :
  def __init__(self, reducedfluxstd_fits, reducedstd_fits, specphot_fits):
    self.reducedfluxstd = PlotableExtractedScience(reducedfluxstd_fits)
    self.reducedstd     = PlotableExtractedScience(reducedstd_fits)
    self.specphot       = PlotableSpecPhot(specphot_fits)
    self.tabstddisp     = ScatterDisplay()
    self.stdreddisp     = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #This will select the brightest spectrum, which is the criteria
    #used to extract the standard star
    self.reducedfluxstd.loadFromFits()
    self.reducedstd.loadFromFits()
    self.specphot.loadFromFits()
    self.std_ref_flux_nonnull = self.specphot.std_ref_flux[self.specphot.raw_response > 0]
    self.std_ref_flux_used = self.specphot.std_ref_flux[self.specphot.used_fit > 0]
    

  def plotStdTabRedFlux(self, subplot, title, tooltip) :
    self.tabstddisp.setLabels('Angstrom','10^ (-16) erg/(cm^ (-2) Angstrom)')
    self.tabstddisp.setLimits(self.reducedfluxstd.wave[0],
                              self.reducedfluxstd.wave[len(self.reducedfluxstd.wave)-1],
                              0.,
                              numpy.max(self.reducedfluxstd.flux) * 1.1)
    self.tabstddisp.setColor('red') 
    self.tabstddisp.display(subplot, title, tooltip, 
                            self.reducedfluxstd.wave, self.reducedfluxstd.flux)
    self.tabstddisp.setColor('darkblue') 
    self.tabstddisp.setPointSize(20) 
    subplot.scatter(self.specphot.wave, self.specphot.std_ref_flux)
    subplot.scatter(self.specphot.wave_used, self.std_ref_flux_used, color='lightgreen')


  def plotStdRed(self, subplot, title, tooltip) :
    self.stdreddisp.setLabels('Angstrom','10^ (-16) erg/(cm^ (-2) Angstrom)')
    self.stdreddisp.setLimits(self.reducedstd.wave[0],
                              self.reducedstd.wave[len(self.reducedstd.wave)-1],
                              0.,
                              numpy.max(self.reducedstd.flux) * 1.1) 
    self.stdreddisp.setColor('red') 
    self.stdreddisp.display(subplot, title, tooltip, 
                            self.reducedstd.wave, self.reducedstd.flux)



