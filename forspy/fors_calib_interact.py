# import the needed modules
try:
  import reflex
  import_sucess = True

#NOTE for developers: 
# -If you want to modify the current script to cope
#  with different parameters, this is the function to modify:
#  setInteractiveParameters()
# -If you want to modify the current script to read different data from
#  the input FITS, this is the function to modify:
#  readFitsData()                  (from class DataPlotterManager) 
# -If you want to modify the current script to modify the plots (using the same
#  data),  this is the function to modify:
#  plotProductsGraphics()          (from class DataPlotterManager)
# -If you want to modify the text that appears in the "Help" button,
#  this is the function to modify:
#  setWindowHelp()
# -If you want to modify the title of the window, modify this function:
#  setWindowTitle()


  #This class deals with the specific details of data reading and final plotting.
  class DataPlotterManager:
    # This function will read all the columns, images and whatever is needed
    # from the products. The variables , self.plot_x, self.plot_y, etc...
    # are used later in function plotProductsGraphics().
    # Add/delete these variables as you need (only that plotProductsGraphics()
    # has to use the same names).
    # You can also create some additional variables (like statistics) after
    # reading the files.
    # If you use a control variable (self.xxx_found), you can modify 
    # later on the layout of the plotting window based on the presence of 
    # given input files. 
    # sof contains all the set of frames
    def readFitsData(self, fitsFiles):
      #Initialise the objects to read/display
      self.lamp_reduced    = None
      self.slit_map        = None
      self.flat_norm       = None
      self.flat_mscreen    = None
      self.disp_residuals  = None
      self.detected_lines  = None
      self.excluded_lines  = list()

      #Read all the products
      frames = dict()
      for frame in fitsFiles:
        if frame == '' :
          continue
        category = frame.category
        frames[category] = frame


      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'REDUCED_LAMP_'+inst_mode in frames:
          self.lamp_reduced = PlotableReducedArc(frames["REDUCED_LAMP_"+inst_mode])

      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'SPATIAL_MAP_'+inst_mode in frames:
          self.slit_map = PlotableSpatialMap(frames["SPATIAL_MAP_"+inst_mode])

      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'DISP_RESIDUALS_TABLE_'+inst_mode in frames:
          self.disp_residuals = PlotableDispResiduals(frames["DISP_RESIDUALS_TABLE_"+inst_mode])

      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'DETECTED_LINES_'+inst_mode in frames:
          self.detected_lines = PlotableDetectedLines(frames["DETECTED_LINES_"+inst_mode])

      curv_frame = None
      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'CURV_COEFF_'+inst_mode in frames :
          curv_frame = frames['CURV_COEFF_'+inst_mode]
                  
      flat_norm = None
      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'MASTER_NORM_FLAT_'+inst_mode in frames :
          flat_norm = frames['MASTER_NORM_FLAT_'+inst_mode]
          self.flat_norm  = PlotableNormFlat(flat_norm, curv_frame)
                  
      for inst_mode in ('MXU', 'MOS', 'LSS', 'LONG_MOS', 'LONG_MXU') :
        if 'SCREEN_FLAT_'+inst_mode in frames :
          self.flat_mscreen  = PlotableRawFlat(frames['SCREEN_FLAT_'+inst_mode],
                                               flat_norm, curv_frame)
                  
    # This function creates all the subplots. It is responsible for the plotting 
    # layouts. 
    # There can different layouts, depending on the availability of data
    # Note that subplot(I,J,K) means the Kth plot in a IxJ grid 
    # Note also that the last one is actually a box with text, no graphs.
    def addSubplots(self, figure):
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None \
         and self.detected_lines is not None:
        if self.slit_map is not None:
          self.subplot_slit_map       = figure.add_subplot(3,2,2)
          self.subplot_lamp_reduced   = figure.add_subplot(3,2,1)
        else:
          self.subplot_lamp_reduced   = figure.add_subplot(3,1,1)
        self.subplot_line_x_y       = figure.add_subplot(3,2,3)
        self.subplot_flat_raw       = figure.add_subplot(3,2,4)
        self.subplot_line_res_wave  = figure.add_subplot(3,2,5)
        self.subplot_flat_norm      = figure.add_subplot(3,2,6)
      else : 
        self.subtext_nodata      = figure.add_subplot(1,1,1)
          
    # This is the function that makes the plots.
    # Add new plots or delete them using the given scheme.
    # The data has been already stored in self.plot_x, self.plot_xdif, etc ...
    # It is mandatory to add a tooltip variable to each subplot.
    # One might be tempted to merge addSubplots() and plotProductsGraphics().
    # There is a reason not to do it: addSubplots() is called only once at
    # startup, while plotProductsGraphics() is called always there is a resize.
    def plotProductsGraphics(self):
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None \
         and self.detected_lines is not None :

        #Reduced lamp
        self.plotReducedLamp()

        #Spatial map
        self.plotSpatialMap()

        #Master screen flat 
        self.plotScreenFlat()

        #Master flat normalise
        self.plotNormalisedFlat()
        
        #Dispersion residuals vs wave
        self.plotResiduals()

        #Line positions X vs Y
        self.plotLinePositions()
        
      else :
        #Data not found info
        self.showNoData()
 
    def plotReducedLamp(self) :
      title_lamp_reduced   = 'Wavelength-calibrated arc lamp frame'
      tooltip_lamp_reduced ="""Wavelength-calibrated arc lamp frame.
Arc lines should be vertical without scatter or tilt. 
There should be no regions with no arc lines at all. 
Arc line coverage may vary with slit position."""
      self.lamp_reduced.plot(self.subplot_lamp_reduced, title_lamp_reduced,
                             tooltip_lamp_reduced)

    def plotSpatialMap(self) :
      if self.slit_map is not None:
        title_slit_map   = 'Slit spatial map'
        tooltip_slit_map ="""Reduced arc lamp frame."""
        self.slit_map.plot(self.subplot_slit_map, title_slit_map,
                           tooltip_slit_map)

    def plotScreenFlat(self) :
      title_flat_mscreen   = 'First raw flat field frame'
      tooltip_flat_mscreen ="""First raw flat field frame.
This frame is the first raw flat frame. It serves to verify that for MOS/MXU data all slits have been detected (compare with spatial map displayed above)."""
      self.flat_mscreen.plot(self.subplot_flat_raw, title_flat_mscreen,
                         tooltip_flat_mscreen)

    def plotNormalisedFlat(self) :
      title_flat_norm   = 'Normalised master flat frame'
      tooltip_flat_norm ="""Normalised master flat frame.
This is the result of the flat field normalisation. 
For LSS data it may make sense to keep the relative
flux variation along the slit to correct for slit illumination. 
For MOS/MXU data this generally does not work well."""
      self.flat_norm.plot(self.subplot_flat_norm, title_flat_norm,
                         tooltip_flat_norm)

    def plotResiduals(self) :
      title_line_res_wave   = 'Residuals of wavelength calibration'
      tooltip_line_res_wave ="""Residuals of wavelength calibration.
The residuals should not show any trends.  
Outliers may be tolerated, but at the edges 
they may indicate a questionable result. 
Clicking with middle button will add that line 
to the --ignore_lines recipe parameter."""
      self.subplot_line_res_wave.clear()
      self.detected_lines.plotResVsWave(self.subplot_line_res_wave,
                                        title_line_res_wave,
                                        tooltip_line_res_wave,
                                        self.excluded_lines)

    def plotLinePositions(self) :
      title_line_x_y   = 'Detected / Identified arc lines'
      tooltip_line_x_y ="""Detected/identified arc lines.
Black points are detected lines. 
Green points are identified lines.
Light green points are identified lines after pattern-matching iteration
(only if wradius >0).
Most detected lines should also be identified. 
If lines at the edges are not identified this means that
the dispersion relation there will be extrapolated."""
      self.detected_lines.plotXVsY(self.subplot_line_x_y,
                                   title_line_x_y, tooltip_line_x_y)

    def showNoData(self) :
      self.subtext_nodata.set_axis_off()
      self.text_nodata = 'Calibrations not found in the products'
      self.subtext_nodata.text(0.1, 0.6, self.text_nodata, color='#11557c', fontsize=18,
                               ha='left', va='center', alpha=1.0)
      self.subtext_nodata.tooltip='Calibrations  not found in the products'

  
    def plotWidgets(self) :
      widgets = list()
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None \
         and self.detected_lines is not None :
        # Clickable subplot
        self.clickableresidual = InteractiveClickableSubplot(
          self.subplot_line_res_wave, self.setExcludedLine)
        widgets.append(self.clickableresidual)
      return widgets

    def setExcludedLine(self, point) :
      excluded_line = self.disp_residuals.getClosestLine(point.xdata)
      if excluded_line in self.excluded_lines :
        self.excluded_lines.remove(excluded_line)
      else :
        self.excluded_lines.append(excluded_line)
      self.plotResiduals()

      param_value = str()
      for line in self.excluded_lines :
        param_value = param_value + '%.4f, ' % line
      if len(param_value) > 1 :
        param_value = param_value[:-2]
      new_params = list()
      new_params.append(reflex.RecipeParameter('fors_calib','ignore_lines',
                                               value=param_value))
      return new_params
                
    # This function specifies which are the parameters that should be presented
    # in the window to be edited.
    # Note that the parameter has to be also in the in_sop port (otherwise it 
    # won't appear in the window) 
    # The descriptions are used to show a tooltip. They should match one to one
    # with the parameter list 
    # Note also that parameters have to be prefixed by the 'recipe name:'
    def setInteractiveParameters(self):
      paramList = list()
      paramList.append(reflex.RecipeParameter('fors_calib','sradius',group='flat field norm',description='Smooth box radius for flat field along spatial direction (used if s_knots < 0)'))
      paramList.append(reflex.RecipeParameter('fors_calib','dradius',group='flat field norm',description='Smooth box radius for flat field along dispersion direction (if d_knots < 0)'))
      paramList.append(reflex.RecipeParameter('fors_calib','s_degree',group='flat field norm',description='Polynomial degree for the flat field fitting along spatial direction'))
      paramList.append(reflex.RecipeParameter('fors_calib','d_nknots',group='flat field norm',description='Number of knots in flat field fitting splines along dispersion direction'))
      paramList.append(reflex.RecipeParameter('fors_calib','fit_threshold',group='flat field norm',description='Threshold percentage for flat spline fitting with respect to the maximum'))
      paramList.append(reflex.RecipeParameter('fors_calib','stack_method',group='master flat',description='Frames combination method. <sum | mean | median | ksigma>'))
      paramList.append(reflex.RecipeParameter('fors_calib','kiter',group='master flat',description='Max number of iterations in ksigma method'))
      paramList.append(reflex.RecipeParameter('fors_calib','klow',group='master flat',description='Low threshold in ksigma method'))
      paramList.append(reflex.RecipeParameter('fors_calib','khigh',group='master flat',description='High threshold in ksigma method'))
      paramList.append(reflex.RecipeParameter('fors_calib','slit_ident',group='master flat',description='Attempt slit identification for MOS or MXU'))
      paramList.append(reflex.RecipeParameter('fors_calib','startwavelength',group='wave calib / distorsions',description='Start wavelength in spectral extraction'))
      paramList.append(reflex.RecipeParameter('fors_calib','endwavelength',group='wave calib / distorsions',description='End wavelength in spectral extraction'))
      paramList.append(reflex.RecipeParameter('fors_calib','wdegree',group='wave calib / distorsions',description='Degree of wave calib polynomial'))
      paramList.append(reflex.RecipeParameter('fors_calib','wradius',group='wave calib / distorsions',description='Search radius if iterating pattern-matching with first-guess method'))
      paramList.append(reflex.RecipeParameter('fors_calib','wreject',group='wave calib / distorsions',description='Rejection threshold in dispersion relation fit (pixel)'))
      paramList.append(reflex.RecipeParameter('fors_calib','wmode',group='wave calib / distorsions',description='Interpolation mode of wavelength solution applicable to LSS-like data (0 = no interpolation, 1 = fill gaps, 2 = global model)'))
      paramList.append(reflex.RecipeParameter('fors_calib','wmosmode',group='wave calib / distorsions',description='Interpolation mode of wavelength solution (0 = no interpolation, 1 = local (slit) solution, 2 = global model)'))
      paramList.append(reflex.RecipeParameter('fors_calib','dispersion',group='wave calib / distorsions',description='Expected spectral dispersion (Angstrom/pixel)'))
      paramList.append(reflex.RecipeParameter('fors_calib','peakdetection',group='wave calib / distorsions',description='Initial peak detection threshold (ADU)'))
      paramList.append(reflex.RecipeParameter('fors_calib','ignore_lines',group='wave calib / distorsions',description='Catalog lines nearest to wavelengths in this list will be ignored for wavelength calibration')) 
      paramList.append(reflex.RecipeParameter('fors_calib','used_linesets',group='wave calib / distorsions',description='Linesets to use. Valid are standard and extended (see column LINE_SET in the line catalogue)')) 
      return paramList

    def setWindowHelp(self):
      help_text = """
In this window, the user will interact with the Fors calibration (flat, slit traces and wavelength calibration)"""
      return help_text

    def setWindowTitle(self):
      title = 'Fors Interactive Calibration'
      return title

except ImportError:
  import_sucess = 'false'
  print "Error importing modules pyfits, wx, matplotlib, numpy"

#This is the 'main' function
def main(args=None):

  # import reflex modules
  from reflexy.base import reflex_interactive_app

  import sys

  # import UVES reflex modules
  from fors_plot_common import *

  # Create interactive application
  interactive_app = reflex_interactive_app.PipelineInteractiveApp()

  #Check if import failed or not
  if import_sucess == 'false' :
    interactive_app.setEnableGUI('false')

  #Open the interactive window if enabled
  if interactive_app.isGUIEnabled() :
    #Get the specific functions for this window
    dataPlotManager = DataPlotterManager()
    interactive_app.setPlotManager(dataPlotManager)
    interactive_app.showGUI()
  else :
    interactive_app.passProductsThrough()

  # print outputs
  interactive_app.print_outputs()

  sys.exit()
