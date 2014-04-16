# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import json
from pyworkflow.em import *
from views_util import *
from views_protocol import updateProtocolParams
from pyworkflow.manager import Manager
from pyworkflow.em import emProtocolsDict
from django.http import HttpResponse

from pyworkflow.gui.plotter import Plotter

from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
from pyworkflow.viewer import WEB_DJANGO

# XMIPP
from viewers.xmipp_ml2d import *
from viewers.xmipp_cl2d import *
#from viewers.xmipp_ml3d import *
from viewers.xmipp_nma import *
from viewers.xmipp_nma_align import *
# SPIDER
from viewers.spider_capca import *
from viewers.spider_ward import *
# BRANDEIS
from viewers.brandeis_frealign import *
# RELION
from viewers.relion import *


############## 1ST STEP: LAUNCH VIEWER METHODS ##############

def launch_viewer(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        project = loadProject(projectName)
        protId = request.GET.get('protocolId', None)
        protocol = project.mapper.selectById(int(protId))
        
        viewers = findViewers(protocol.getClassName(), WEB_DJANGO)
        
        if len(viewers) == 0:
            urls = []
            for _, output in protocol.iterOutputAttributes(EMObject):
                urls.append("/visualize_object/?objectId=%s" % output.getObjId())
                
            #msg = "There is not viewer for protocol: <strong>" + protocol.getClassName() +"</strong>"
            ioDict = {'urls': urls}
        else:
            viewer = viewers[0]()
            # Lets assume if there is a viewer for the protocol
            # it will be a ProtocolViewer with a Form
            if isinstance(viewer, ProtocolViewer):
                ioDict = viewerForm(project, protocol, viewer)
            else:
                views = viewer.visualize(protocol)
                urls =  [viewToUrl(request, v) for v in views]
                ioDict = {'urls': urls}
                
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
        
    return HttpResponse(jsonStr, mimetype='application/javascript')


def viewerForm(project, protocol, viewer):
    protId = protocol.getObjId()
    viewerClassName = viewer.getClassName()
    
    return {"url_form": "/form/?protocolClass="+ viewerClassName +
              "&protRunIdViewer="+ str(protId) +
              "&action=visualize"}
    
############## 2ND STEP: VIEWER FUNCTION METHODS ##############

def viewToUrl(request, view):
    
    if isinstance(view, Plotter):
        url = savePlot(request, view)
    if isinstance(view, DataView):
        url = "/visualize_object/?path="+ view.getPath()
    
    return url

def viewerElement(request):
    project, protocolViewer = loadProtocolProject(request)
    protId = request.POST.get('protRunIdViewer', None)
    viewerParam = request.POST.get('viewerParam', None)
    
    protocol = project.mapper.selectById(int(protId))
    protocolViewer.setProtocol(protocol)
    updateProtocolParams(request, protocolViewer, project)
    
    views = protocolViewer._getVisualizeDict()[viewerParam]()
    
    urls = [viewToUrl(request, v) for v in views]
    
    jsonStr = json.dumps({'urls':urls}, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')

############## AUX METHODS ##############
    
def view_plot_xmipp(request):
    projectName = request.session['projectName']
    project = loadProject(projectName)
    protId = request.GET.get('protocolId', None)
    protocol = project.mapper.selectById(int(protId))
    
    mdFn = getattr(protocol.outputParticles, '_xmippMd', None)
    if mdFn:
        fn = mdFn.get()
    else:
        fn = project.getTmpPath(protocol.outputParticles.getName() + '_images.xmd')
        writeSetOfParticles(protocol.outputParticles, fn, protocol.getTmpPath())
        
    md = xmipp.MetaData(fn)
    if md.containsLabel(xmipp.MDL_ZSCORE):
               
#        print "MD contains ZSCORE"
        xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
        xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
        xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
        figFn = fn.replace('.xmd', '.png')
        canvas = xplotter.getCanvas()
        
        response= HttpResponse(content_type='image/png')
        canvas.print_png(response)
        return response
    
