#
import os
import vtk
from vtk.util import numpy_support
import numpy as np
#
def para2d(nelx,nely,x,u,g,s):
#
    n_nds=(nelx+1)*(nely+1)
    n_elm=nelx*nely
#
    ply=vtk.vtkPolyData()
    pts=vtk.vtkPoints()
    cel=vtk.vtkCellArray()
#
    pts.SetNumberOfPoints(n_nds)
#
    scl=vtk.vtkFloatArray()
    scl.SetName('Decision')
    scl.SetNumberOfComponents(1)
    scl.SetNumberOfTuples(n_elm)
    scl1=vtk.vtkFloatArray()
    scl1.SetName('Sensitivity')
    scl1.SetNumberOfComponents(1)
    scl1.SetNumberOfTuples(n_elm)
#
    vec=vtk.vtkFloatArray()
    vec.SetName('Displacement')
    vec.SetNumberOfComponents(3)
    vec.SetNumberOfTuples(n_nds)
#
    c=0
    minu=0e0
    for i in range(nelx+1):
        for j in range(nely+1):
            pts.SetPoint(c,i,j,0)
            vec.SetTuple3(c,u[2*c],-u[2*c+1],0.)
            minu=min(minu, -u[2*c+1])
            c=c+1
#
    c=0;c1=0; c2=1; c3=nely+2; c4=nely+1
    for i in range(nelx):
        for j in range(nely):
            gon=vtk.vtkQuad()
            gon.GetPointIds().SetId(0,c1); gon.GetPointIds().SetId(1,c2)
            gon.GetPointIds().SetId(2,c3); gon.GetPointIds().SetId(3,c4)
            scl.SetTuple1(c,x[c])
            scl1.SetTuple1(c,s[c])
            cel.InsertNextCell(gon)
            c=c+1; c1=c1+1; c2=c2+1; c3=c3+1; c4=c4+1
        c1=c1+1; c2=c2+1; c3=c3+1; c4=c4+1
            
#
    ply.SetPoints(pts)
    ply.SetPolys(cel)
    ply.SetPoints(pts)
    ply.GetCellData().AddArray(scl)
    ply.GetCellData().AddArray(scl1)
    ply.GetPointData().AddArray(vec)
#
    c=0
    while os.path.isfile('topo_%d.vtp'%c): c=c+1
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('topo_%d.vtp'%c)
    writer.SetInputData(ply)
    writer.Write()
#
