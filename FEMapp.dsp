# Microsoft Developer Studio Project File - Name="FEMapp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=FEMapp - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "FEMapp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "FEMapp.mak" CFG="FEMapp - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "FEMapp - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "FEMapp - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "FEMapp - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MTd /W3 /GX /O2 /I "Algorithms" /I "ControlCenter" /I "Elements" /I "Materials" /I "Matrix" /I "Properties" /I "Output" /I "String" /I "../FemApp" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x408 /d "NDEBUG"
# ADD RSC /l 0x408 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "FEMapp - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "../FemApp/Debug"
# PROP Intermediate_Dir "../FemApp/Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "Algorithms" /I "ControlCenter" /I "Elements" /I "Materials" /I "Matrix" /I "Properties" /I "Output" /I "../FemApp" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x408 /d "_DEBUG"
# ADD RSC /l 0x408 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /pdb:none

!ENDIF 

# Begin Target

# Name "FEMapp - Win32 Release"
# Name "FEMapp - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\FEM.cpp
# End Source File
# Begin Source File

SOURCE=.\GaussPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\InputFileReader.cpp
# End Source File
# Begin Source File

SOURCE=.\Node.cpp
# End Source File
# Begin Source File

SOURCE=.\Statheres.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\FEM.h
# End Source File
# Begin Source File

SOURCE=.\FemIncludes.h
# End Source File
# Begin Source File

SOURCE=.\GaussPoint.h
# End Source File
# Begin Source File

SOURCE=.\Node.h
# End Source File
# Begin Source File

SOURCE=.\StageCommand.h
# End Source File
# Begin Source File

SOURCE=.\Statheres.h
# End Source File
# Begin Source File

SOURCE=.\Output\stopwatch.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Group "Algorithms"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Algorithms\Analysis.cpp
# End Source File
# Begin Source File

SOURCE=.\Algorithms\Analysis.h
# End Source File
# Begin Source File

SOURCE=.\Algorithms\ArcLength.cpp
# End Source File
# Begin Source File

SOURCE=.\Algorithms\ArcLength.h
# End Source File
# Begin Source File

SOURCE=.\Algorithms\ElasticSolver.cpp
# End Source File
# Begin Source File

SOURCE=.\Algorithms\ElasticSolver.h
# End Source File
# End Group
# Begin Group "Elements"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Elements\Beam2D.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\Beam2D.h
# End Source File
# Begin Source File

SOURCE=.\Elements\EightNode2DElasticityElement.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\EightNode2DElasticityElement.h
# End Source File
# Begin Source File

SOURCE=.\Elements\Element.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\Element.h
# End Source File
# Begin Source File

SOURCE=.\Elements\FourNode2DElasticityElement.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\FourNode2DElasticityElement.h
# End Source File
# Begin Source File

SOURCE=.\Elements\IsoparametricElement.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\IsoparametricElement.h
# End Source File
# Begin Source File

SOURCE=.\Elements\Spring2D.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\Spring2D.h
# End Source File
# Begin Source File

SOURCE=.\Elements\SpringXYZ2D.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\SpringXYZ2D.h
# End Source File
# Begin Source File

SOURCE=.\Elements\Truss2D.cpp
# End Source File
# Begin Source File

SOURCE=.\Elements\Truss2D.h
# End Source File
# End Group
# Begin Group "Materials"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Materials\InhomSerialSprings.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\InhomSerialSprings.h
# End Source File
# Begin Source File

SOURCE=.\Materials\Material.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\Material.h
# End Source File
# Begin Source File

SOURCE=.\Materials\PlaneStress.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\PlaneStress.h
# End Source File
# Begin Source File

SOURCE=.\Materials\SpringMaterial.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\SpringMaterial.h
# End Source File
# Begin Source File

SOURCE=.\Materials\SpringMaterialS.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\SpringMaterialS.h
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussMaterial.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussMaterial.h
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussMaterialS.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussMaterialS.h
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussPolyonimial.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\TrussPolyonimial.h
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DModS.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DModS.h
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DS.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DS.h
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DT.cpp
# End Source File
# Begin Source File

SOURCE=.\Materials\Veccio2DT.h
# End Source File
# End Group
# Begin Group "Matrix"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Matrix\BANDMAT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\CHOLESKY.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\EVALUE.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\FFT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\HHOLDER.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\JACOBI.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\MYEXCEPT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWFFT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT1.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT2.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT3.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT4.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT5.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT6.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT7.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT8.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMAT9.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMATEX.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMATNL.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\NEWMATRM.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\SOLUTION.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\SORT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\STR.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\SUBMAT.CPP
# End Source File
# Begin Source File

SOURCE=.\Matrix\SVD.CPP
# End Source File
# End Group
# Begin Group "Properties"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Properties\BeamProperty.cpp
# End Source File
# Begin Source File

SOURCE=.\Properties\BeamProperty.h
# End Source File
# Begin Source File

SOURCE=.\Properties\PlaneElementProperty.cpp
# End Source File
# Begin Source File

SOURCE=.\Properties\PlaneElementProperty.h
# End Source File
# Begin Source File

SOURCE=.\Properties\Property.cpp
# End Source File
# Begin Source File

SOURCE=.\Property.h
# End Source File
# Begin Source File

SOURCE=.\Properties\TrussProperty.cpp
# End Source File
# Begin Source File

SOURCE=.\Properties\TrussProperty.h
# End Source File
# End Group
# Begin Group "Output"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Output\Proccessing.cpp
# End Source File
# Begin Source File

SOURCE=.\Output\Proccessing.h
# End Source File
# End Group
# Begin Group "ControlCenter"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ControlCenter\Constrains.cpp
# End Source File
# Begin Source File

SOURCE=.\ControlCenter\Constrains.h
# End Source File
# Begin Source File

SOURCE=.\ControlCenter\ControlCenter.cpp
# End Source File
# Begin Source File

SOURCE=.\ControlCenter\ControlCenter.h
# End Source File
# Begin Source File

SOURCE=.\ControlCenter\ParametersCenter.cpp
# End Source File
# Begin Source File

SOURCE=.\ControlCenter\ParametersCenter.h
# End Source File
# End Group
# End Target
# End Project
