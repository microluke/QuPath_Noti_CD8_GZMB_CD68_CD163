/**
 * @author Luca Noti
 */

// Import
import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.objects.TMACoreObject
import qupath.lib.roi.interfaces.ROI
import qupath.lib.objects.hierarchy.DefaultTMAGrid
import qupath.lib.roi.RectangleROI
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.Geometry
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathAnnotationObject

// User Input (ImageTypes: "ellipse", "rectangle", null (= none of the above, e.g. full TMA))
String imageType = "ellipse";

// Image data
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

def cal = server.getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
    print ("WARNING: Missing pixel size!")
    return
}

if (!GeneralTools.almostTheSame(cal.getPixelWidthMicrons(), cal.getPixelHeightMicrons(), 0.0001)) {
  print ("WARNING: The pixel width & height are different!")
  return
}
def pixelSize = cal.getAveragedPixelSizeMicrons()

// TMA Preparations
if (imageType == "ellipse" || imageType == "rectangle"){
	int widthPixels = server.getWidth()
	int heightPixels = server.getHeight()
	int diameterNewCore
	PathObject c = new TMACoreObject()
	if (imageType == "ellipse"){
		diameterNewCore = Math.max(widthPixels, heightPixels)
    	c = PathObjects.createTMACoreObject((widthPixels / 2), (heightPixels / 2), diameterNewCore, false)
	}
	if (imageType == "rectangle"){
		ROI roi = new RectangleROI(0, 0, widthPixels, heightPixels)
		c = new TMACoreObject(roi, false)
	}
	
        def tmaGrid = new DefaultTMAGrid([c], 1)
	hierarchy.setTMAGrid(tmaGrid)
	
}

if (!isTMADearrayed()) {
    print ("Please dearray TMA first!")
    return
}

def cores = hierarchy.getTMAGrid().getTMACoreList()

// Defining PathClasses
def tissue = getPathClass("Tissue")
def tumor = getPathClass("Tumor")
def stroma = getPathClass("Stroma")
def stromaZone0 = getPathClass("Stroma: Zone 0")
def stromaZone1 = getPathClass("Stroma: Zone 1")
def stromaZone2 = getPathClass("Stroma: Zone 2")
def stromaZone3 = getPathClass("Stroma: Zone 3")
def cell = getPathClass("Cell")
def other = getPathClass("Other")
def marker1 = getPathClass("Marker 1 positive")
def marker2 = getPathClass("Marker 2 positive")
def marker1and2 = getPathClass("Marker 1 and 2 positive")

// Remove any preexisting objects and TMA measurements
removeObjects(getDetectionObjects(), false)
removeObjects(getAnnotationObjects(), false)
clearTMACoreMeasurements()
print("Removed all objects and measurements. Starting to process " + cores.size() + " TMA cores...")

// Setting default staining vectors
setImageType('FLUORESCENCE');

// Directory for PC folder
String pathProject = buildFilePath(PROJECT_BASE_DIR)
String pathClassifiers = pathProject.substring(0, pathProject.lastIndexOf("\\")) + "/classifiers_project";
String pathPixelClassifiers = pathClassifiers + "/pixel_classifiers/";
String pathObjectClassifiers = pathClassifiers + "/object_classifiers/";

// Create tissue annotation
selectTMACores()
createAnnotationsFromPixelClassifier(pathPixelClassifiers + "PC_Tissue.json", 1000.0, 1000.0)
for (c in cores){
    def annotationTissue = c.getDescendantObjects().findAll{it.isAnnotation() == true && it.getPathClass() == tissue}[0]
    if (annotationTissue != []) {   
        def plane = annotationTissue.getROI().getImagePlane()
        def roiTissue = annotationTissue.getROI()
        def geomTissue = roiTissue.getGeometry()
        def geomTissueShrinked = geomTissue.buffer(-5)
        geomTissueShrinked = GeometryTools.homogenizeGeometryCollection(geomTissueShrinked)
        def roiTissueShrinked = GeometryTools.geometryToROI(geomTissueShrinked, plane)
        def annotationTissueShrinked = PathObjects.createAnnotationObject(roiTissueShrinked)
        annotationTissueShrinked.setPathClass(tissue)
        annotationTissueShrinked.setLocked(true)
        removeObject(annotationTissue, true)
        hierarchy.insertPathObject(annotationTissueShrinked, true)
                         
    }   
}
print("Tissue Areas Annotated and shrinked by 1 pixel!")

// Cell detection in tissue
selectObjects(getAnnotationObjects().findAll{it.getPathClass() == tissue})
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "Channel 1",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 10.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 2.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

// Shape features calculation
selectCells();
addShapeMeasurements("LENGTH", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER")
print("Shape features calculated!")

// Classify Cells vs Artefacts, delete the latters
runObjectClassifier(pathObjectClassifiers + "GrzB_AE1AE3_CD163_OC_CellOther(15_10_IIIB).json");
int numArtefacts = getDetectionObjects().findAll{it.getPathClass() == other}.size();
//removeObjects(getDetectionObjects().findAll{it.getPathClass() == other}, true);
print("Number artefact detections classified: " + numArtefacts);

// Remove measurements (performance...)
removeMeasurements(qupath.lib.objects.PathCellObject, "Nucleus: Length µm", "Nucleus: Solidity", "Nucleus: Max diameter µm", "Nucleus: Min diameter µm", "Cell: Length µm", "Cell: Solidity", "Cell: Max diameter µm", "Cell: Min diameter µm");
print("Measurements removed!")

// Create tumor annotation (DAB)
selectObjects {
   return it.isAnnotation() && it.getPathClass() == tissue
}
createAnnotationsFromPixelClassifier(pathPixelClassifiers + "PC_Tumor_GrzB_AE1AE3_CD163.json", 100.0, 1000.0, "SPLIT")
print("Tumor Areas Annotated!")

// Create stroma annotation
selectObjects {
   return it.isAnnotation() && it.getPathClass() == tissue
}
createAnnotationsFromPixelClassifier(pathPixelClassifiers + "PC_Stroma_GrzB_AE1AE3_CD163.json", 100.0, 0.0)

// Remove tumor annotations without nucleus
removeObjects(getAnnotationObjects().findAll{it.getPathClass() == tissue}, true)
resolveHierarchy()
for (t in getAnnotationObjects().findAll{it.getPathClass() == tumor}){
    if (t.getDescendantObjects().findAll{it.isDetection() == true && it.getPathClass() == cell} == []){
        removeObject(t, true)
    }
}
mergeAnnotations(getAnnotationObjects().findAll{it.getPathClass() == tumor})
annotationTumor = getAnnotationObjects().findAll{it.getPathClass() == tumor}[0]
if (annotationTumor != null) {
    removeObject(annotationTumor, true)
    hierarchy.insertPathObject(annotationTumor, true)
}

// Operations with annotations: Subtract tumor from stroma (different hole settings), crypt necrosis-subtraction from stroma(???), tumor neighborhood zones (10um, 25um, 50um, 100um)
int radiusTN0 = 10
int radiusTN1 = 25
int radiusTN2 = 50
int radiusTN3 = 100
radiusesTN = [radiusTN0, radiusTN1, radiusTN2, radiusTN3]

PathAnnotationObject annotationStromaBase

for (c in cores){
    // Cores with stroma AND tumor annotations
    if ((c.getDescendantObjects().findAll{it.getPathClass() == tumor && it.isAnnotation() == true} != []) && (c.getDescendantObjects().findAll{it.getPathClass() == stroma && it.isAnnotation() == true} != [])) {
        annotationStromaOriginal = c.getDescendantObjects().findAll{it.getPathClass() == stroma && it.isAnnotation() == true}[0]
        annotationTumor = c.getDescendantObjects().findAll{it.getPathClass() == tumor && it.isAnnotation() == true}[0]
       
        def roiStroma = annotationStromaOriginal.getROI()        
        geomStroma = roiStroma.getGeometry()
        removeObject(annotationStromaOriginal, true)
        
        def plane = roiStroma.getImagePlane()
       
        def roiTumor = annotationTumor.getROI()
        geomTumor = roiTumor.getGeometry()
        geomTumor = geomTumor.buffer(-4)
        geomTumor = GeometryTools.homogenizeGeometryCollection(geomTumor)
        roiTumor = GeometryTools.geometryToROI(geomTumor, plane)
        annotationTumorNew = PathObjects.createAnnotationObject(roiTumor)
        annotationTumorNew.setPathClass(tumor)
        annotationTumorNew.setLocked(true)
        removeObject(annotationTumor, true)            
        hierarchy.insertPathObject(annotationTumorNew, true)
        
        geomTumorFilled = GeometryTools.fillHoles(geomTumor)
        
        geomStroma = geomStroma.difference(geomTumorFilled)
        geomStroma = geomStroma.buffer(4)
        geomStromaBase = GeometryTools.homogenizeGeometryCollection(geomStroma)
        
        geomZonesFilled = []
        
        boolean allZonesPresent = true
        
        // Somehow complex architecture, but important for performance...
        for(radiusTN in radiusesTN){
            if (geomZonesFilled == []){
                geomStromaZoneFilled = geomTumorFilled.buffer(radiusTN / pixelSize)
            }
            else{
                geomStromaZoneFilled = geomZonesFilled[radiusesTN.indexOf(radiusTN)-1].buffer((radiusTN - (radiusesTN[radiusesTN.indexOf(radiusTN)-1]))/ pixelSize)
            }
            
            geomStromaZoneFilled = geomStromaZoneFilled.intersection(geomStromaBase)
            geomStromaZoneFilled = GeometryTools.homogenizeGeometryCollection(geomStromaZoneFilled)
            
            if (geomZonesFilled == []) {
                geomZonesFilled.add(geomStromaZoneFilled)
            }
            
            else if (geomStromaZoneFilled != geomZonesFilled[-1]) {
                geomZonesFilled.add(geomStromaZoneFilled)
            }
            
            else {
                allZonesPresent = false
            }
        }
        
        for (geomStromaZone in geomZonesFilled){
            String pathClassZone = "Stroma: Zone " + geomZonesFilled.indexOf(geomStromaZone).toString()
            pathClassZone = getPathClass(pathClassZone)
            print("Computing Annotation for " + pathClassZone + "...")
                         
            if (geomZonesFilled.indexOf(geomStromaZone) == 0){
                geomStromaZone = geomStromaZone.difference(geomTumorFilled)
            }
            else{
                geomStromaZone = geomStromaZone.difference(geomTumorFilled)
                //geomStromaZone = geomStromaZone.difference(geomZonesFilled[geomZonesFilled.indexOf(geomStromaZone)-1])
            }
            
            geomStromaZone = GeometryTools.homogenizeGeometryCollection(geomStromaZone)         
            roiStromaZone = GeometryTools.geometryToROI(geomStromaZone, plane)
            annotationStromaZone = PathObjects.createAnnotationObject(roiStromaZone)
            annotationStromaZone.setPathClass(getPathClass(pathClassZone))
            annotationStromaZone.setLocked(true)

            hierarchy.insertPathObject(annotationStromaZone, true)
            
            if (geomStromaZone == geomStromaBase) {
                allZonesPresent = false
            }

        }

        roiStromaBase = GeometryTools.geometryToROI(geomStromaBase, plane)
        annotationStromaBase = PathObjects.createAnnotationObject(roiStromaBase)
        annotationStromaBase.setPathClass(stroma)
        annotationStromaBase.setLocked(true)
        
        if (allZonesPresent){
            hierarchy.insertPathObject(annotationStromaBase, true)
        }              
    }
    
    // Cores with stroma and without tumor annotations
    else if ((c.getDescendantObjects().findAll{it.getPathClass() == tumor && it.isAnnotation() == true} == []) && (c.getDescendantObjects().findAll{it.getPathClass() == stroma && it.isAnnotation() == true} != [])) {
        annotationStromaBase = c.getDescendantObjects().findAll{it.getPathClass() == stroma && it.isAnnotation() == true}[0]

        removeObject(annotationStromaBase, true)    
        annotationStromaBase.setPathClass(stroma)    
        hierarchy.insertPathObject(annotationStromaBase, true)
        
    }
    
    // Cores without stroma and with tumor annotations
    else if ((c.getDescendantObjects().findAll{it.getPathClass() == tumor && it.isAnnotation() == true} != []) && (c.getDescendantObjects().findAll{it.getPathClass() == stroma && it.isAnnotation() == true} == [])) {
        annotationTumor = c.getDescendantObjects().findAll{it.getPathClass() == tumor && it.isAnnotation() == true}[0]
        
        removeObject(annotationTumor, true)                   
        hierarchy.insertPathObject(annotationTumor, true)      
    }
}
print("Stroma Areas Annotated!")

// Cell classification
selectObjects(getCellObjects())
runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1.0,  "detection[Channel 2]": -1.0,  "detection[Channel 3]": 50.0,  "detection[Channel 4]": -1.0,  "doSmoothing": true,  "splitByIntensity": false,  "splitByShape": false,  "spotSizeMicrons": 1.0,  "minSpotSizeMicrons": 0.5,  "maxSpotSizeMicrons": 20.0,  "includeClusters": false}');

selectObjects(getDetectionObjects().findAll{it.getPathClass() == getPathClass("Subcellular spot: Channel 3 object")})
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 0.25,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": true,  "channel2": true,  "channel3": true,  "channel4": true,  "doMean": true,  "doStdDev": false,  "doMinMax": false,  "doMedian": false,  "doHaralick": false,  "haralickMin": 0,  "haralickMax": 0,  "haralickDistance": 1,  "haralickBins": 32}');

objectsToRemove = []

for (i in getDetectionObjects().findAll{it.getPathClass() == getPathClass("Subcellular spot: Channel 3 object")}){
    def ml = i.getMeasurementList()
    float m1 = ml.getMeasurementValue("ROI: 0.25 µm per pixel: Channel 1: Mean")
    float m2 = ml.getMeasurementValue("ROI: 0.25 µm per pixel: Channel 2: Mean")
    float m3 = ml.getMeasurementValue("ROI: 0.25 µm per pixel: Channel 3: Mean")
    float m4 = ml.getMeasurementValue("ROI: 0.25 µm per pixel: Channel 4: Mean")
    
    boolean m3Specifity = ((m3 > m1*2) && (m3 > m2*2) && (m3 > m4*2))
    
    if (!m3Specifity){
        objectsToRemove.add(i)
    }
}

removeObjects(objectsToRemove, false)
objectsToRemove = []

print("DAB-specific subcellular detections added!")

def m2 = "Cell: Channel 4 mean";
def th2 = 30;

for (i in getDetectionObjects().findAll{it.getPathClass() == cell}){
    // Remove Cells within crypt necrosis (= tumor annotation holes of any size)
    if(i.getLevel() == 2){
        objectsToRemove.add(i)
        continue
    }
    parentClass = i.getParent().getPathClass();
    i.setPathClass(parentClass);
    
    if (i.getDescendantObjects() != [] && measurement(i, m2) >= th2){
            i.setPathClass(marker1and2);
    }
    
    else if (i.getDescendantObjects() != []){
            i.setPathClass(marker1);
    }
    
    else if (measurement(i, m2) >= th2){
            i.setPathClass(marker2);
    }
}
removeObjects(objectsToRemove, false)

// Intensity classification
print("Basic cells and marker-positive cells classified!")

// Calculating values
int tumorAreaMicrons = 0
int tumorNumberBase = 0
int m1TumorPositive = 0
int m2TumorPositive = 0

int stromaTotalAreaMicrons = 0
int stromaTotalNumberBase = 0
int m1StromaTotalPositive = 0
int m2StromaTotalPositive = 0

int stromaZone0AreaMicrons = 0
int stromaZone0NumberBase = 0
int m1StromaZone0Positive = 0
int m2StromaZone0Positive = 0

int stromaZone1AreaMicrons = 0
int stromaZone1NumberBase = 0
int m1StromaZone1Positive = 0
int m2StromaZone1Positive = 0

int stromaZone2AreaMicrons = 0
int stromaZone2NumberBase = 0
int m1StromaZone2Positive = 0
int m2StromaZone2Positive = 0

int stromaZone3AreaMicrons = 0
int stromaZone3NumberBase = 0
int m1StromaZone3Positive = 0
int m2StromaZone3Positive = 0

int stromaZoneTDAreaMicrons = 0
int stromaZoneTDNumberBase = 0
int m1StromaZoneTDPositive = 0
int m2StromaZoneTDPositive = 0


annotationTumor = getAnnotationObjects().findAll{it.getPathClass() == tumor}[0]
if (annotationTumor != null) {
    tumorAreaMicrons = Math.round(annotationTumor.getROI().getArea() * pixelSize * pixelSize)
    tumorNumberBase = getCellObjects().findAll{it.getPathClass() == tumor}.size()
    m1TumorPositive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == tumor}.size()
    m2TumorPositive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == tumor}.size()
}

stromaTotalAreaMicrons = Math.round(annotationStromaBase.getROI().getArea() * pixelSize * pixelSize)

annotationStromaZone0 = getAnnotationObjects().findAll{it.getPathClass() == stromaZone0}[0]
if (annotationStromaZone0 != null) {
    stromaZone0AreaMicrons = Math.round(annotationStromaZone0.getROI().getArea() * pixelSize * pixelSize)
    stromaZone0NumberBase = getCellObjects().findAll{it.getPathClass() != other && it.getParent().getPathClass() == stromaZone0}.size()
    m1StromaZone0Positive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone0}.size()
    m2StromaZone0Positive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone0}.size()
}

annotationStromaZone1 = getAnnotationObjects().findAll{it.getPathClass() == stromaZone1}[0]
if (annotationStromaZone1 != null) {
    stromaZone1AreaMicrons = Math.round(annotationStromaZone1.getROI().getArea() * pixelSize * pixelSize) - stromaZone0AreaMicrons
    stromaZone1NumberBase = getCellObjects().findAll{it.getPathClass() != other && it.getParent().getPathClass() == stromaZone1}.size()
    m1StromaZone1Positive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone1}.size()
    m2StromaZone1Positive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone1}.size()
}

annotationStromaZone2 = getAnnotationObjects().findAll{it.getPathClass() == stromaZone2}[0]
if (annotationStromaZone2 != null) {
    stromaZone2AreaMicrons = Math.round(annotationStromaZone2.getROI().getArea() * pixelSize * pixelSize) - stromaZone0AreaMicrons - stromaZone1AreaMicrons
    stromaZone2NumberBase = getCellObjects().findAll{it.getPathClass() != other && it.getParent().getPathClass() == stromaZone2}.size()
    m1StromaZone2Positive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone2}.size()
    m2StromaZone2Positive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone2}.size()
}

annotationStromaZone3 = getAnnotationObjects().findAll{it.getPathClass() == stromaZone3}[0]
if (annotationStromaZone3 != null) {
    stromaZone3AreaMicrons = Math.round(annotationStromaZone3.getROI().getArea() * pixelSize * pixelSize) - stromaZone0AreaMicrons - stromaZone1AreaMicrons - stromaZone2AreaMicrons
    stromaZone3NumberBase = getCellObjects().findAll{it.getPathClass() != other && it.getParent().getPathClass() == stromaZone3}.size()
    m1StromaZone3Positive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone3}.size()
    m2StromaZone3Positive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stromaZone3}.size()
}

stromaZoneTDAreaMicrons = stromaTotalAreaMicrons - stromaZone0AreaMicrons - stromaZone1AreaMicrons - stromaZone2AreaMicrons - stromaZone3AreaMicrons
stromaZoneTDNumberBase = getCellObjects().findAll{it.getPathClass() != other && it.getParent().getPathClass() == stroma}.size()
m1StromaZoneTDPositive = getCellObjects().findAll{(it.getPathClass() == marker1 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stroma}.size()
m2StromaZoneTDPositive = getCellObjects().findAll{(it.getPathClass() == marker2 || it.getPathClass() == marker1and2) && it.getParent().getPathClass() == stroma}.size()

stromaTotalNumberBase = stromaZone0NumberBase + stromaZone1NumberBase + stromaZone2NumberBase + stromaZone3NumberBase + stromaZoneTDNumberBase
m1StromaTotalPositive = m1StromaZone0Positive + m1StromaZone1Positive + m1StromaZone2Positive + m1StromaZone3Positive + m1StromaZoneTDPositive
m2StromaTotalPositive = m2StromaZone0Positive + m2StromaZone1Positive + m2StromaZone2Positive + m2StromaZone3Positive + m2StromaZoneTDPositive


// Preparing export file
String imageName = getProjectEntry().getImageName();
String imageNameParent = imageName.substring(0, imageName.lastIndexOf("CoreID")-1);
String coreID = imageName.substring(imageName.lastIndexOf("_")+1, imageName.lastIndexOf("."))

String pathExport = buildFilePath(PROJECT_BASE_DIR + "/export")
String pathCSVExportFile = pathExport + "/Summary_" + imageNameParent + ".csv";

File fileCSVExport = new File(pathCSVExportFile)
if (fileCSVExport.exists() == false) {
    print("WARNING: No File for export found.")
    return
}
FileWriter fw = new FileWriter(fileCSVExport, true)
BufferedWriter exportWriter = new BufferedWriter(fw)

// Writing values to export file
StringBuilder sb = new StringBuilder();
sb.append(imageNameParent)
sb.append(";")
sb.append(coreID)
sb.append(";")
sb.append("False")
sb.append(";")

sb.append(tumorAreaMicrons)
sb.append(";")
sb.append(tumorNumberBase)
sb.append(";")
sb.append(m1TumorPositive)
sb.append(";")
sb.append(m2TumorPositive)
sb.append(";")

sb.append(stromaTotalAreaMicrons)
sb.append(";")
sb.append(stromaTotalNumberBase)
sb.append(";")
sb.append(m1StromaTotalPositive)
sb.append(";")
sb.append(m2StromaTotalPositive)
sb.append(";")

sb.append(stromaZone0AreaMicrons)
sb.append(";")
sb.append(stromaZone0NumberBase)
sb.append(";")
sb.append(m1StromaZone0Positive)
sb.append(";")
sb.append(m2StromaZone0Positive)
sb.append(";")

sb.append(stromaZone1AreaMicrons)
sb.append(";")
sb.append(stromaZone1NumberBase)
sb.append(";")
sb.append(m1StromaZone1Positive)
sb.append(";")
sb.append(m2StromaZone1Positive)
sb.append(";")

sb.append(stromaZone2AreaMicrons)
sb.append(";")
sb.append(stromaZone2NumberBase)
sb.append(";")
sb.append(m1StromaZone2Positive)
sb.append(";")
sb.append(m2StromaZone2Positive)
sb.append(";")

sb.append(stromaZone3AreaMicrons)
sb.append(";")
sb.append(stromaZone3NumberBase)
sb.append(";")
sb.append(m1StromaZone3Positive)
sb.append(";")
sb.append(m2StromaZone3Positive)
sb.append(";")

for (i = 0; i < 0; i++){
    sb.append("NaN;")
}

sb.append(stromaZoneTDAreaMicrons)
sb.append(";")
sb.append(stromaZoneTDNumberBase)
sb.append(";")
sb.append(m1StromaZoneTDPositive)
sb.append(";")
sb.append(m2StromaZoneTDPositive)
sb.append("\n")


exportWriter.write(sb.toString());
exportWriter.flush();
exportWriter.close();

print("CSV exported!")

print("Done!")
